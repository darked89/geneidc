/*************************************************************************
*                                                                        *
*   Module: ReadExonsGFF                                                 *
*                                                                        *
*   Reading exons (GFF format) from file                                 *
*                                                                        *
*   This file is part of the geneid 1.4 distribution                     *
*                                                                        *
*     Copyright (C) 2006 - Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           *
*                          Tyler   ALIOTO                                *
*                                                                        *
*  This program is free software; you can redistribute it and/or modify  *
*  it under the terms of the GNU General Public License as published by  *
*  the Free Software Foundation; either version 2 of the License, or     *
*  (at your option) any later version.                                   *
*                                                                        *
*  This program is distributed in the hope that it will be useful,       *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*  GNU General Public License for more details.                          *
*                                                                        *
*  You should have received a copy of the GNU General Public License     *
*  along with this program; if not, write to the Free Software           *
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.             *
*************************************************************************/

/*  $Id: ReadExonsGFF.c,v 1.14 2011-01-13 11:06:16 talioto Exp $  */

#include "geneid.h"
extern float EvidenceEW;
extern float EvidenceFactor;
extern int   FWD;
extern int   RVS;
extern long  LOW;
extern long  HI;

/* According to contig_name, select a group of annotations */
packEvidence *SelectEvidence(packExternalInformation *external,
                             char                    *contig_name){
    int          a;
    packEvidence *p;

    /* 1. Select the position in the array corresponding to contig_name */
    a = getkeyDict(external->locusNames, contig_name);

    if (a == NOTFOUND) {
        p = NULL;
    }
    else {
        p = external->evidence[a];

        /* 2. Init counters */
        external->i1vExons = 0;
        external->i2vExons = 0;
        external->ivExons  = 0;
    }

    return(p);
}

/* Read annotations (exons) to improve or fixed some gene prediction */
/* GFF format: tab "\t" is the field separator and # for comments */
/* Name  Source  Type  Begin  End  Score  Strand  Frame  [group] */
long ReadExonsGFF(char                    *exons_gff_fn,
                  packExternalInformation *external,
                  dict                    *d){

    /* File handle */
    FILE *exons_gff_fptr;

    /* Final number of exons loaded from file (including copies) */
    long i;
    long j;

    /* Split every input line into several tokens (gff records) */
    char line[MAXLINE];
    char lineCopy[MAXLINE];

    char *column_1;
    char *column_2;
    char *column_3;
    char *column_4;
    char *column_5;
    char *column_6;
    char *column_7;
    char *column_8;
    char *column_9;

    /* Control of good sorting property: starting position, increasing */
    long lastAcceptor[MAXNSEQUENCES];
    long currAcceptor;

    /* If frame = '.' then make three copies of current exon (3 frames) */
    int three;

    /* Identifier for a sequence (from dictionary) */
    int     a;

    char    contig_name[CONTIG_NAME_MAX_LENGTH];
    char    saux[MAXTYPE];
    char    c;
    int     slen;
    char    mess[MAXSTRING];

    int     acceptorclass = U2;
    int     donorclass    = U2;
    char    groupCopy[MAXLINE];
    char    *k;
    char    *v;

    char    *utrintrontypes[] = {sUTR5INTRON, sUTR3INTRON};
    int     introncopy;
    int     isIntron          = 0;
    exonGFF *original;

    /* 0. Open exons file to read the information */
    if ((exons_gff_fptr = fopen(exons_gff_fn, "r")) == NULL) {
        printError("The exonsGFF file can not be opened to read");
        exit(EXIT_FAILURE);
    }

    /* 1. Reset counters */
    i     = 0;
    three = 0;
    for (i = 0; i < MAXNSEQUENCES; i++) {
        lastAcceptor[i] = -INFI;
    }

    /* 2. Read while there are exons left in the input file */
    while (fgets(line, MAXLINE, exons_gff_fptr) != NULL) {
        /* 2.a. Comment or empty line: forget it */
        if (line[0] == '#' || line[0] == '\n') {
            printMess("Skipping comment line in evidences file");
        }
        else {
            /* 2.b. Processing exon (annotation) */
            /* Backup copy of the line to show error messages */
            strcpy(lineCopy, line);

            /* Extracting GFF features */
            /* Extracting GFF features */
            column_1 = (char *) strtok(line, "\t");
            column_2 = (char *) strtok(NULL, "\t");
            column_3 = (char *) strtok(NULL, "\t");
            column_4 = (char *) strtok(NULL, "\t");
            column_5 = (char *) strtok(NULL, "\t");
            column_6 = (char *) strtok(NULL, "\t");
            column_7 = (char *) strtok(NULL, "\t");
            column_8 = (char *) strtok(NULL, "\t");
            column_9 = (char *) strtok(NULL, "\n");

            /* There are 8 mandatory columns and the last one is optional */
            if (column_1 == NULL || column_2 == NULL || column_3 == NULL
                || column_4 == NULL || column_5 == NULL || column_6 == NULL
                || column_7 == NULL || column_8 == NULL) {
                sprintf(mess, "Wrong GFF format in HSPs (number of records):\n-->%s\n", lineCopy);
                printError(mess);
                exit(EXIT_FAILURE);
            }

            /* 1. contig_name: leave the exon into the correct array */
            if (sscanf(column_1, "%s", contig_name) != 1) {
                sprintf(mess, "Wrong GFF format in annotations (locusname):\n-->%s\n", lineCopy);
                printError(mess);
                exit(EXIT_FAILURE);
            }

            /* 2. Look-up the ID for that sequence */
            a = setkeyDict(external->locusNames, contig_name);

            if (a >= MAXNSEQUENCES) {
                printError("Too many DNA sequences: increase MAXNSEQUENCES parameter");
                exit(EXIT_FAILURE);
            }

            /* 3. Exon feature: Single, First, Internal, Terminal, ... */
            if (sscanf(column_3, "%s",
                       (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Type) != 1) {
                sprintf(mess, "Wrong GFF format in annotations (feature):\n-->%s\n", lineCopy);
                printError(mess);
                exit(EXIT_FAILURE);
            }

            /* 4. Starting position */
            if (sscanf(column_4, "%ld",
                       &((external->evidence[a]->vSites + external->evidence[a]->nvSites)->Position)) != 1) {
                sprintf(mess, "Wrong GFF format in annotations (starting position):\n-->%s\n", lineCopy);
                printError(mess);
                exit(EXIT_FAILURE);
            }

            /* 5. Finishing position */
            if (sscanf(column_5, "%ld",
                       &((external->evidence[a]->vSites + external->evidence[a]->nvSites + 1)->Position)) != 1) {
                sprintf(mess, "Wrong GFF format in annotations (finishing position):\n-->%s\n", lineCopy);
                printError(mess);
                exit(EXIT_FAILURE);
            }

            /* 6. Score = double value or '.'(infinitum) */
            if (sscanf(column_6, "%f", &((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Score)) != 1) {
                if ((sscanf(column_6, "%c", &c) != 1) || (c != '.')) {
                    sprintf(mess, "Wrong GFF format in annotations (score):\n-->%s\n", lineCopy);
                    printError(mess);
                    exit(EXIT_FAILURE);
                }

                (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Score = MAXSCORE;
            }
            else {
                (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Score = ((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Score * EvidenceFactor) + EvidenceEW;
            }

            /* 7. Strand (reading sense) [+|-] */
            if ((sscanf(column_7, "%c",
                        &((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Strand)) != 1)
                ||
                (((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Strand != '+')
                 &&
                 ((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Strand != '-'))) {
                sprintf(mess, "Wrong GFF format in annotations (strand):\n-->%s\n", lineCopy);
                printError(mess);
                exit(EXIT_FAILURE);
            }

            if (((((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Strand == '-') && RVS)
                 || (((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Strand == '+') && FWD))
                && (
                    (LOW ? ((external->evidence[a]->vSites + external->evidence[a]->nvSites)->Position + 1 >= LOW) : 1)
                    && (HI ? ((external->evidence[a]->vSites + external->evidence[a]->nvSites + 1)->Position < HI) : 1)
                    )) {
                /* 8. Frame = integer or '.' */
                if (sscanf(column_8, "%hd",
                           &((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Frame)) != 1) {
                    /* Is it a dot? */
                    if ((sscanf(column_8, "%c", &c) != 1) || (c != '.')) {
                        sprintf(mess, "Wrong GFF format in annotations (frame):\n-->%s\n", lineCopy);
                        printError(mess);
                        exit(EXIT_FAILURE);
                    }

                    /* make three copies */
                    three = 1;
                }
                else {
                    /* Checking input frame between 0..2 */
                    if (((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Frame < 0)
                        ||
                        ((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Frame > 2)) {
                        sprintf(mess, "Wrong GFF value in annotations (frame not between 0..2):\n-->%s\n", lineCopy);
                        printError(mess);
                        exit(EXIT_FAILURE);
                    }
                }

                /* 9. Group: optional, string */
                if (column_9 != NULL) {
                    if (sscanf(column_9, "%s",
                               (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Group) != 1) {
                        sprintf(mess, "Wrong GFF value in annotations (group):\n-->%s\n", lineCopy);
                        printError(mess);
                        exit(EXIT_FAILURE);
                    }

                    if (!(strcmp((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Group, "."))) {
                        /* strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Group, ""); */
                        strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Group, NOGROUP);
                    }
                    else {
                        /* Handle GFF3 format here*/
                        /* Parse group field into tokens delimited by '=;' */
                        /* copy the field first so as not to screw it up */
                        strcpy(groupCopy, column_9);
                        k = (char *) strtok(groupCopy, "=");

                        if (k != NULL) {
                            strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Group, NOGROUP);
                        }

                        while (k != NULL) {

                            v = (char *) strtok(NULL, ";\n");

/*                    sprintf(mess,"Tags\n%s-->%s\n",k,v); */
/*                    printMess(mess); */
                            if (!(strcmp(k, "Parent"))) {
                                if (!(strcmp(v, ".")) || v == NULL) {
                                    strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Group, NOGROUP);
                                }
                                else {
                                    strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Group, v);
                                }
                            }

                            if (
                                (
                                    (
                                        (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Strand == '+')
                                    && strcmp((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Type, sINTRON)
                                    && strcmp((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Type, sUTRINTRON)
                                )
                                ||
                                (
                                    (
                                        (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Strand == '-')
                                    && (!strcmp((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Type, sINTRON)
                                        || !strcmp((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Type, sUTRINTRON)
                                        )
                                )
                                ) {
                                if (!(strcmp(k, "donor"))) {
                                    if (!(strcmp(v, sU2))) {
                                        donorclass = U2;
                                    }

                                    if (!(strcmp(v, sU12gtag))) {
                                        donorclass = U12gtag;
                                    }

                                    if (!(strcmp(v, sU12atac))) {
                                        donorclass = U12atac;
                                    }
                                }

                                if (!(strcmp(k, "acceptor"))) {
                                    if (!(strcmp(v, sU2))) {
                                        acceptorclass = U2;
                                    }

                                    if (!(strcmp(v, sU12gtag))) {
                                        acceptorclass = U12gtag;
                                    }

                                    if (!(strcmp(v, sU12atac))) {
                                        acceptorclass = U12atac;
                                    }

                                }
                            }
                            else {
                                if (!(strcmp(k, "acceptor"))) {
                                    if (!(strcmp(v, sU2))) {
                                        donorclass = U2;
                                    }

                                    if (!(strcmp(v, sU12gtag))) {
                                        donorclass = U12gtag;
                                    }

                                    if (!(strcmp(v, sU12atac))) {
                                        donorclass = U12atac;
                                    }
                                }

                                if (!(strcmp(k, "donor"))) {
                                    if (!(strcmp(v, sU2))) {
                                        acceptorclass = U2;
                                    }

                                    if (!(strcmp(v, sU12gtag))) {
                                        acceptorclass = U12gtag;
                                    }

                                    if (!(strcmp(v, sU12atac))) {
                                        acceptorclass = U12atac;
                                    }

                                }
                            }

/*                    sprintf(mess,"Tags\n%s-->%s\n",k,v); */
/*                    printMess(mess); */

                            k = strtok(NULL, "=");

                        }
/*                sprintf(mess,"donorclass=%i\tacceptorclass=%i\n",donorclass,acceptorclass); */
/*                    printMess(mess); */

                        /* if only one token, then leave it -- it's the group in GFF1 format */
                        /* process them in pairs, looking for 'Parent', 'donor', 'acceptor' */
                        /* valid values for donor and acceptor are 'U2', 'U12gtag', and 'U12atac' */
                        /* later we can add tags for stop codon values to ensure non-creation of stops across splice junctions */
                    }
                }
                else {
                    /* This exon will be allowed to join to ab initio predictions */
                    strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Group, NOGROUP);
                }

                /* 2.c. Process current exon */
                /* (A). Checking exon feature (gene model): type.strand */
                saux[0]      = '\0';
                strcpy(saux, (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Type);
                slen         = strlen(saux);
                saux[slen++] = (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Strand;
                saux[slen]   = '\0';

                if (getkeyDict(d, saux) == NOTFOUND) {
                    /* Forget it: exon type is not in current gene model */
                    sprintf(mess, "Wrong GFF feature in annotations (unknown):\n-->%s\n", lineCopy);
                    printMess(mess);
                }
                else {
                    /* (B). Well-sorted (by start position) list of read exons */
                    currAcceptor = (external->evidence[a]->vSites + external->evidence[a]->nvSites)->Position;

                    if (lastAcceptor[a] > currAcceptor) {
                        sprintf(mess, "Order violation: annotations (starting position %ld):\n-->%s\n",
                                lastAcceptor[a],
                                lineCopy);
                        printError(mess);
                        exit(EXIT_FAILURE);
                    }

                    else {
                        lastAcceptor[a] = currAcceptor;

                        if (!strcmp((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Type, sINTRON)) {
                            isIntron = 1;
                        }
                        else {
                            isIntron = 0;
                        }

                        original                                                                    = external->evidence[a]->vExons + external->evidence[a]->nvExons;
                        /* (C). Setting evidence splice sites to U2 class */
                        (external->evidence[a]->vSites + external->evidence[a]->nvSites)->class     = acceptorclass;
                        (external->evidence[a]->vSites + external->evidence[a]->nvSites + 1)->class = donorclass;
/*                  sprintf(mess,"donorclass=%i\tacceptorclass=%i\n",(external->evidence[a]->vSites + external->evidence[a]->nvSites + 1)->class,(external->evidence[a]->vSites + external->evidence[a]->nvSites)->class); */
/*                    printMess(mess); */
                        /* (C). Setting dummy sites to this exon */
                        (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Acceptor
                            = (external->evidence[a]->vSites + external->evidence[a]->nvSites);
                        (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Donor
                            = (external->evidence[a]->vSites + external->evidence[a]->nvSites + 1);

                        /* Updating information about sites to the range 0..L-1 */
                        (external->evidence[a]->vExons + external->evidence[a]->nvExons)->offset1 = -COFFSET;
                        (external->evidence[a]->vExons + external->evidence[a]->nvExons)->offset2 = -COFFSET;

                        /* (D). Making three (two more) copies if needed */
                        if (three) {
                            /* Creating three exons (3 frames): sharing sites */
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Acceptor
                                = (external->evidence[a]->vSites + external->evidence[a]->nvSites);
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Donor
                                = (external->evidence[a]->vSites + external->evidence[a]->nvSites + 1);
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Acceptor
                                = (external->evidence[a]->vSites + external->evidence[a]->nvSites);
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Donor
                                = (external->evidence[a]->vSites + external->evidence[a]->nvSites + 1);

                            /* Updating information about sites to the range 0..L-1 */
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->offset1 = -COFFSET;
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->offset2 = -COFFSET;
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->offset1 = -COFFSET;
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->offset2 = -COFFSET;

                            /* Setting frame values */
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Frame     = 0;
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Frame = 1;
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Frame = 2;

                            /* Copy some exon attributes */
                            strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Type,
                                   (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Type);
                            strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Type,
                                   (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Type);
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Score  = original->Score;
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Score  = original->Score;
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Strand = original->Strand;
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Strand = original->Strand;

                            /* Computing remainder from frame value for copies 1,2 */
                            if (!strcmp((external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Type, sINTRON)) {
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Remainder = (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Frame;
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Frame
                                                                                                                = (3 - (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Remainder) % 3;
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Remainder = (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Frame;
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Frame
                                                                                                                = (3 - (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Remainder) % 3;
                            }
                            else {
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Remainder
                                    = ((3 - (((external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Donor->Position
                                              - (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Acceptor->Position
                                              - (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Frame + 1) % 3)) % 3);

                                (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Remainder
                                    = ((3 - (((external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Donor->Position
                                              - (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Acceptor->Position
                                              - (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Frame + 1) % 3)) % 3);
                            }

                            /* The same group */
                            strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Group,
                                   original->Group);
                            strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Group,
                                   original->Group);

                            /* Evidence flag activated */
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->evidence = 1;
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->evidence = 1;

                        } /* End of if(three) */

                        /* (E). Doing the same for the original exon: */

                        /* Computing remainder from frame value */
                        if (!strcmp((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Type, sINTRON)) {
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Remainder = (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Frame;
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Frame
                                                                                                        = (3 - (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Remainder) % 3;
                        }
                        else {
                            (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Remainder
                                = ((3 - (((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Donor->Position
                                          - (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Acceptor->Position
                                          - (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Frame + 1) % 3)) % 3);
                        }

                        /* Evidence flag activated */
                        (external->evidence[a]->vExons + external->evidence[a]->nvExons)->evidence = 1;

                        /* Updating and checking loop counters */
                        external->evidence[a]->nvExons = (three) ?
                                                         external->evidence[a]->nvExons + 3 :
                                                         external->evidence[a]->nvExons + 1;

                        /* (D). Making necessary intron types */
                        if (isIntron) {
                            /* printMess("Making UTR3 and UTR5 Introns"); */
                            for (introncopy = 0; introncopy < 2; introncopy++) {
                                if (three) {
                                    /* Creating three exons (3 frames): sharing sites */
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Acceptor
                                        = (external->evidence[a]->vSites + external->evidence[a]->nvSites);
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Donor
                                        = (external->evidence[a]->vSites + external->evidence[a]->nvSites + 1);
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Acceptor
                                        = (external->evidence[a]->vSites + external->evidence[a]->nvSites);
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Donor
                                        = (external->evidence[a]->vSites + external->evidence[a]->nvSites + 1);

                                    /* Updating information about sites to the range 0..L-1 */
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->offset1 = -COFFSET;
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->offset2 = -COFFSET;
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->offset1 = -COFFSET;
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->offset2 = -COFFSET;

                                    /* Setting frame values */
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Frame     = 0;
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Frame = 1;
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Frame = 2;

                                    /* Copy some exon attributes */
                                    strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Type, utrintrontypes[introncopy]);
                                    strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Type, utrintrontypes[introncopy]);
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Score  = original->Score;
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Score  = original->Score;
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Strand = original->Strand;
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Strand = original->Strand;

                                    /* Computing remainder from frame value for copies 1,2 */
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Remainder = (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Frame;
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Frame
                                                                                                                    = (3 - (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Remainder) % 3;
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Remainder = (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Frame;
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Frame
                                                                                                                    = (3 - (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Remainder) % 3;

                                    /* The same group */
                                    strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->Group, original->Group);
                                    strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->Group, original->Group);

                                    /* Evidence flag activated */
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 1)->evidence = 1;
                                    (external->evidence[a]->vExons + external->evidence[a]->nvExons + 2)->evidence = 1;
                                }

                                /* (E). Doing the same for the original exon: */
                                /* Creating copy exon: sharing sites */
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Acceptor
                                    = (external->evidence[a]->vSites + external->evidence[a]->nvSites);
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Donor
                                    = (external->evidence[a]->vSites + external->evidence[a]->nvSites + 1);

                                /* Updating information about sites to the range 0..L-1 */
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons)->offset1 = -COFFSET;
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons)->offset2 = -COFFSET;

                                /* Setting frame values */
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Frame = 0;

                                /* Copy some exon attributes */
                                strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Type, utrintrontypes[introncopy]);
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Score  = original->Score;
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Strand = original->Strand;

                                /* The same group */
                                strcpy((external->evidence[a]->vExons + external->evidence[a]->nvExons)->Group, original->Group);

                                /* Computing remainder from frame value */
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Remainder = (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Frame;
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Frame
                                                                                                            = (3 - (external->evidence[a]->vExons + external->evidence[a]->nvExons)->Remainder) % 3;

                                /* Evidence flag activated */
                                (external->evidence[a]->vExons + external->evidence[a]->nvExons)->evidence = 1;
                                /* Updating and checking loop counters */
                                external->evidence[a]->nvExons                                             = (three) ?
                                                                                                             external->evidence[a]->nvExons + 3 :
                                                                                                             external->evidence[a]->nvExons + 1;
                            }
                        } /* End of if(Intron) */

                        external->evidence[a]->nvSites = external->evidence[a]->nvSites + 2;
                        three                          = 0;

                        if ((i + FRAMES) > MAXEVIDENCES) {
                            printError("Too many annotations: increase MAXEVIDENCES definition");
                            exit(EXIT_FAILURE);
                        }

                        if ((external->evidence[a]->nvSites + (2 * FRAMES)) > MAXSITESEVIDENCES) {
                            printError("Too many site annotations: increase MAXEVIDENCES definition");
                            exit(EXIT_FAILURE);
                        }

                    } /* End of sorting checkpoint */

                } /* End of feature_in_gene_model checkpoint */

            }/* End of Strand checkpoint */

        } /* End of checkpoint for comments (#) */

    } /* End of while (input exons) */
    fclose(exons_gff_fptr);

    /* Obtain the number of different sequences */
    external->nSequences = external->locusNames->nextFree;

    /* Return the number of created exons (including replications) */
    for (i = 0, j = 0; j < external->nSequences; j++) {
/*     sprintf(mess,"nvExons: %ld",external->evidence[j]->nvExons);printMess(mess); */
        i = i + external->evidence[j]->nvExons;
    }

    return(i);
}
