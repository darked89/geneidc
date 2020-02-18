/*************************************************************************
*                                                                        *
*   Module: geneid                                                       *
*                                                                        *
*   geneid main program                                                  *
*                                                                        *
*   This file is part of the geneid 1.4 distribution                     *
*                                                                        *
*     Copyright (C) 2006 - Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           *
*                          Tyler   ALIOTO                                *
*     with contributions from:                                           *
*                          Moises  BURSET ALVAREDA                       *
*                          Genis   PARRA FARRE                           *
*                          Xavier  MESSEGUER PEYPOCH                     *
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

/* $Id: geneid.c,v 1.27 2011-01-13 11:06:16 talioto Exp $ */

#include <assert.h>
#include "geneid.h"
/* #include <mcheck.h> */

/* geneid setup flags */
/* sites to print */
int    SFP      = 0;
int    SDP      = 0;
int    SAP      = 0;
int    STP      = 0;
/* exons to print */
int    EFP      = 0;
int    EIP      = 0;
int    ETP      = 0;
int    EXP      = 0;
int    ESP      = 0;
int    EOP      = 0;
/* introns to print */
int    PRINTINT = 0;
/* Partial or full prediction engine */
int    GENAMIC  = 1;
int    GENEID   = 1;
/* Only forward or reverse prediction engine */
int    FWD      = 1;
int    RVS      = 1;
/* switch ORF prediction on */
int    scanORF  = 0;
/* Input annotations or homology to protein information/reads to UTR prediction */
int    EVD      = 0;
int    SRP      = 0;
int    UTR      = 0;
/* Output formats */
int    GFF      = 0;
int    GFF3     = 0;
int    X10      = 0;
int    XML      = 0;
int    cDNA     = 0;
int    PSEQ     = 0;
int    tDNA     = 0;
/* Verbose flag (memory/processing information) */
int    BEG      = 0;
int    VRB      = 0;
/* Score for regions not-supported by protein homology */
int   NO_SCORE;
/* Force single prediction: 1 gene */
int    SGE                        = 0;
/* Detection of PolyPTracts in Acceptors */
int    PPT                        = 0;
/* Detection of BranchPoints in Acceptors */
int    BP                         = 0;
/* Detection of recursive splice sites */
int    RSS                        = 0;
/* Detection of U12 introns */
int U12                        = 0;
/* Detection of U12gtag sites (acceptor uses BranchPoint)*/
int U12GTAG                    = 0;
/* Detection of U12atac sites (acceptor uses BranchPoint)*/
int U12ATAC                    = 0;
/* Detection of U2gcag sites */
int U2GCAG                     = 0;
/* Detection of U2gta donor sites */
int U2GTA                      = 0;
/* Detection of U2gtg donor sites */
int U2GTG                      = 0;
/* Detection of U2gty donor sites */
int U2GTY                      = 0;
/* Detection of PolyA Signal */
int PAS                        = 0;
/* Length of flank around exons to subtract background RNA-seq signal */
int BKGD_SUBTRACT_FLANK_LENGTH = 0;

/* Splice classes: the number of compatible splice site combinations used in genamic for joining exons */
unsigned short SPLICECLASSES = 1;

/* User defined lower limit */
long LOW = 0;
/* User defined upper limit */
long HI  = 0;

/* Millions of reads mapped */
float MRM = 15.0;

/* Optional Predicted Gene Prefix */
char GenePrefix[MAXSTRING] = "";

/* Increase/decrease exon weight value (exon score) */
float EW                      = NOVALUE;
float U12EW                   = 0;
float EvidenceEW              = 0;
float EvidenceFactor          = 1;
float U12_SPLICE_SCORE_THRESH = -1000;
float U12_EXON_SCORE_THRESH   = -1000;

/* Detection of recursive splice sites */
float RSSMARKOVSCORE = 0;
float RSSDON         = RDT;
float RSSACC         = RAT;

/* Generic maximum values: sites, exons and backup elements */
long NUMSITES;
long NUMEXONS;
long MAXBACKUPSITES;
long MAXBACKUPEXONS;
long NUMU12SITES;
long NUMU12EXONS;
long NUMU12U12EXONS;

/* Accounting time and results */
account *m;

/************************************************************************
                            geneid MAIN program
************************************************************************/

int main(int    argc,
         char  *argv[]){

    /* DNA sequence data structures */
    FILE       *fasta_fptr;
    char       *Sequence;
    char       *RSequence;
    long        LengthSequence;

    /* Current split ends */
    long        l1;
    long        l2;
    long        upperlimit;
    long        lowerlimit;

    /* Forward sense data structures */
    packSites  *allSites;
    packExons  *allExons;

    /* Reverse sense data structures */
    packSites  *allSites_r;
    packExons  *allExons_r;

    /* Structures for sorting sites */
    site       *donorsites;
    site       *acceptorsites;
    site       *tssites = NULL;
    site       *tesites = NULL;

    /* Table to sort predicted exons by acceptor */
    exonGFF    *exons;
    long        nExons;

    /* External information: reannotation */
    packExternalInformation  *external;
    packEvidence             *evidence;
    packHSP                  *hsp;

    /* Best partial predicted genes */
    packGenes  *genes;

    /* Dumpster for backup operations between splits */
    packDump   *dumpster;

    /* Amino acid dictionary (genetic code) */
    dict       *dAA;

    /* geneid prediction parameters: data structures */
    gparam     *gp = NULL;
    gparam    **isochores;

    /* Input Filenames */
    char        fasta_fn[FILENAMELENGTH]     = "";
    char        param_fn[FILENAMELENGTH]     = "";
    char        exons_gff_fn[FILENAMELENGTH]  = "";
    char        blastHSP_gff_fn[FILENAMELENGTH]      = "";

    /* contig_name sequence name */
    char        contig_name[CONTIG_NAME_MAX_LENGTH];
    char        next_contig_name[CONTIG_NAME_MAX_LENGTH];

    /* Measure of C+G content to select the isochore */
    packGC     *GCInfo;
    packGC     *GCInfo_r;
    int         inigc;
    int         endgc;
    float        percentGC;
    int         currentIsochore;
    int         nIsochores;
    int         reading;
    int         lastSplit;
    char        mess[MAXSTRING];

    /* Start memory trace -- for debugging memory leaks */
    /* mtrace(); */

    /** 0. Starting and reading options, parameters and sequence... **/
    nExons   = 0;
    evidence = NULL;
    hsp      = NULL;

    /* 0.a. Previous checkpoint about length in splits and overlapping */
    if (LENGTHSi <= OVERLAP) {
        printError("LENGTHSi must be greater than OVERLAP parameter (geneid.h)");
        exit(EXIT_FAILURE);
    }

    /* 0.b. Initializing stats and time counters */
    m = (account *) InitAcc();

    /* 0.c. Read setup options */
    readargv(argc,
             argv,
             param_fn,
             fasta_fn,
             exons_gff_fn,
             blastHSP_gff_fn,
             GenePrefix);

    printRes("\n\n\t\t\t** Running geneid 1.4.5+ 2020 geneid@crg.es **\n\n");

    /* 0.d. Prediction of DNA sequence length to request memory */
    LengthSequence = get_fasta_size(fasta_fn);
    sprintf(mess, "DNA sequence file size = %ld bytes", LengthSequence);
    printMess(mess);

    /* 0.e. Computing ratios for every type of signal and exons */
    printMess("Computing Ratios");
    SetRatios(&NUMSITES,
              &NUMEXONS,
              &MAXBACKUPSITES,
              &MAXBACKUPEXONS,
              LengthSequence);

    /* Estimation of memory required to execute geneid */
    if (BEG) {
        beggar(LengthSequence);
    }

    /** 1. Allocating main geneid data structures **/
    printMess("Request Memory to Operating System\n");

    /* 1.a. Mandatory geneid data structures */
    printMess("Request Memory Sequence\n");
    Sequence   = (char *) RequestMemorySequence(LengthSequence);
    RSequence  = (char *) RequestMemorySequence(LengthSequence);

    printMess("Request Memory Sites\n");
    allSites   = (packSites *) RequestMemorySites();
    allSites_r = (packSites *) RequestMemorySites();

    printMess("Request Memory Exons\n");
    allExons   = (packExons *) RequestMemoryExons();
    allExons_r = (packExons *) RequestMemoryExons();

    printMess("Request Memory Sort Exons\n");

    exons         = (exonGFF *) RequestMemorySortExons();

    printMess("Request Memory Sort Sites\n");
    donorsites    = (site *) RequestMemorySortSites();      /* Temporary structure for sorting donor sites */
    acceptorsites = (site *) RequestMemorySortSites();      /* Temporary structure for sorting acceptor sites */

    if (UTR) {
        tssites = (site *) RequestMemorySortSites();    /* Temporary structure for sorting acceptor sites */
        tesites = (site *) RequestMemorySortSites();    /* Temporary structure for sorting acceptor sites */
    }

    printMess("Request Memory Isochores, etc.\n");
    isochores = (gparam **) RequestMemoryIsochoresParams();
    GCInfo    = (packGC *) RequestMemoryGC();
    GCInfo_r  = (packGC *) RequestMemoryGC();
    dAA       = (dict *) RequestMemoryAaDictionary();
    printMess("Request Memory External\n");
    external  = (packExternalInformation *) RequestMemoryExternalInformation();

    /* 1.b. Backup information might be necessary between splits */
    if (LengthSequence > LENGTHSi) {
        printMess("Request Memory Dumpster\n");
        dumpster = (packDump *) RequestMemoryDumpster();
    }
    else {
        dumpster = NULL;
    }

    /** 2. Reading statistical model parameters file **/
    printMess("Reading parameters...");
    nIsochores = readparam(param_fn, isochores);

    if (U12) {
        if ((!U12GTAG) && (!U12ATAC)) {
            U12GTAG = 0;
            U12ATAC = 0;
        }
    }
    else {
        U12GTAG = 0;
        U12ATAC = 0;
    }

    /** 1. Allocating genes data structure (after the number of splice classes has been determined **/
    genes = (packGenes *) RequestMemoryGenes();

    /** 3. Starting processing: complete or partial prediction ? **/
    if (GENEID) {
        /* A. Predicting signals, exons and genes in DNA sequences */
        /* A.1. Reading external information I: annotations */
        if (EVD) {
            printMess("Reading evidence (annotations)...");
            external->nvExons
                = ReadExonsGFF(exons_gff_fn, external, isochores[0]->D);
            sprintf(mess, "%ld annotations acquired from file\n",
                    external->nvExons);
            printMess(mess);
        }

        /* A.2. Reading external information II: homology information */
        if (SRP) {
            printMess("Reading homology information...");
            external->nHSPs = ReadHSP(blastHSP_gff_fn, external);
            sprintf(mess, "%ld HSPs acquired from file",
                    external->nHSPs);
            printMess(mess);
        }

        if (EVD || SRP) {
            sprintf(mess, "External information acquired from %ld sequences\n",
                    external->nSequences);
            printMess(mess);
        }

        /** A.3. Input DNA sequences (perhaps more than one) **/
        if ((fasta_fptr = fopen(fasta_fn, "rb")) == NULL) {
            printError("The input sequence file can not be accessed");
        }

        /* reading the locusname of sequence (in Fasta format) */
        reading = IniReadSequence(fasta_fptr, contig_name);

        while (reading != EOF) {
            printMess("Loading DNA sequence");
            reading = ReadSequence(fasta_fptr, Sequence, next_contig_name);

            /* A.3. Prepare sequence to work on */
            printMess("Processing DNA sequence");
            LengthSequence = FetchSequence(Sequence, RSequence);
            OutputHeader(contig_name, LengthSequence);

            /* A.4. Prepare external information */
            if (SRP) {
                printMess("Select homology information");
                hsp = (packHSP *) SelectHSP(external, contig_name, LengthSequence);

                if (hsp == NULL) {
                    sprintf(mess, "No information has been provided for %s\n",
                            contig_name);
                }
                else {
                    sprintf(mess, "Using %ld HSPs in %s\n",
                            hsp->nTotalSegments, contig_name);
                }

                printMess(mess);
            }

            if (EVD) {
                printMess("Select annotations");
                evidence = (packEvidence *) SelectEvidence(external, contig_name);

                if (evidence == NULL) {
                    sprintf(mess, "No information has been provided for %s\n",
                            contig_name);
                }
                else {
                    sprintf(mess, "Using %ld annotations in %s\n",
                            evidence->nvExons, contig_name);
                }

                printMess(mess);
            }

            /* A.5. Processing sequence into several fragments if required */
            /* l1 is the left end and l2 is the right end in Sequence */
            /* The arguments HI and LO are converted into lower and upper limit coordinates and l1 and l2 are adjusted */
            upperlimit = LengthSequence - 1;
            lowerlimit = 0;

            if ((HI > 0) && (HI >= LOW) && (HI < LengthSequence)) {
                upperlimit = HI - 1;
            }
            else {
                upperlimit = LengthSequence - 1;
            }

            if ((LOW > 0) && (LOW <= upperlimit)) {
                lowerlimit = LOW - 1;
            }
            else {
                lowerlimit = 0;
            }

            l1        = lowerlimit;
            l2        = MIN(l1 + LENGTHSi - 1, LengthSequence - 1);
            l2        = MIN(l2, upperlimit);
            /* Check to see if we are on last split */
            lastSplit = (l2 == upperlimit);
            sprintf(mess, "Running on range %ld to %ld\n",
                    lowerlimit, upperlimit);
            printMess(mess);

            while ((l1 < (upperlimit + 1 - OVERLAP)) || (l1 == 0) || (l1 == lowerlimit)) {
                /** B.1. Measure G+C content in the current fragment: l1,l2 **/
                GCScan(Sequence, GCInfo, l1, l2);
                GCScan(RSequence, GCInfo_r, LengthSequence - 1 - l2, LengthSequence - 1 - l1);

                /* G+C range: from 0 (l1) to l2 -l1 (l2) */
                /*
                 * inigc     = l1 - l1; //error, 20200218
                 */
                inigc     = 0;
                endgc     = l2 - l1;
                percentGC = ComputeGC(GCInfo, inigc, endgc);
                sprintf(mess, "G+C content in [%ld-%ld] is %f", l1, l2, percentGC);
                printMess(mess);

                /* Choose the isochore to predict sites according to the GC level */
                currentIsochore = SelectIsochore(percentGC, isochores);
                gp              = isochores[currentIsochore];
                sprintf(mess, "Selecting isochore %d", currentIsochore + COFFSET);
                printMess(mess);

                /* B.2. Prediction of sites and exons construction/filtering */
                if (FWD) {
                    /* Forward strand predictions */
                    sprintf(mess, "Running FWD  %s: %ld - %ld", contig_name, l1, l2);
                    printMess(mess);
                    manager(Sequence,
                            LengthSequence,
                            allSites,
                            allExons,
                            l1,
                            l2,
                            lowerlimit,
                            upperlimit,
                            FORWARD,
                            external,
                            hsp,
                            gp,
                            isochores,
                            nIsochores,
                            GCInfo,
                            acceptorsites,
                            donorsites,
                            tssites,
                            tesites);
                }

                if (RVS) {
                    /* Reverse strand predictions */
                    sprintf(mess,
                            "Running Reverse  %s: %ld - %ld(%ld - %ld)",
                            contig_name, 
                            LengthSequence - 1 - l2,
                            LengthSequence - 1 - l1, l1, l2);
                    printMess(mess);

                    manager(RSequence, LengthSequence,
                            allSites_r, allExons_r,
                            LengthSequence - 1 - l2,
                            LengthSequence - 1 - l1,
                            LengthSequence - 1 - upperlimit, LengthSequence - 1 - lowerlimit,
                            REVERSE,
                            external, hsp, gp,
                            isochores, nIsochores,
                            GCInfo_r, acceptorsites, donorsites, tssites, tesites);

                    /* normalised positions: according to forward sense reading */
                    RecomputePositions(allSites_r, LengthSequence);

                    /* exchange acc and donor sites to preserve Acc < Don */
                    SwitchPositions(allExons_r);
                }

                /* B.3. Sort all of exons by left (minor) position */
                if (EVD && evidence != NULL) {
                    /* Searching evidence exons in this fragment */
                    printMess("Searching annotations to be used in this fragment");
                    SearchEvidenceExons(external,
                                        evidence,
                                        (lastSplit) ? l2 : l2 - OVERLAP);

                    /* Unused annotations: out of range (info) */
                    if (lastSplit) {
                        sprintf(mess, "Leaving out last %ld evidences (out of range)",
                                evidence->nvExons - external->i2vExons);
                        printMess(mess);
                    }
                }

                nExons = allExons->nExons + allExons_r->nExons;

                if (EVD && evidence != NULL) {
                    nExons = nExons + external->ivExons;
                }

                /* BEGIN artificial exon: + and - */
                if (l1 == lowerlimit) {
                    nExons = nExons + 2;
                }

                /* END artitificial exon: + and - */
                if (l2 == upperlimit) {
                    nExons = nExons + 2;
                }

/*            sprintf(mess,"l1: %ld   ll:%ld   l2: %ld   ul: %ld\n", l1,lowerlimit,l2,upperlimit); */
/*            printMess(mess); */
/*            /\* B.4. Printing current fragment predictions (sites and exons) *\/ */
/*            Output(allSites, allSites_r, allExons, allExons_r,  */
/*                   exons, nExons, contig_name, l1, l2, lowerlimit, Sequence, gp, dAA, GenePrefix);  */

                sprintf(mess, "Sorting %ld exons\n", nExons);
                printMess(mess);

                /* Merge predicted exons with some evidence exons */
                SortExons(allExons, allExons_r,
                          external,
                          evidence,
                          exons,
                          l1,
                          l2,
                          lowerlimit,
                          upperlimit);
                sprintf(mess, "Finished sorting %ld exons\n", nExons);
                printMess(mess);

                /* Next block of annotations to be processed */
                if (EVD && evidence != NULL) {
                    SwitchCounters(external);
                }

                /* B.4. Printing current fragment predictions (sites and exons) */
                Output(allSites, allSites_r, allExons, allExons_r,
                       exons, nExons, contig_name, l1, l2, lowerlimit, Sequence, gp, dAA, GenePrefix);

                /* recompute stats about splice sites and exons */
                updateTotals(m, allSites, allSites_r, allExons, allExons_r);

                /* B.5. Calling to genamic for assembling the best gene */
                if (GENAMIC && nExons) {

                    genamic(exons, nExons, genes, gp);

                    if (upperlimit - lowerlimit + 1 > LENGTHSi) {/*  if (LengthSequence > LENGTHSi) */
                        /* clean hash table of exons */
                        cleanDumpHash(dumpster->h);
                    }

                    /* B.6. Backup operations of genes for the next split */
                    if (!lastSplit) {
                        /* backup of unused genes */
                        printMess("Back-up of d-genes");
                        BackupArrayD(genes, l2 - OVERLAP, gp, dumpster);

                        /* back-up best partial genes */
                        printMess("Back-up of best partial genes\n");
                        BackupGenes(genes, gp->nclass, dumpster);
                    }
                }

                /* Computing new boundaries: next fragment in current sequence */
                l1       += LENGTHSi - OVERLAP;
                l2        = MIN(l1 + LENGTHSi - 1, upperlimit);
                lastSplit = (l2 == upperlimit);
            } /* processing next fragment */

            /* A.6. Full sequence processed: displaying best predicted gene */
            if (GENAMIC) {
                /* Printing gene predictions */
                OutputGene(genes,
                           (EVD && evidence != NULL) ?
                           m->totalExons + evidence->nvExons :
                           m->totalExons,
                           contig_name, 
                           Sequence, 
                           gp, 
                           dAA, 
                           GenePrefix);

                /* Reset best genes data structures for next input sequence */
                printMess("Cleaning gene structures and dumpster");
                cleanGenes(genes, gp->nclass, dumpster);
            }

            /* showing global stats about last sequence predicted */
            OutputStats(contig_name);

            /* Reset evidence temporary counters */
            if (EVD && evidence != NULL) {
                resetEvidenceCounters(external);
            }

            cleanAcc(m);
            strcpy(contig_name, next_contig_name);
        } /* endwhile(reading): next sequence to be processed... */
    } /*endifgeneid*/
    else {
        /* B. Only assembling genes from input exons */

        /* B.0. Reading DNA sequence to make the translations */
        /* open the Sequence File */
        if ((fasta_fptr = fopen(fasta_fn, "rb")) == NULL) {
            printError("The Sequence file can not be open for read");
            exit(EXIT_FAILURE);
        }

        printMess("Reading DNA sequence");
        reading = IniReadSequence(fasta_fptr, contig_name);

        if (reading != EOF) {
            reading        = ReadSequence(fasta_fptr, Sequence, next_contig_name);
            LengthSequence = FetchSequence(Sequence, RSequence);
        }

        /* Header Output */
        OutputHeader(contig_name, LengthSequence);

        /* B.1. Reading exons in GFF format */
        printMess("Reading exonsGFF from file");
        external->nvExons
            = ReadExonsGFF(exons_gff_fn, external, isochores[0]->D);
        sprintf(mess, "%ld exons acquired from file\n",
                external->nvExons);
        printMess(mess);

        if (external->nSequences > 1) {
            sprintf(mess, "Exons in more than one different locus were detected (%ld sequences)\n", external->nSequences);
            printError(mess);
        }

        if (external->nvExons > 0) {
            /* B.2. Calling to genamic for assembling the best gene */
            genamic(external->evidence[0]->vExons,
                    external->evidence[0]->nvExons,
                    genes, isochores[0]);
        }

        /* B.3. Printing gene predictions */
        OutputGene(genes, external->evidence[0]->nvExons,
                   contig_name, Sequence, isochores[0], dAA, GenePrefix);
    } /* end only gene assembling from exons file */

    /* CHECK_LEAKS(); */

    /* 4. The End */
    OutputTime();
    exit(EXIT_SUCCESS);

}
