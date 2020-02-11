
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/043eeacbb72e462d9604bfaed06d8ca8)](https://www.codacy.com/manual/darked89/geneidc?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=darked89/geneidc&amp;utm_campaign=Badge_Grade)


** Caveat: unofficial geneid repository**



# geneid (dev branch)
=======

* version:  1.4.5+

- [Synopsis](#synopsis)
- [Installation](#installation)
	- [Requirements](#requirements)
		- [Platforms](#platforms)
		- [Programs and libraries](#programs-and-libraries)
	- [Obtaining the sources](#obtaining-the-sources)
	- [Compilation](#compilation)
	- [Testing](#testing)
- [Usage](#usage)
	- [Inputs](#inputs)
	- [Outputs](#outputs)
	- [Basic](#basic)
	- [Typical](#typical)
	- [Advanced](#advanced)
	- [Best practices](#best-practices)
	- [Parameter files](#parameter-files)
	- [Command line options](#command-line-options)
- [Help](#help)
- [Legal](#legal)
	- [License](#license)
	- [Authors](#authors)
	- [Citation](#citation)




## Synopsis

Geneid is a gene prediction program in eukaryotic genomes using Hidden Markov Models to detect signals in the DNA sequence. Integrating predictions from multiple sources is also supported, but quite basic at this stage.


Installation, setup and basic usage of geneid is fairly easy. The command line options control the program output and the program behaviour.


## Installation
### Requirements
#### Platforms/compilers

The program is written in ANSI C, but started moving towards the C11 C standard. 
It compiles on multiple Linux/Unix system with default compilers. The lists below contain only a subset of platforms /compiler versions, just to give an idea about recent distributions on which the builds were tested.

- Linux: 
    - Manjaro 18.1 (M_18.1)
    - Ubuntu 18.04.1 (U_18.0)
#### Programs and libraries
Tested with following compilers:
- gcc: 
    - 9.1.0 (M_18.1)
    - 7.4.0 (U_18.0)
- clang:  8.0.1  (M_18.1)
- cc (from Oracle Solaris Studio): 12.6 (M_18.1)
- icc (Intel): 19.0.5.281 (U_18.0)

TODO
- Python for testing (future)

### Obtaining the sources

Options:
- git clone https://github.com/darked89/geneidc.git


TODO: 

- Release

```
tar -xvfz geneid.tar.gz
```

### Compilation
Go to geneid directory and type:

```
#gcc
make

#clang
make --file=Makefile.clang

#Oracle CC
make --file=Makefile.oracle_cc
```

Executable will be created as bin/geneid

### Testing


```
#help
geneid -h

#running with predictions

./bin/geneid  -GP param/dictyostelium.param ./samples/dict_1chr.fa

```


## Usage
### Memory requirements
With the default configuration, geneid requires:
* 31 Gb for human sized genome
* 9 Gb for D.melanogaster genome (140 Mbp)

For estimating the virtual memory needed for geneid:

`geneid -B genome.fa
`
If the required memory size exceeds the RAM + swap size, one can split the genome into individual chromosomes/contigs.

### Inputs
TODO
### Outputs
TODO
### Basic
To run geneid type:
`geneid -P parameter_filename genome.fasta`


### Typical
#### Setting environment variables

XXX 
Alternatively you can set the parameter file using the environment
variable GENEID. 
XXX not explained

#### [Using existing parameter file on a new genome]

TBD

### Advanced
TODO
- Training
- MEME
- splice branches
- other


### Best practices

TODO

### Parameter files
For gene prediction `geneid` needs various information about splice sites, nucleotide/various kmer frequencies, codon biases, etc. This is provided in species specific .param files. One (??) .param file is used for an individual geneid prediction run. Using the already created parameters from species close to the one with which are working often gives satisfactory results.

To download individual speciess parameter files (compressed with XXX):
TODO
XXX dir with .param files XXX (if possible rename all of these according to some rule. Plus get more comments inside how & when these were created)

To download all/reviewed? current parameter files:
XXX link to geneid_all_param.tar.gz

#### Creating param file for new species

TODO

### Command line options

#### Short help
```
	geneid	[-bdaefitnxszru]
		[-TDAZU]
		[-p gene_prefix]
		[-G] [-3] [-X] [-M] [-m]
		[-WCF] [-o]
		[-j lower_bound_coord]
		[-k upper_bound_coord]
		[-N numer_nt_mapped]
		[-O <gff_exons_file>]
		[-R <gff_annotation-file>]
		[-S <gff_homology_file>]
		[-P <parameter_file>]
		[-E exonweight]
		[-V evidence_exonweight]
		[-Bv] [-h]
		<locus_seq_in_fasta_format>
```


#### Long help
```
geneid [flags] <locus_seq_in_fasta_format>

-b: Output Start codons
-d: Output Donor splice sites
-a: Output Acceptor splice sites
-e: Output Stop codons
-f: Output Initial exons
-i: Output Internal exons
-t: Output Terminal exons
-n: Output introns
-s: Output Single genes
-x: Output all predicted exons
-z: Output Open Reading Frames

-T: Output genomic sequence of exons in predicted genes
-D: Output genomic sequence of CDS in predicted genes
-A: Output amino acid sequence derived from predicted CDS
-p: Prefix this value to the names of predicted genes, peptides and CDS

-G: Use GFF format to print predictions
-3: Use GFF3 format to print predictions
-X: Use extended-format to print gene predictions
-M: Use XML format to print gene predictions
-m: Show DTD for XML-format output

-j <coord>: Begin prediction at this coordinate
-k <coord>: End prediction at this coordinate
-N <num_reads>: Millions of reads mapped to genome
-W: Only Forward sense prediction (Watson)
-C: Only Reverse sense prediction (Crick)
-U: Allow U12 introns (Requires appropriate U12 parameters to be set in the parameter file)
-r: Use recursive splicing
-F: Force the prediction of one gene structure
-o: Only running exon prediction (disable gene prediction)
-O <exons_filename>: Only running gene prediction (not exon prediction)
-Z: Activate Open Reading Frames searching
-R <exons_filename>: Provide annotations to improve predictions
-S <HSP_filename>: Using information from protein sequence alignments to improve predictions

-u: Turn on UTR prediction. Only valid with -S option: HSP/EST/short read ends are used to determine UTR ends
-E: Add this value to the exon weight parameter (see parameter file)
-V: Add this value to the score of evidence exons
-P <parameter_file>: Use other than default parameter file (human)

-B: Display memory required to execute geneid given a sequence
-v: Verbose. Display info messages
-h: Show this help


```


## Help
- program homepage: FFF
- help desk (technical: installation, bugs, etc.): please use GitHub issue tracker
- help with running geneid: please use Biostars https://www.biostars.org/
- email: geneid@imim.es

## Legal

### License
GNU General Public License

### Authors
Code:
- Enrique Blanco
- Roderic Guigo
- Tyler Alito

Parameter files:
- Genis Parra
- Francisco Camara

Contributions:
- Josep F.Abril
- Moises Burset
- Xavier Messegue

Code cleanup:
- darked89

### Citation
If used for published results, please cite:
TODO

### Geneid publications
TODO
