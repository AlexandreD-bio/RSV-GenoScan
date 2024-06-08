# RSV-GenoScan

RSV-GenoScan is a tool for analysing RSV genomes from clinical samples, sequenced by Illumina (single-read or Paired-end) or Nanopore platforms.
The algorithm is able to generate the whole genome sequence, to determine the viral group (and subgroup), and to identify mutations in F and G glycoproteins (including drug resistance mutations).

# Prerequisites
- Linux Ubuntu 20.04 LTS (or higher)
- Python v3.9 (or higher)

# Installation
To make the program easier to use and install, once you have downloaded all the files from the repository, you will need to run the ``packets_installation.sh`` script. It will check and install, if necessary, the various packages/software the program needs to function correctly.
You can easily download RSV-GenoScan and install all dependensies with the following commands:

Via HTTPS link :
```
git clone https://github.com/AlexandreD-bio/RSV-GenoScan.git
```

```
cd RSV-GenoScan/Script/
bash packets_installation.sh
```

# Usage
Once installed, you can easily run the program with the following command:

```
bash Bash_struct.sh
```

Then, just follow the instructions until the program ends.

Note: For the phylogenetic analysis, additional reference genomes can be added to the existing panel by pasting their fasta files to the ``references_phylogeny`` repository (files ref_A.fasta and ref_B.fasta). 

# Outputs
The pipeline:
- extracts the percentage of genome coverage (at 10x)
- generates graphs showing the distribution of reads along the genome
- determines the viral group (A or B)
- determines the viral subgroup (based on the classification made by Goya et al.) using a phylogenetic approach
- generates the genome consensus sequence in .fasta format
- extracts coverage data for F and G glycoproteins
- generates the F and G glycoprotein sequences in .fasta format
- detects the presence or absence of G protein duplication
- lists the mutations in both F and G glycoproteins
- detects mutations in F and G known to confer drug resistance

# Tutorial
A video tutorial has been created for easy installation and execution of RSV-GenoScan:  https://youtu.be/8LQlHOGjkfI 

# Citation
RSV-GenoScan: An automated pipeline for whole-genome human respiratory syncytial virus (RSV) sequence analysis
Alexandre Dosbaa, Romane Guilbaud, Anna-Maria Franco Yusti, Valentine Marie Ferré, Charlotte Charpentier, Diane Descamps, Quentin Le Hingrat, Romain Coppée 
https://doi.org/10.1016/j.jviromet.2024.114938
