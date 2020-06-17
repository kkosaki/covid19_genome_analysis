Variant calling pipeline for amplicon-based sequencing of the SARS-Cov-2 viral genome
===========

Following article has been has posted on medRxiv.

Clinical Utility of SARS-CoV-2 Whole Genome Sequencing in Deciphering Source of Infection

Toshiki Takenouchi, Yuka W. Iwasaki, Sei Harada, Hirotsugu Ishizu, Yoshifumi Uwamino, Shunsuke Uno, Asami Osada, Naoki Hasegawa, Mitsuru Murata, Toru Takebayashi, Koichi Fukunaga, Hideyuki Saya, Yuko Kitagawa, Masayuki Amagai, Haruhiko Siomi, Kenjiro Kosaki

doi: https://doi.org/10.1101/2020.05.21.20107599

(https://medrxiv.org/cgi/content/short/2020.05.21.20107599v1)
# System requirements

The system was tested on Cent OS 6.3 and Ubuntu 16.04LTS and 18.06LTS.

Conda may be helpful during installation of the required packages

https://github.com/conda/conda

## Reference SARS-Cov-2 sequence

MN908947.3.fasta was used as the reference

https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3

Wu,F., Zhao,S., Yu,B., Chen,Y.M., Wang,W., Song,Z.G., Hu,Y.,Tao,Z.W., Tian,J.H., Pei,Y.Y., Yuan,M.L., Zhang,Y.L., Dai,F.H.,Liu,Y., Wang,Q.M., Zheng,J.J., Xu,L., Holmes,E.C. and Zhang,Y.Z. A new coronavirus associated with human respiratory disease in China.
Nature 2020; 579:265-269

## Raw data generation

### PCR amplified using the ARITC primer set version 3

https://artic.networ/ncov-2019
https://artic.network/resources/ncov/ncov-amplicon-v3.pdf

Amplicon sequencing by Illumina MiSeq

## Software components used in the pipeline
### Resampling and quality control
#### seqtk

Resample fastq data when depth is excessively high

https://github.com/lh3/seqtk

Wei Shen, Shuai Le, Yan Li, and Fuquan Hu
SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q  File Manipulation
PLoS One. 2016; 11: e0163962.

#### fastp
Trim adapters

https://github.com/OpenGene/fastp

Chen S, Zhou Y, Chen Y, Gu J.
fastp: an ultra-fast all-in-one FASTQ preprocessor.
Bioinformatics. 2018;34:i884-i890.

#### ivar

Soft clip PCR primer sequences

https://github.com/andersen-lab/ivar

Grubaugh ND, Gangavarapu K, Quick J, Matteson NL, De Jesus JG, Main BJ, Tan AL, Paul LM, Brackney DE, Grewal S, Gurfield N, Van Rompay KKA, Isern S, Michael SF, Coffey LL, Loman NJ, Andersen KG.

An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar.
Genome Biol. 2019 Jan 8;20:8.

### Alignment

#### bwa version 0.7.17-r1188 or above

https://github.com/lh3/bwa

Li H, Durbin R.
Fast and accurate short read alignment with Burrows-Wheeler transform.
Bioinformatics. 2009;25:1754-60.

#### variant calling

samtools (version1.9 or above)

bcftools (version 1.9 or above)

Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup.
The Sequence Alignment/Map format and SAMtools.
Bioinformatics. 2009;25:2078-9.

#### ivar

https://github.com/andersen-lab/ivar

See above

### Annotation
#### SnpEff

Predicted effects on translated protein

http://snpeff.sourceforge.net/

Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM.
A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.
Fly (Austin). 2012 ;6:80-92.

Requires to build custom database for SARS-COV-2

Custom GFF3 file for building database

See atached GFF3 file.

Binary database

covid19

#### water

Pairwise sequence alignment of the sample sequence and the reference sequence base on the Smith-Waterman algorithm.
From the emboss package

https://github.com/kimrutherford/EMBOSS

Carver T, Bleasby A: The design of Jemboss: a graphical user interface to EMBOSS. Bioinformatics. 2003; 19: 1837-1843.

