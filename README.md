Welcome to the **LO**ng-read **Gen**eral github repo!

This is a GitHub repo of (mostly independent) Python/R scripts that I developed to analyse data from long-read sequencing experiments. Purpose of scripts vary from generating txt files to run community tools (example pipelines), generating plots post-SQANTI, running differential expression analyses to more custom applications.  

#
###  [Processing ONT raw data](https://github.com/SziKayLeung/LOGen/wiki/Processing-ONT-raw-data) 
A pipeline for processing raw ONT reads from transcriptome cDNA processing, using research community tools (i.e. _Porechop_,_Minimap2_,_SQANTI3_) and own custom scripts. 

#
###  Data exploration post-SQANTI
Below listed are features that can be explored on `<sample>_classification.txt` generated from [SQANTI](https://github.com/ConesaLab/SQANTI3). 

*  [number of isoforms by structural category](https://github.com/SziKayLeung/LOGen/blob/master/transcriptome_stats/plot_basic_stats.R) 
*  [correlate exon number, gene length with isoform number](https://github.com/SziKayLeung/LOGen/blob/master/transcriptome_stats/corr_exon_length_num.R)
*  [identify long-non-coding RNA isoforms](https://github.com/SziKayLeung/LOGen/blob/master/transcriptome_stats/identify_lncRNA.R)
*  [plot and test the number of isoforms with/without certain features](https://github.com/SziKayLeung/LOGen/blob/master/transcriptome_stats/plot_hist_cage_SS_peaks.R) (i.e. within/without 50bp of CAGE peak/TSS/TTS)

To run functions, read in `<sample>_classification.txt` file using:
* `SQANTI_class_preparation(<sample>_classification.txt, standard)` if expression columns are included in the file (after running --FL_count in SQANTI)
* `SQANTI_class_preparation(<sample>_classification.txt, nstandard)` if expression is not included
#

###  [Characterize merged datasets](https://github.com/SziKayLeung/LOGen/wiki/Characterize-merged-datasets)
* `subset_targetgenes_classfiles.py`: Subset SQANTI classification file based on genes and reads
* `colour_transcripts_by_countandpotential.py`: Colour bed file by abundance and coding potential
* `extract_fasta_bestorf.py`: Create a fasta file based on best ORF defined from CPAT
#

###  Differential expression analysis
Current script dump to maintain. Scripts to input results after running tappAS, running linear regression etc...
#
 
###  [Miscellaneous](https://github.com/SziKayLeung/LOGen/wiki/Miscellaneous-Scripts)
* `replace_filenames_with_csv.py`: Replace multiple file names in a directory using reference csv file  
* `search_fasta_by_sequence.py`: Subset fasta based on sequence 
* `subset_fasta_gtf.py`: Subset gtf, fasta and bed files based on list of transcript IDs 
