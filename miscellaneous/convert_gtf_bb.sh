#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## convert gtf to bigbed file for upload to UCSC (assume hg38)
## convert bed to bigBed file for upload to UCSC (not assume hg38, require chrom input
## only necessary for big gtf files 
## --------------------------------

source activate Spatial

# convertGTF2BigBed <full_path_to_gtf>
convertGTF2BigBed(){

	
	gtf=$1
	outputdir=$(dirname $gtf)
	name=$(basename $gtf .gtf)
	
	echo "Converting $gtf"
	echo "Output $name"
	echo "Output dir $outputdir"
	
	cd $outputdir
	# Convert Gtf to genePred
	gtfToGenePred ${gtf} ${name}.genePred
	# Convert genPred to bed12
	genePredToBed ${name}.genePred ${name}.bed12
	# sort bed12
	sort -k1,1 -k2,2n ${name}.bed12 > ${name}.sorted.bed12
	# fetch chrom sizes
	fetchChromSizes hg38 > hg38.chrom.sizes
	# Convert sorted bed12 to bigBed (useful for trackhubs)
	bedToBigBed ${name}.sorted.bed12 hg38.chrom.sizes ${name}.sorted.bb
	
	rm hg38.chrom.sizes

}

# convertBed2BigBed <input_file> <chrom> <output_name> <output_dir>
# recreate bb file with different name
convertBed2BigBed(){

	echo "Converting $1"
	
	cd $4
  # sort bed file
	sort -k1,1 -k2,2n $1 > $3.sorted.bed
  # Convert sorted bed12 to bigBed (useful for trackhubs) 
	bedToBigBed $3.sorted.bed $2 $3.bb
}

