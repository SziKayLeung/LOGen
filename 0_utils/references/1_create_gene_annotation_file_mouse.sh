cat gencode.vM22.annotation.gtf | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' |sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_type "//'| sed 's/gene_name "//' | sed 's/"//g' | awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' | sed "1i\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength"  >  gencode.vM22.annotation.geneannotation.txt

# command to extract ENSEBML gene id and gene names from annotation.gtf
paste <(grep -o -P 'gene_id "([^"]+)"' gencode.v40.annotation.gtf | awk '{gsub(/"/, ""); print $NF}') \ <(grep -o -P 'gene_name "([^"]+)"' gencode.v40.annotation.gtf | awk '{gsub(/"/, "");
print $NF}') > gencode.v40.annotation_ENSEMBL_geneNames.txt
sort -u gencode.v40.annotation_ENSEMBL_geneNames.txt -o gencode.v40.annotation_ENSEMBL_geneNames.txt
