cat gencode.v40.annotation.gtf | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' |sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_type "//'| sed 's/gene_name "//' | sed 's/"//g' | awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' | sed "1i\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength"  >  gencode.v40.annotation.geneannotation.txt

grep "^>" homo_sapiens_uniprot.fasta | awk '{for(i=1;i<=NF;i++) if($i ~ /^GN=/) print $1, $i}' | tr "|" "\t" | tr " " "\t" > homo_sapiens_uniprot.txt

awk 'BEGIN{FS="\t"} 
{
    if($3 == "transcript") {
        split($9, a, "; ");
        gene_id = "";
        transcript_id = "";
        for(i in a) {
            if(a[i] ~ /gene_id/) {
                split(a[i], b, " ");
                gene_id = b[2];
                gsub(/"/, "", gene_id);
            }
            if(a[i] ~ /transcript_id/) {
                split(a[i], c, " ");
                transcript_id = c[2];
                gsub(/"/, "", transcript_id);
            }
        }
        if(gene_id != "" && transcript_id != "") {
            print transcript_id "\t" gene_id;
        }
    }
}' gencode.v40.annotation.gtf > gencode.v40.annotation.gene_transcript_ids.txt