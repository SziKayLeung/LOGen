#### -------------------- post SQANTI -------------------

## ---------- run_cpat -----------------

# run_cpat <input_fasta> <output_name> <output_dir>
# Aim: 
  # call ORF from fasta file using CPAT (determine whether isoforms are protein-coding or non-protein-coding)
# Input:
  # input_fasta = input fasta for ORF to be called from
  # output_name = prefix output name
  # output_dir = path of output root directory to create CPAT folder directory
# Pre-requisite:
  # ${HEXAMER} = CPAT hexamer file (called from config file)
  # ${LOGITMODEL} = CPAT logit model (called from config file)
# Output
  # CPAT output files
  # CPAT log file

run_cpat(){
  
  mkdir -p $3/CPAT; cd $3/CPAT
  
  cpat.py --version
  cpat.py -x ${HEXAMER} -d ${LOGITMODEL} -g $1 --min-orf=50 --top-orf=50 -o $2 2> $2"_cpat.e"

  
}



## ---------- extract_best_orf -----------------

# extract_best_orf <output_name> <input/output_dir>
# Aim: 
  # extract the best ORF from CPAT for further analysis of ORF predictions for predicted NMD
# Input:
  # output_name = input and output name used from CPAT analysis
  # input/output_dir = directory of CPAT files
# Pre-requisite:
  # run_cpat to generate CPAT output files
# Output:
  # 

extract_best_orf(){
  
  cd $2/CPAT
  extract_fasta_bestorf.py --fa $1".ORF_seqs.fa" --orf $1".ORF_prob.best.tsv" --o_name $1"_bestORF" --o_dir $2 &> orfextract.log

}


## ---------- convert_gtf_bed12 -----------------

# convert_gtf_bed12 <input_gtf> 
# Aim:
  # convert gtf to bed12 file for downstream 
# Input: 
  # input_gtf = input gtf to be converted 
# Output:
  # path/to/original/directory/<sample>_sorted.bed12
convert_gtf_bed12(){
  
  # variables 
  output_dir="$(dirname $1)" 
  sample=${1%.gtf} # removes .gtf
  
  cd ${output_dir}
  
  gtfToGenePred $1 $sample.genePred
  genePredToBed $sample.genePred > $sample.bed12
  sort -k1,1 -k2,2n $sample.bed12 > $sample"_sorted.bed12"
  rm $sample.genePred $sample.bed12

}

## ---------- colour_by_abundance -----------------

# colour_by_abundance <cpat_name> <input_gtf> <input_abundance> <output_dir> <species=mouse/human>
# Aim:
  # generate multiple abundance file using input CPAT and expression using custom script 
  # custom script: ${LOGEN_ROOT}/merge_characterise_dataset/colour_transcripts_by_countandpotential.py
# Input:
  # cpat_name = CPAT input prefix names 
  # input_gtf = input gtf for conversion to bed12
  # input_abundance = path of input abundance file
  # output_dir = root of output directory for characterisation
  # species = mouse/human for determining CPAT threshold in script
# Output:
  # bed files
  
colour_by_abundance(){
  
  mkdir -p $4/bed12Files
  
  # convert gtf to bed12
  convert_gtf_bed12 $2
  
  # variables
  bed12=${2%.gtf}_sorted.bed12
  sample="$(basename $2)" 
  outputname=${sample%.gtf}
  echo $outputname
  
  colour_transcripts_by_countandpotential.py \
    --bed $bed12 \
    --cpat $4/CPAT/$1".ORF_prob.best.tsv" \
    --noORF $4/CPAT/$1".no_ORF.txt" \
    --a $3 \
    --o $outputname \
    --dir $4/bed12Files/ \
    --species $5
  
}


# subset_gene_reference <root_dir>
subset_gene_reference(){
  
  mkdir -p $1/TargetGenesRef
  
  source activate sqanti2_py3
  subset_reference_by_gene.py --r=${GENOME_GTF} --glist ${TGENES[@]} --o $1/TargetGenesRef
  
}


# run_transdecoder <name> <root_dir>
run_transdecoder(){
  source deactivate
  
  mkdir -p $2/Transdecoder; cd $2/Transdecoder
  
  TransDecoder.LongOrfs -t $2/CPAT/$1"_bestORF.fasta" &> transdecoder_longorf.log
  
  source activate sqanti2_py3
  hmmscan --cpu 8 --domtblout pfam.domtblout $PFAM_REF $1"_bestORF.fasta.transdecoder_dir"/longest_orfs.pep &> hmmscan.log
  
  source activate nanopore
  TransDecoder.Predict -t $2/CPAT/$1"_bestORF.fasta" --retain_pfam_hits pfam.domtblout --no_refine_starts &> transdecoder_predict.log
  sed '/^#/ d' < pfam.domtblout > pfam.domtblout.read
  
}
