#!/bin/bash

collate_longread_processed(){
  mkdir -p $WKD_ROOT/2_longread_processed; cd $WKD_ROOT/2_longread_processed
  cp $ISO_CLASSFILE $WKD_ROOT/2_longread_processed/
  cp $ISO_FASTA $WKD_ROOT/2_longread_processed/
  cp $ISO_GTF $WKD_ROOT/2_longread_processed/  
}

# prepare reference tables from reference files in config file
# generate files for downstream usage
prepare_reference_tables(){
  mkdir -p $WKD_ROOT/3_reference_tables; cd $WKD_ROOT/3_reference_tables

  echo "#*************************************** Preparing reference tables"
  python $LREAD/generate_reference_tables/src/prepare_reference_tables.py \
    --gtf $GENOME_GTF \
    --fa $GENOME_TRANSCRIPT_FASTA \
    --ensg_gene ensg_gene.tsv   \
    --enst_isoname enst_isoname.tsv \
    --gene_ensp gene_ensp.tsv \
    --gene_isoname gene_isoname.tsv \
    --isoname_lens isoname_lens.tsv \
    --gene_lens gene_lens.tsv \
    --protein_coding_genes protein_coding_genes.txt

  echo "Make gencode databases"
  python $LREAD/make_gencode_database/src/make_gencode_database.py \
    --gencode_fasta $GENOME_TRANSLATION_FASTA \
    --protein_coding_genes protein_coding_genes.txt \
    --output_fasta gencode_protein.fasta \
    --output_cluster gencode_isoname_clusters.tsv
}

summarise_longread_data(){
  mkdir -p $WKD_ROOT/4_longread; cd $WKD_ROOT/4_longread

  echo "#*************************************** Generate fasta file with six nucleotide frame translatation"
  python $LREAD/six_frame_translation/src/six_frame_translation.py \
    --iso_annot $ISO_CLASSFILE \
    --ensg_gene $WKD_ROOT/3_reference_tables/ensg_gene.tsv \
    --sample_fasta $ISO_FASTA \
    --output_fasta $NAME.6frame.fasta

  echo "#*************************************** Summarise long-read transcriptome"
  python $LREAD/transcriptome_summary/src/transcriptome_summary_ck.py \
    --sq_out $ISO_CLASSFILE \
    --ensg_to_gene $WKD_ROOT/3_reference_tables/ensg_gene.tsv \
    --enst_to_isoname $WKD_ROOT/3_reference_tables/enst_isoname.tsv \
    --len_stats $WKD_ROOT/3_reference_tables/gene_lens.tsv

}

call_orf(){
  echo "#*************************************** Calling open reading frames from long-read transcriptome data"
  mkdir -p $WKD_ROOT/5_calledOrfs; cd $WKD_ROOT/5_calledOrfs
  
  source activate sqanti2_py3 
  cpat.py \
    -x $HEXAMER \
    -d $LOGITMODEL \
    -g $ISO_FASTA \
    --min-orf=50 \
    --top-orf=50 \
    -o $NAME \
    1> $NAME"_cpat.output" \
    2> $NAME"_cpat.error"

  python $LREAD/orf_calling/src/orf_calling.py \
    --orf_coord $NAME".ORF_prob.tsv" \
    --orf_fasta $NAME".ORF_seqs.fa" \
    --gencode $GENOME_GTF \
    --sample_gtf $ISO_GTF \
    --pb_gene $WKD_ROOT/4_longread/pb_gene.tsv \
    --classification $ISO_CLASSFILE \
    --sample_fasta $ISO_FASTA \
    --num_cores 2 \
    --output $NAME"_best_orf.tsv"
  
  source deactivate 
  source activate lrp
}

refine_calledorf(){
  echo "#*************************************** Filter on called open reading frames"
  mkdir -p $WKD_ROOT/6_refined_database; cd $WKD_ROOT/6_refined_database

  python $LREAD/refine_orf_database/src/refine_orf_database.py \
    --name $NAME \
    --orfs $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv"  \
    --pb_fasta $ISO_FASTA \
    --coding_score_cutoff $coding_score_cutoff &> refine_org.log

  python $LREAD/visualization_track/src/make_pacbio_cds_gtf.py \
    --name $NAME \
    --sample_gtf $ISO_GTF \
    --refined_database $NAME"_orf_refined.tsv" \
    --called_orfs $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv" \
    --pb_gene $WKD_ROOT/4_longread/pb_gene.tsv \
    --include_transcript yes &> make_cds_gtf1.log

  # modified make_pacbio_cds_gtf.py to remove cpm
  python $LREAD/visualization_track/src/make_pacbio_cds_gtf.py \
    --name $NAME"_no_transcript" \
    --sample_gtf $ISO_GTF \
    --refined_database $NAME"_orf_refined.tsv" \
    --called_orfs $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv" \
    --pb_gene $WKD_ROOT/4_longread/pb_gene.tsv \
    --include_transcript no &> make_cds_gtf1.log   
}


classify_protein(){
  mkdir -p $WKD_ROOT/7_classified_protein; cd $WKD_ROOT/7_classified_protein

  python $LREAD/rename_cds_to_exon/src/rename_cds_to_exon.py \
    --sample_gtf $WKD_ROOT/6_refined_database/$NAME"_with_cds.gtf" \
    --sample_name $NAME \
    --reference_gtf $GENOME_GTF \
    --reference_name gencode \
    --num_cores 2 &> rename_cds_to_exon.log

  python $LREAD/sqanti_protein/src/sqanti3_protein.py $NAME.transcript_exons_only.gtf \
    $NAME.cds_renamed_exon.gtf \
    $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv" \
    gencode.transcript_exons_only.gtf \
    gencode.cds_renamed_exon.gtf \
    -d ./ \
    -p $NAME &> sqanti3_protein.log

  # modified script with capture_output=TRUE as only works for python 3.7 and above
  python $LREAD/5p_utr_status/src/1_get_gc_exon_and_5utr_info.py \
    --gencode_gtf $GENOME_GTF \
    --odir ./

  python $LREAD/5p_utr_status/src/2_classify_5utr_status.py \
    --gencode_exons_bed gencode_exons_for_cds_containing_ensts.bed \
    --gencode_exons_chain gc_exon_chain_strings_for_cds_containing_transcripts.tsv \
    --sample_cds_gtf $WKD_ROOT/6_refined_database/$NAME"_with_cds.gtf" \
    --odir ./

  python $LREAD/5p_utr_status/src/3_merge_5utr_info_to_pclass_table.py \
    --name $NAME \
    --utr_info pb_5utr_categories.tsv \
    --sqanti_protein_classification $NAME.sqanti_protein_classification.tsv \
    --odir ./

  python $LREAD/protein_classification/src/protein_classification_add_meta.py \
    --protein_classification $NAME.sqanti_protein_classification_w_5utr_info.tsv \
    --best_orf $WKD_ROOT/5_calledOrfs/$NAME"_best_orf.tsv" \
    --refined_meta $WKD_ROOT/6_refined_database/$NAME"_orf_refined.tsv" \
    --ensg_gene $WKD_ROOT/3_reference_tables/ensg_gene.tsv \
    --name $NAME \
    --dest_dir ./

  python $LREAD/protein_classification/src/protein_classification.py \
    --sqanti_protein $NAME.protein_classification_w_meta.tsv \
    --name $NAME"_unfiltered" \
    --dest_dir ./

  python $LREAD/protein_gene_rename/src/protein_gene_rename.py \
    --sample_gtf $WKD_ROOT/6_refined_database/$NAME"_with_cds.gtf" \
    --sample_protein_fasta $WKD_ROOT/6_refined_database/$NAME"_orf_refined.fasta" \
    --sample_refined_info $WKD_ROOT/6_refined_database/$NAME"_orf_refined.tsv" \
    --pb_protein_genes $NAME"_genes.tsv" \
    --name $NAME

  python $LREAD/protein_filter/src/protein_filter.py \
    --protein_classification $NAME"_unfiltered.protein_classification.tsv" \
    --gencode_gtf $GENOME_GTF \
    --protein_fasta $NAME.protein_refined.fasta \
    --sample_cds_gtf $NAME"_with_cds_refined.gtf" \
    --min_junctions_after_stop_codon $min_junctions_after_stop_codon \
    --name $NAME
}

run_hybrid_annotation(){
  mkdir -p $WKD_ROOT/8_hybrid_annotation; cd $WKD_ROOT/8_hybrid_annotation
  python $LREAD/make_hybrid_database/src/make_hybrid_database_ck.py \
    --protein_classification $WKD_ROOT/7_classified_protein/$NAME".classification_filtered.tsv" \
    --gene_lens $WKD_ROOT/3_reference_tables/gene_lens.tsv \
    --pb_fasta $WKD_ROOT/7_classified_protein/$NAME".filtered_protein.fasta" \
    --gc_fasta $WKD_ROOT/3_reference_tables/gencode_protein.fasta \
    --refined_info $WKD_ROOT/7_classified_protein/$NAME"_orf_refined_gene_update.tsv" \
    --pb_cds_gtf $WKD_ROOT/7_classified_protein/$NAME"_with_cds_filtered.gtf" \
    --name $NAME \
    --lower_kb $lower_kb \
    --upper_kb $upper_kb \
    --lower_cpm $lower_cpm
}

# run_metamorpheus <input_fasta>
run_metamorpheus(){
  mkdir -p $WKD_ROOT/9_metamorpheus; cd $WKD_ROOT/9_metamorpheus
  mkdir -p $WKD_ROOT/9_metamorpheus/$1; cd $WKD_ROOT/9_metamorpheus/$1

  if [ $1 == "refined" ]; then
    input_fasta=$WKD_ROOT/7_classified_protein/$NAME".protein_refined.fasta"
  elif [ $1 == "filtered" ]; then
    input_fasta=$WKD_ROOT/7_classified_protein/$NAME".filtered_protein.fasta"
  elif [ $1 == "hybrid" ]; then
    input_fasta=$WKD_ROOT/8_hybrid_annotation/$NAME"_hybrid.fasta"
  elif [ $1 == "gencode" ]; then
    input_fasta=$WKD_ROOT/3_reference_tables/gencode_protein.fasta
  elif [ $1 == "uniprot" ]; then
    input_fasta=$UNIPROT_FASTA
  else
    echo "Input paramater: <refined> <filtered> <hybrid> <gencode> <uniprot"
  fi

  echo "Processing $input_fasta"
  protein_data=$(for i in $PROTEIN_RAW/*raw*; do echo $i; done)
  echo $protein_data

  dotnet $METAMORPHEUS/CMD.dll -g -o ./toml --mmsettings ./settings
  yes | dotnet $METAMORPHEUS/CMD.dll -d $input_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $protein_data -t $METAMORPHEUS_TOML -v normal --mmsettings settings -o ./$NAME"_"$1"_search_results" &> $NAME"_"$1"_metamorpheus.log"
  mv $NAME"_"$1"_search_results"/Task1SearchTask/AllPeptides.psmtsv $NAME"_"$1"_search_results"/Task1SearchTask/AllPeptides.$1".psmtsv"
  mv $NAME"_"$1"_search_results"/Task1SearchTask/AllQuantifiedProteinGroups.tsv $NAME"_"$1"_search_results"/Task1SearchTask/AllQuantifiedProteinGroups.$1".tsv"
}

run_metamorpheus_rescue(){
  cd $WKD_ROOT/9_metamorpheus
  echo "Processing MetaMorpheus Rescue on hybrid fasta"
  dotnet $METAMORPHEUS/CMD.dll -d $WKD_ROOT/8_hybrid_annotation/$NAME"_hybrid.fasta" settings/Contaminants/MetaMorpheusContaminants.xml -s $protein_data -t $METAMORPHEUS_RESCUE_TOML -v normal --mmsettings settings -o ./$NAME"_rescue_search_results" --orf $WKD_ROOT/8_hybrid_annotation/$NAME"_refined_high_confidence.tsv" --cpm 25
  mv $NAME"_rescue_search_results"/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.$NAME.rescue_resolve.psmtsv
  mv $NAME"_rescue_search_results"/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.$NAME.rescue_resolve.tsv
}

run_peptide_analysis(){
  cd $WKD_ROOT/9_metamorpheus
  python $LREAD/peptide_analysis/src/peptide_analysis.py \
    -gmap $WKD_ROOT/3_reference_tables/gene_isoname.tsv \
    --gencode_peptides gencode/$NAME"_gencode_search_results"/Task1SearchTask/AllPeptides.gencode.psmtsv \
    --pb_refined_fasta $WKD_ROOT/7_classified_protein/$NAME".protein_refined.fasta" \
    --pb_filtered_fasta $WKD_ROOT/7_classified_protein/$NAME".filtered_protein.fasta" \
    --pb_hybrid_fasta $WKD_ROOT/8_hybrid_annotation/$NAME"_hybrid.fasta" \
    --pb_gene $WKD_ROOT/4_longread/pb_gene.tsv \
    -odir ./
}

# convert_gtf_bed12 <sample>.gtf
convert_gtf_bed12(){
  echo "Processing $1.gtf for conversion to bed12"
  gtfToGenePred $1.gtf $1.genePred
  genePredToBed $1.genePred $1.bed12
  if [ $3 == "make_region" ]; then
      # Squish intronic regions
      python $LREAD/visualization_track/src/make_region_bed_for_ucsc.py --name $2 --sample_gtf $1.gtf --reference_gtf $GENOME_GTF
  fi
}

generate_cds_tracks(){
  mkdir -p $WKD_ROOT/10_ucsc_tracks
  mkdir -p $WKD_ROOT/10_ucsc_tracks/cds; cd $WKD_ROOT/10_ucsc_tracks/cds
  
  echo "#*************************************** Generate CDS tracks"
  # Gencode
  python $LREAD/visualization_track/src/gencode_filter_protein_coding.py --reference_gtf $GENOME_GTF
  convert_gtf_bed12 gencode.filtered NA NA
  
  python $LREAD/visualization_track/src/gencode_add_rgb_to_bed.py \
    --gencode_bed gencode.filtered.bed12 \
    --rgb 0,0,140 \
    --version V35

  # Long Read
  convert_gtf_bed12 $WKD_ROOT/7_classified_protein/$NAME"_with_cds_refined" $NAME"_refined" make_region
  convert_gtf_bed12 $WKD_ROOT/7_classified_protein/$NAME"_with_cds_filtered" $NAME"_filtered" make_region
  convert_gtf_bed12 $WKD_ROOT/8_hybrid_annotation/$NAME"_cds_high_confidence" $NAME"_high_confidence" make_region
}

generate_peptide_tracks(){
  mkdir -p $WKD_ROOT/10_ucsc_tracks/peptides; cd $WKD_ROOT/10_ucsc_tracks/peptides

  # make_peptide <refined/filtered>
  make_peptide(){
    if [ $1 == "hybrid" ]; then
      refined_fasta_input=$WKD_ROOT/8_hybrid_annotation/$NAME"_hybrid.fasta"
      sample_gtf_input=$WKD_ROOT/7_classified_protein/$NAME"_with_cds_filtered.gtf"
    else
      refined_fasta_input=$WKD_ROOT/7_classified_protein/$NAME".protein_refined.fasta"
      sample_gtf_input=$WKD_ROOT/7_classified_protein/$NAME"_with_cds_"$1".gtf"
    fi

    python $LREAD/visualization_track/src/make_peptide_gtf_file.py \
      --name $NAME"_"$1 \
      --sample_gtf $sample_gtf_input \
      --reference_gtf $GENOME_GTF \
      --peptides $WKD_ROOT/9_metamorpheus/$1/$NAME"_"$1"_search_results"/Task1SearchTask/AllPeptides.$1".psmtsv" \
      --pb_gene $WKD_ROOT/7_classified_protein/$NAME"_genes.tsv" \
      --gene_isoname $WKD_ROOT/3_reference_tables/gene_isoname.tsv \
      --refined_fasta $refined_fasta_input

    convert_gtf_bed12 $WKD_ROOT/10_ucsc_tracks/peptides/$NAME"_"$1"_peptides" NA NA
  }

  make_peptide refined
  make_peptide filtered
  make_peptide hybrid
}

compare_protein_groups(){
  mkdir -p $WKD_ROOT/11_compare_proteins; cd $WKD_ROOT/11_compare_proteins
  python $LREAD/accession_mapping/src/accession_mapping.py \
    --gencode_fasta $WKD_ROOT/3_reference_tables/gencode_protein.fasta \
    --pacbio_fasta $WKD_ROOT/7_classified_protein/$NAME".protein_refined.fasta" \
    --uniprot_fasta $UNIPROT_FASTA

  python $LREAD/protein_groups_compare/src/protein_groups_compare.py \
    --pg_fileOne $WKD_ROOT/9_metamorpheus/gencode/$NAME"_gencode_search_results"/Task1SearchTask/AllQuantifiedProteinGroups.gencode.tsv \
    --pg_fileTwo $WKD_ROOT/9_metamorpheus/hybrid/$NAME"_hybrid_search_results"/Task1SearchTask/AllQuantifiedProteinGroups.hybrid.tsv \
    --mapping accession_map_gencode_uniprot_pacbio.tsv \
    --output ./

  python $LREAD/protein_groups_compare/src/protein_groups_compare.py \
    --pg_fileOne $WKD_ROOT/9_metamorpheus/gencode/$NAME"_gencode_search_results"/Task1SearchTask/AllQuantifiedProteinGroups.gencode.tsv \
    --pg_fileTwo $WKD_ROOT/9_metamorpheus/uniprot/$NAME"_uniprot_search_results"/Task1SearchTask/AllQuantifiedProteinGroups.uniprot.tsv \
    --mapping accession_map_gencode_uniprot_pacbio.tsv \
    --output ./
}

identify_novel_peptides(){
  mkdir -p $WKD_ROOT/12_novel_peptides; cd $WKD_ROOT/12_novel_peptides

  run_peptide_novelty_analysis(){
    python $LREAD/peptide_novelty_analysis/src/peptide_novelty_analysis.py \
      --pacbio_peptides $WKD_ROOT/9_metamorpheus/$1/$NAME"_"$1"_search_results"/Task1SearchTask/AllPeptides.$1".psmtsv" \
      --gencode_fasta $WKD_ROOT/3_reference_tables/gencode_protein.fasta \
      --uniprot_fasta $UNIPROT_FASTA \
      --name $NAME"_"$1
  }

  run_peptide_novelty_analysis refined
  run_peptide_novelty_analysis filtered
  run_peptide_novelty_analysis hybrid

}

