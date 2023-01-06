#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: characterise/compare Iso-Seq vs ONT targeted datasets 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## split_iso_ont_abundance
## IsoSeq_vs_ONT_descriptive
## characterise_iso_ont_merge


## ------------------- split_iso_ont_characteristics

# Aim: summarise the characteristics of isoforms that are unique/common to merged Iso-Seq&ONT dataset 
# extract the original characteristics of isoforms before they were merged into one dataset
# Note: 
  # Requires input to have been processed through dataset_identifer.R to identify the isoform Id of merged datasets
# Input:
  # merged.class.files = classification file of merged datsaet with "ONT_isoform" and "IsoSeq_isoform" columns
  # iso.class.files.1 = original Iso-Seq classification file with "FL" column populated
  # ont.class.files.2 = original ONT classification file with "FL" column populated
# Output:
  # merged.class.files.exp = df with columns: Dataset, isoform, Length, Exons, FL, and LogFL for each dataset

split_iso_ont_characteristics <- function(merged.class.files, iso.class.files.1, ont.class.files.2){
  
  merged.class.files.exp <- merged.class.files[,c("Dataset","ONT_isoform","IsoSeq_isoform")] %>%  
    full_join(iso.class.files.1[,c("isoform","length","exons","FL")], by = c("IsoSeq_isoform" = "isoform")) %>% 
    full_join(ont.class.files.2[,c("isoform","length","exons","FL")], by = c("ONT_isoform" = "isoform")) %>%
    `colnames<-`(c("Dataset","ONT_isoform","IsoSeq_isoform", "IsoSeq_Length","IsoSeq_Exons",
                   "IsoSeq_FL","ONT_Length","ONT_Exons","ONT_FL")) %>% 
    mutate(IsoSeq_LogFl = log10(IsoSeq_FL), ONT_LogFL = log10(ONT_FL))
  
  return(merged.class.files.exp)
}


## ------------------- IsoSeq_vs_ONT_descriptive

# Aim: wrapper to plot Iso-Seq vs ONT descriptive characteristics
# Note:   
  # Function is to compare to different datasets that have been merged using rbind 
  # the two datasets are not merged therefore there may be isoforms that are commonly detected in both
  # this is just to test what are the characteristics of the two different datasets
# Input: 
  # bound.class.files = df: SQANTI classification file of rbind datasets with "Dataset" column as identifier
# Functions:
  # geom_density_plot() from draw_density.R
  # plot_stats_feature_2datasets() from identify_within_cage_SS_peaks.R
  # plot_polyA_motifs() from plot_hist_cage_SS_peaks.R
  # corr_gene_exp_toknown() from base_comparison.R
  # length_exon_description() from plot_basic_stats.R
# Output: plots between Iso-Seq vs ONT
  # p1: histogram of transcript length
  # p2: histogram of number of exons
  # p3: bar-plot of the percentage of isoforms within 50bp of CAGE peak
  # p4: bar-plot of the percentage of isoforms with 50bp of TTS (transcription termination site)
  # p5: bar-plot of the percentage of isoforms with 50bp of TSS (transcription start site)
  # p6: bar-plot of the percentage of isoforms with 3'ISM (degraded products)
  # p7: bar-plot of the percentage of isoforms within polyA site (split by structural category)
  # p8: bar-plot of the percentage of isoforms within 5-bp of CAGE peak (split by structural category)
  # p9: bar-plot of the number of isoforms with polyA motif 
  # p10: density plot of the gene expression of Iso-Seq vs RNA-Seq expression
  # p11: density plot of the gene expression of ONT vs RNA-Seq expression
  # p12: density plot of the number of isoforms between Iso-Seq vs ONT

IsoSeq_vs_ONT_descriptive <- function(bound.class.files){
  
  p1 <- geom_density_plot(bound.class.files, "length","Transcript Length (bp)") + 
    scale_fill_manual(values = c(label_colour("IsoSeq"),label_colour("ONT")), labels = c("Iso-Seq","ONT"))
  
  p2 <- geom_density_plot(bound.class.files, "exons","Number of Exons")  +
    scale_fill_manual(values = c(label_colour("IsoSeq"),label_colour("ONT")), labels = c("Iso-Seq","ONT"))
  
  p3 <- plot_stats_feature_2datasets(bound.class.files, "within_cage_peak", "Cage Peak")   
  
  p4 <- plot_stats_feature_2datasets(bound.class.files, "within_50_TTS", "TTS")
  
  p5 <- plot_stats_feature_2datasets(bound.class.files, "within_50_TSS", "TSS")
  
  p6 <- plot_stats_feature_2datasets(bound.class.files, "ISM_3prime", "3'ISM")
  
  p7 <- plot_features_byisocate(bound.class.files, "within_polya_site", "PolyA Site")
  
  p8 <- plot_features_byisocate(bound.class.files, "within_cage_peak", "Cage Peak") 
  
  p9 <- plot_polyA_motifs(bound.class.files)
  
  p10 <- corr_gene_exp_toknown(IsoSeq_filtered_class, "Iso-Seq")
  
  p11 <- corr_gene_exp_toknown(ONT_retained_class, "ONT")
  
  p12 <- bound.class.files %>% group_by(Dataset, associated_gene) %>% tally() %>% 
    spread(., Dataset, n) %>% 
    mutate(IsoSeq = log10(`Iso-Seq`), ONT = log10(ONT)) %>%
    density_plot(., "IsoSeq","ONT", "Number of Iso-Seq Isoforms detected (Log10)", "Number of ONT Isoforms detected (Log10)", "")
  
  # stats paper reporting
  length_exon_description(bound.class.files[bound.class.files$Dataset == "Iso-Seq",])
  
  length_exon_description(bound.class.files[bound.class.files$Dataset == "ONT",])
  
  return(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12))
  
}


## ------------------- characterise_iso_ont_merge

# Aim: characterise the merged transcriptome 
# compare characteristics of isoforms that are common to both Iso-Seq and ONT datasets, or unique to each dataset
# Note:
  # 3 output plots generated rely on RNA-Seq input to SQANTI classification files (i.e input to min_cov column)
# Functions:
  # plot_iso_length_mdatasets() from base_comparison.R
  # density_plot() from draw_density.R
  # rnaseq_support_mdatasets() from base_comparison.R
  # common_vs_unique_iso() from base_comparison.R
# Input:
  # merged.class.files = df: classification file of merged dataset using gffcompare for comparison, hence no redundant isoforms
  # merged.class.files.exp = df: df with expression of the isoforms, "FL", and "LogFL" 
# Output:
  # p1: box-plot of the Iso-Seq/ONT isoform expression of isoforms common or unique to datasets
  # p2: box-plt of the isoform length of isoforms common or unique to datasets
  # p3: box-plt of the number of exons in isoforms common or unique to datasets
  # p4: density plot (correlation) of the transcript length of common isoforms in both datasets
  # p5: density plot (correlation) of the exon number of common isoforms in both datasets
  # p6: density plot (correlation) of the abundance (FL) of common isoforms in both datasets
  # p7: bar-plot of the number of isoforms with RNA-Seq support in both or unique to dataset
  # p8: bar-plot of the percentage of isoforms with RNA-Seq support in both or unique to dataset
  # p9: bar-plot/stats of number of isoforms unique to dataset
  # p10: box-plot of isoform expression (Log10FL) of common and unique isoforms supported by RNA-Seq
  # p11: box-plot of isoform expression (Log10FL) of common and unique isoforms supported by CAGE peak

characterise_iso_ont_merge <- function(merged.class.files, merged.class.files.exp){
  
  # stats tests 
  print(wilcox.test(IsoSeq_FL ~ Dataset, data = merged.class.files.exp, exact = FALSE))
  print(wilcox.test(ONT_FL ~ Dataset, data = merged.class.files.exp, exact = FALSE))
  wilcox.test(merged.class.files.exp[merged.class.files.exp$Dataset == "Both",c("IsoSeq_FL")],
              merged.class.files.exp[merged.class.files.exp$Dataset == "Iso-Seq",c("IsoSeq_FL")])
  
  for(i in c("ONT","Iso-Seq")){
    cat("#### Stats test for", i, "\n")
    df = merged.class.files.exp %>% filter(Dataset %in% c("Both",i))
    cat("#### Length ##### \n")
    res = wilcox.test(length ~ Dataset, data = df, exact = FALSE)
    print(res)
    print(res$p.value)
    cat("#### Exons ##### \n")
    res = wilcox.test(exons ~ Dataset, data = df, exact = FALSE)
    print(res)
    print(res$p.value)
  }
  
  # isoforms detected in both datasets
  common_iso = merged.class.files.exp %>% filter(Dataset == "Both")
  cat("Number of isoforms commonly detected:",nrow(common_iso),"\n")
  
  p1 <- merged.class.files.exp %>% dplyr::select(Dataset,IsoSeq_FL,ONT_FL) %>% reshape2::melt(id = "Dataset") %>% 
    filter(!is.na(Dataset)) %>% filter(!is.na(value)) %>%
    ggplot(., aes(x = Dataset, y = log10(value), fill = Dataset)) + geom_boxplot() + ylim(0,6) + 
    facet_grid(~variable, scales = "free_x", labeller = as_labeller(c(`IsoSeq_FL` = "Iso-Seq Expression",`ONT_FL` = "ONT Expression"))) +
    mytheme + theme(legend.position = "none") + labs(x = "", y = "Transcript Expression") +
    scale_fill_manual(values = c(label_colour ("BothTech"),label_colour("IsoSeq"),label_colour("ONT")))
  
  p2 <- plot_iso_length_mdatasets(merged.class.files, "length", "Length (bp)") +
    scale_fill_manual(values = c(label_colour ("BothTech"),label_colour("IsoSeq"),label_colour("ONT"))) 
  
  p3 <- plot_iso_length_mdatasets(merged.class.files, "exons", "Exon Number") +
    scale_fill_manual(values = c(label_colour ("BothTech"),label_colour("IsoSeq"),label_colour("ONT"))) 
  
  p4 <- density_plot(common_iso, "IsoSeq_Length","ONT_Length","Iso-Seq Transcript Length","ONT Transcript Length (bp)", "")
  
  p5 <- density_plot(common_iso, "IsoSeq_Exons","ONT_Exons","Iso-Seq Transcript Exons","ONT Transcript Exon Number", "")
  
  p6 <- density_plot(common_iso, "IsoSeq_LogFl","ONT_LogFL","Iso-Seq Transcript Expression","ONT Transcript Expression", "")
  
  p7 <- rnaseq_support_mdatasets(merged.class.files)[[1]]
  
  p8 <- rnaseq_support_mdatasets(merged.class.files)[[2]]
  
  p9 <- rnaseq_support_mdatasets(merged.class.files)[[3]]
  
  p10 <- common_vs_unique_iso(merged.class.files,merged.class.files.exp, "RNASeq_supported", "RNA-Seq","ONT") 
  
  p11 <- common_vs_unique_iso(merged.class.files,merged.class.files.exp, "within_cage_peak", "CAGE Peak","ONT") 
  
  return(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11))
}