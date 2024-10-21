
tx2gene = "~/Documents/phd/research_lines/rbp_bed/gencode.v42.tx2gene.csv"
column_name = "condition"
baseline = 0
contrast = 1
#controls_name = "EV"
#contrast_name = "USP10_OE"
#column_name = opt$column_name

#output_path=paste0(outputdir,controls_name,"-",contrast_name,".")

# ============================ section 1: import data ===========================

tx2gene <- data.table::fread(tx2gene,header=FALSE)
tx2gene <- tx2gene[,c(1:2)]
colnames(tx2gene) = c("TXNAME", "GENEID")

#(1) First read in the metadata. if only a subset of the files are used, the opt$pattern option will be taken.

salmon_deseq2 <- function(sample_dir, metadata_dir, experiment_dir) {
  metadata <- read.csv(metadata_dir)
metadata2 = metadata %>% 
  dplyr::filter(experiment == experiment_dir) %>% 
  #dplyr::select(sample, !!as.symbol(column_name)) %>% 
  mutate(comparison_condition = case_when(!!as.symbol(column_name) == baseline ~ 'baseline',
                                          !!as.symbol(column_name) == contrast ~ 'contrast',
                                          TRUE ~ NA_character_)) %>% 
  dplyr::filter(!is.na(comparison_condition))

##
    
#(2) Generate a vector of the wanted file names.
files = unique(paste0(sample_dir,"/salmon_quant_",experiment_dir,"/",metadata2$sample,"/quant.sf")) 
names(files) = unique(metadata2$sample)


#(3) To check if all the files exist
if(all(file.exists(files)) == FALSE) {
  stop("It seems that I cannot find those files...Please check if your directory is correct.")
}

# ====================== section 3: import salmon files ==============================
# files is a vector of directory where quant.sf file locates.
# just ignore the version... to make it easier for following steps.
txi.tx <- tximport(files, 
                   type="salmon", 
                   tx2gene=tx2gene,
                   ignoreTxVersion = TRUE,
                   ignoreAfterBar = TRUE,
                   txOut = TRUE, dropInfReps = T)

txi.sum <- summarizeToGene(txi.tx, tx2gene)

keep <- rowSums(edgeR::cpm(txi.sum$counts) > 5) >= 2
print("Filtered Genes by CPM greater than 5 in a least 2 samples")
print(table(keep))
txi.sum$counts <- txi.sum$counts[keep, ]
txi.sum$abundance <- txi.sum$abundance[keep, ]
txi.sum$length <- txi.sum$length[keep, ]

# make it csv
TPM_transcripts = as.data.frame(txi.tx$abundance) %>% 
  tibble::rownames_to_column(.,var="transcript_id")
TPM_gene = as.data.frame(txi.sum$abundance) %>% 
  tibble::rownames_to_column(.,var="gene_id")

#write.csv(TPM_transcripts,"TPM_transcripts.csv")
#write.csv(TPM_gene,"TPM_gene.csv")


# ========================================== section 4: RUN A DEFAULT DESEQ 2 (optional) =============================================================


dds = DESeqDataSetFromTximport(txi.sum,
                               colData = metadata2,
                               design = ~ comparison_condition) 



# 'Note that the tximport-to-DESeq2 approach uses estimated gene counts from the transcript abundance quantifiers, but not normalized counts' -- <Deseq2 vignette> (just a note - Deseq() wraps the normalization step inside)
# perform the Deseq function
dds = DESeq(dds)

# Now, extract the result and named them by their contrast group
results_table <<- results(dds) %>%
  as.data.frame() %>%
  tibble::rownames_to_column('Geneid') %>%
  mutate(ensgene = gsub("\\..*", "", Geneid)) %>%
  left_join(annotables::grch38[c(1,3)])
#assign(paste0("results_table_", experiment_dir), results_table)
#results_table ->> paste0("results_table",experiment)  
# Now, extract the DESeq2 normed counts
normed_counts <<- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column('Geneid') %>%
  mutate(ensgene = gsub("\\..*", "", Geneid)) %>%
  left_join(annotables::grch38[c(1,3)])
#summary(normed_counts)
#assign(paste0("normed_counts_", experiment_dir), normed_counts)
#return(results_table)

}

#write.csv(results_table,"DESeq2_results.csv")
#write.csv(normed_counts,"DESeq2_normalized_counts.csv")

