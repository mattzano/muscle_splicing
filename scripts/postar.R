#Biostrings::getSeq(BSgenome, transcripts)...
#Biostrings::translate...
library(data.table)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
#devtools::install_github("dzhang32/ggtranscript")
library(ggtranscript)
library(patchwork)

postar_bed <- fread("~/Documents/phd/research_lines/rbp_bed/human_RBP_binding_sites_sorted.bed", 
                    col.names = c("seqnames", "start", "end", "dataset", "score", "strand", "QC"))

ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
humandb <- biomaRt::getBM(attributes = c("external_gene_name", 
                                         "ensembl_gene_id", "ensembl_transcript_id", "transcript_appris", 
                                         "uniprotswissprot", "uniprotsptrembl", "uniprot_isoform"), 
                          mart = ensembl)
princ <- humandb[which(humandb$transcript_appris == "principal1" |
                       humandb$transcript_appris == "principal2" | 
                       humandb$transcript_appris == "principal3" | 
                       humandb$transcript_appris == "principal4" |
                       humandb$transcript_appris == "principal5"), ]

cds_regions = cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "tx",use.names = TRUE)
cds_regions = unlist(cds_regions)
cds_regions$transcript_id = gsub("\\..*", "", names(cds_regions))

exons_regions = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "tx",use.names = TRUE)
exons_regions = unlist(exons_regions)
exons_regions$transcript_id = gsub("\\..*", "", names(exons_regions))

#my_clean_reader = function(file){
#  as.data.table(janitor::clean_names(fread(file)))
#}
#my_name_fixer = function(tbl){
#  colnames(tbl) = gsub(colnames(tbl)[[9]],"psi",colnames(tbl))
#  return(tbl)
#}

#gene_target = "STMN2"

other_loci <- read.table(here::here("data/other_loci.txt"), sep = "\t", header = T)

splicing_bedss <- splicing_bedparse_invivo[,c(1:7)] %>% 
  distinct(V2, .keep_all = T) %>% 
  filter(V2 != "158062221") %>% 
  rbind(other_loci) %>% 
  filter(V7 != "SLC24A3")

#violin_plot<- list()
transcript_bind_plot <- function(gene_target) {
#for (i in of_interest){  
  cds_parent = as_tibble(cds_regions[cds_regions$transcript_id %in% princ$ensembl_transcript_id[princ$external_gene_name == gene_target]])
  exons_parent = as_tibble(exons_regions[exons_regions$transcript_id %in% princ$ensembl_transcript_id[princ$external_gene_name == gene_target]])
  
  parent_cryptic_delta <- splicing_bedss %>%
    mutate(seqnames = V1,
           start = V2,
           end = V3) %>% 
    dplyr::filter(seqnames %in% exons_parent$seqnames & start > min(exons_parent$start) & end < max(exons_parent$end)) #%>%

parent_postar_bed <- postar_bed %>%
  dplyr::filter(seqnames %in% exons_parent$seqnames & start > min(exons_parent$start) & end < max(exons_parent$end)) %>%
  dplyr::filter(seqnames %in% parent_cryptic_delta$seqnames & start > min(parent_cryptic_delta$start)-100 & end < max(parent_cryptic_delta$end)+100) %>%
  mutate(RBP = paste0(".",word(dataset, 1, sep = "_"))) %>%
  #group_by(RBP) %>% dplyr::filter(n() > 3) %>%
  mutate(tdp = ifelse(RBP == ".TARDBP", "yes", "zno"))
  #dplyr::filter(RBP == ".TARDBP")
intronic <- to_intron(exons_parent, group_var = "transcript_id")

plotz <- ggplot(aes(xstart = start, xend = end, y = gene_target), data = exons_parent) +
  geom_range(data = exons_parent, fill = "white", height = 0.2) +
  geom_range(data = cds_parent, fill = "black", height = 0.4) +
  #geom_intron(data = to_intron(exons_parent), aes(strand = strand)) +
  geom_intron(data = intronic, aes(strand = strand)) +
  geom_junction(data = parent_cryptic_delta, junction.orientation = "top", junction.y.max = 1) +
  geom_range(aes(y=RBP, color = tdp, fill = tdp), data = parent_postar_bed, 
             #color = "#377EB8", fill = "#377EB8", 
             height = 0.3, show.legend = F) +
  ggforce::facet_zoom(xlim = c(parent_cryptic_delta$start-100, parent_cryptic_delta$end+100)) + #xlim = c(17627940, 17642845)) +
  ylab("") +
  scale_y_discrete(expand = c(0,2)) +
  scale_color_manual(values = c("darkred", "navy")) +
  scale_fill_manual(values = c("darkred", "navy")) +
  theme_classic() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
#plot(plotz)

if (nrow(parent_cryptic_delta) > 0 & nrow(parent_postar_bed) > 0) {
  ploty <- plotz + geom_junction(data = parent_cryptic_delta, junction.orientation = "bottom", junction.y.max = 0.6) +
    annotate(geom = 'text', label = gene_target, x = -Inf, y = Inf, hjust = 0, vjust = 1)
  plot(ploty)
} else {
  print(paste0("MAJIQ found no cryptics in ", gene_target))
  #plot(plotz) + annotate(geom = 'text', label = paste0("MAJIQ found no cryptics in ", gene_target), x = -Inf, y = Inf, hjust = 0, vjust = 1)
}

}
#patchwork::wrap_plots(violin_plot, nrow = 3, ncol = 2) -> caietto
#caietto

of_interest[c(1:5,7:18)]
#violin_plot<- list()
for (i in of_interest[c(1:5,7:18)]) {
  transcript_bind_plot(i)
  ggsave(filename = paste0("~/Desktop/", i, ".pdf"))
}

transcript_bind_plot("HDGFL2")


