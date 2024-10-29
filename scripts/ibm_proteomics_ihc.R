##ihc and proteomics

pca_plt <- read.table(here::here("data/Muscle_Proteome_Fratta group/clustering/PCA_labels.tsv"), sep = "\t", header = T)

hd_low <- c("IBM_CTRL_1", "IBM_CTRL_2", "IBM_17", "IBM_CTRL_3", "IBM_CTRL_4", "IBM_16") #6
hd_high <- c("IBM_2", "IBM_5", "IBM_7", "IBM_10", "IBM_9", "IBM_8", "IBM_13", "IBM_3") #8
hd_dk <- c("IBM_14", "IBM_6", "IBM_11", "IBM_15", "IBM_1", "IBM_4", "IBM_12") #7


pca_plt %>% 
  #select(-(c("IBM_1", "IBM_4", "IBM_5", "IBM_7", "6" "11", "12", "14")))
  dplyr::mutate(has_ihc = ifelse(ihc_immune == "", "no", "yes")) %>% 
  dplyr::mutate(HDGFL2_CE = ifelse(sample %in% hd_low, "low", 
                                   ifelse(sample %in% hd_high, "high", "NA"))) %>% 
  ggplot(aes(x = PC1, y = PC2, color = HDGFL2_CE, alpha = has_ihc, shape = Condition)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = sample), color = "black", show.legend = F) +
  scale_color_brewer(palette =  "Set2") +
  scale_alpha_manual(values = c(0.1,1)) +
  theme_bw()

#umap_plt <- read.table(here::here("data/Muscle_Proteome_Fratta group/clustering/UMAP_labels.tsv"), sep = "\t", header = T)

#umap_plt %>% 
#  mutate(alpha_code = ifelse(ihc_immune == "", "low", "high")) %>% 
#  ggplot(aes(x = UMAP1, y = UMAP2, color = ihc_p62, alpha = alpha_code)) +
#  geom_point() +
#  scale_color_brewer(palette =  "Set1") +
#  scale_alpha_manual(values = c(1,0.1)) +
#  theme_bw()

hd_low <- c("IBM_CTRL_1", "IBM_CTRL_2", "IBM_17", "IBM_CTRL_3", "IBM_CTRL_4","IBM_16") #6
hd_high <- c("IBM_2", "IBM_5", "IBM_7", "IBM_10", "IBM_9", "IBM_8", "IBM_13", "IBM_3") #8
hd_dk <- c("IBM_14", "IBM_6", "IBM_11", "IBM_15", "IBM_1", "IBM_4", "IBM_12") #7

main_list <- c("TARDBP", "CD3E", "TRBC2;TRBC1;TRB")
other_list <- c("B2M", "CALR", "ERAP1", "ERAP2", "TAP1", "TAP2", "HLA-A", "HLA-B", "HLA-C") #grepl("THOP", PG.Genes) | grepl("MARCH", PG.Genes) | PG.Genes %in% c("TAP1", "TAP2", "TAPBP") |
#PG.Genes %in% c("TPP1", "TPP2") | PG.Genes == "PDIA3" | PG.Genes == "CALR" | PG.Genes == "B2M" |
#grepl("ERAP", PG.Genes) | PG.Genes == "DERL1" | PG.Genes == "VCP") %>% 
ibm_raw <- read.csv(here::here("data", "Muscle_Proteome_Fratta group/20230807_IBM_protein_Report.csv")) #,
ehi <- unlist(strsplit(names(ibm_raw)[4:24], "..", fixed = T))
ehiii <- unlist(strsplit(ehi[seq(2,42,2)], ".", fixed = T))
names(ibm_raw)[4:24] <- ehiii[seq(1,81,4)]

ibm_raw_long <- ibm_raw %>% 
  filter(PG.Genes %in% main_list | PG.Genes %in% other_list) %>% 
  pivot_longer(cols = c(4:24)) %>% 
  mutate(color_code = ifelse(grepl("CTRL", name), "Control", "IBM")) %>% 
  mutate(color_code_level = #ifelse(name %in% immune_neg, "neg",
                                   ifelse(name %in% hd_low, "Low",
                                          ifelse(name %in% hd_high, "High", "NA"))) %>%
  filter(color_code_level != "NA") %>% 
  mutate(Pg.Genes = factor(PG.Genes, levels = c(main_list,other_list))) %>% 
  mutate(color_code_level = factor(color_code_level, levels = c("Low", "High")),
         PG.Genes = factor(PG.Genes, levels = c("TARDBP", "CD3E", "TRBC2;TRBC1;TRB", "CALR",
                                                "HLA-A", "HLA-B", "HLA-C", "B2M",
                                                "ERAP1", "ERAP2", "TAP1", "TAP2")))
  #mutate(color_code_level = ifelse(name %in% tardbp_data$level, "high", "low")) %>% 
  #mutate(#Genes_f = factor(PG.Genes, levels = c('TARDBP', 'CD3E', 'TRBC2;TRBC1;TRB')),
  #  color_code_f = factor(color_code, levels = c("control", "IBM")))
#ggboxplot(
#  ibm_raw_long, x = "color_code_level", y = "value",
#  add = "mean_se", facet.by = "PG.Genes"
#)
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
ibm_raw_long %>% 
  mutate(alpha_man = ifelse(is.nan(value), "0", "1"),
         value_man = ifelse(is.nan(value), 0, value)) %>% 
  ggplot(aes(x = color_code_level, y = value_man)) +
  geom_boxplot(outliers = F) +
  geom_point(aes(color = color_code), position = position_dodge2(width = 0.5), show.legend = T) +
  ggpubr::stat_compare_means(method = "wilcox", paired = F, symnum.args = symnum.args, 
                     label = "p.signif", label.x.npc = "centre", label.y.npc = "top"
                     ) +
  #stat_kruskal_test(#group.by = 'x.var',
  #  position = position_nudge(y = 100)) +
  scale_alpha_manual(values = c(0.5,1)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  facet_wrap(~PG.Genes, nrow = 3,  scales = "free_y") +
  #scale_color_viridis_d(begin = 0.4, end = 0.9) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "HDGFL2-CE expression", y = "Normalised protein levels", color = "Condition") +
  theme_classic() +
  theme(legend.position = "top")
ggsave("~/Desktop/figure_4e_graded.png", width = 10)

ibm_raw_long_suppl <- ibm_raw %>% 
#grepl("PSMB8", PG.Genes) | grepl("PSMB9", PG.Genes) | grepl("PTPRC", PG.Genes) |
#grepl("TRIM", PG.Genes) |
#grepl("GSK3", PG.Genes) | 
#grepl("GBP2", PG.Genes) | 
#grepl("CYLD", PG.Genes) | 
#grepl("CCR5", PG.Genes) |
#grepl("IRF8", PG.Genes) |
filter(PG.Genes %in% c("TPP1", "TPP2") | PG.Genes == "PDIA3" | PG.Genes == "ACE" | PG.Genes == "TAPBP" |
  grepl("THOP", PG.Genes) | grepl("MARCH", PG.Genes) | PG.Genes == "DERL1" | PG.Genes == "VCP" |#) %>% 
#grepl("METTL", PG.Genes) |
#grepl("CD74", PG.Genes) | ###MHC II filter(grepl("CTP", PG.Genes)) %>% 
grepl("FAS", PG.Genes) |
grepl("HLA", PG.Genes) | 
#grepl("ATXN", PG.Genes) |
grepl("TIA", PG.Genes) |
grepl("GZM", PG.Genes)) %>% 
  filter(PG.Genes != "HLA-A" & PG.Genes != "HLA-B" & PG.Genes != "HLA-C") %>% 
  pivot_longer(cols = c(4:24)) %>% 
  mutate(color_code = ifelse(grepl("CTRL", name), "Control", "IBM")) %>% 
  mutate(color_code_level = #ifelse(name %in% immune_neg, "neg",
           ifelse(name %in% hd_low, "Low",
                  ifelse(name %in% hd_high, "High", "NA"))) %>%
  filter(color_code_level != "NA") %>% 
mutate(Pg.Genes = factor(PG.Genes, levels = c(main_list,other_list))) %>% 
  mutate(color_code_level = factor(color_code_level, levels = c("Low", "High"))) #,
      #   PG.Genes = factor(PG.Genes, levels = c("TARDBP", "CD3E", "TRBC2;TRBC1;TRB", "CALR",
      #                                          "HLA-A", "HLA-B", "HLA-C", "B2M",
      #                                          "ERAP1", "ERAP2", "TAP1", "TAP2")))

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
ibm_raw_long_suppl %>% 
  ggplot(aes(x = color_code_level, y = value)) +
  geom_boxplot(outliers = F) +
  geom_point(aes(color = color_code), position = position_dodge2(width = 0.5), show.legend = T) +
  ggpubr::stat_compare_means(method = "wilcox", paired = F, symnum.args = symnum.args, 
                             label = "p.signif", label.x.npc = "centre", label.y.npc = "top"
  ) +
  #stat_kruskal_test(#group.by = 'x.var',
  #  position = position_nudge(y = 100)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  facet_wrap(~PG.Genes, nrow = 4,  scales = "free_y") +
  #scale_color_viridis_d(begin = 0.4, end = 0.9) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "HDGFL2-CE expression", y = "Normalised protein levels", color = "Condition") +
  theme_classic() +
  theme(legend.position = "top")
ggsave("~/Desktop/suppl_prot.png", width = 12)
