library(tidyverse)
library(reshape2)
library(ggpubr)
library(rstatix)


###t-test, BH adjusted
ibm_de <- read.table(here::here("data", "Muscle_Proteome_Fratta group/differential_intensity/IBM-vs-IBM_CTRL.tsv"),
                     sep = "\t", header = T)

ibm_de %>% 
  mutate(labelling = ifelse(Genes == "TARDBP" | (abs(log2FC) > 5 & `p.adj` < 0.05), Genes, NA)) %>% 
  mutate(alpha_code = ifelse(`p.adj` < 0.05 & abs(log2FC) > 3, "A", "B")) %>% 
  mutate(color_code = ifelse(`p.adj` < 0.05 & alpha_code == "A", "A", "B")) %>% 
  ggplot(aes(x = log2FC, y = -log10(`p.adj`))) +
  geom_point(aes(alpha = alpha_code, color = color_code), show.legend = F) +
  geom_hline(yintercept = -log10(0.05)) +
  ggrepel::geom_text_repel(aes(label = labelling), 
                           min.segment.length = 0.1,
                           max.overlaps = Inf) +
  scale_alpha_manual(values = c(1,0.1)) +
  scale_color_manual(values = c("red", "black")) +
  theme_classic()

ibm_raw <- read.csv(here::here("data", "Muscle_Proteome_Fratta group/20230807_IBM_protein_Report.csv")) #,
ehi <- unlist(strsplit(names(ibm_raw)[4:24], "..", fixed = T))
ehiii <- unlist(strsplit(ehi[seq(2,42,2)], ".", fixed = T))
names(ibm_raw)[4:24] <- ehiii[seq(1,81,4)]


#sum(ibm_raw[,c(4:24)])

df_long <- melt(ibm_raw[,c(4:24)], variable.name = "Sample", value.name = "Expression")
ggplot(df_long, aes(x = Sample, y = Expression)) +
  geom_violin(aes(fill = Sample), show.legend = F) +
  geom_boxplot(alpha = 0.1) +
  labs(title = "Distribution of Expression Levels by Sample",
       y = "Expression Level",
       x = "Sample") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

df_long %>% 
  ggplot(aes(x = log10(Expression))) +
  geom_histogram(binwidth = 0.5) +
  facet_wrap(~ Sample)

ibm_raw %>% 
  filter(PG.Genes %in% c("TARDBP")) %>% 
  pivot_longer(cols = c(4:24)) %>% 
  mutate(color_code = ifelse(grepl("CTRL", name), "control", "IBM")) %>% 
  ggplot(aes(x = reorder(name, value), y = value)) +
  geom_col(aes(fill = color_code),show.legend = T) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = c("","","","","")) +
  xlab("") + ylab("TARDBP normalised protein levels") +
  labs(fill = "Condition") +
  theme_bw() +
  #legend(x = 1,y =350000) +
  theme(axis.text.x = element_text(angle = 90),
        axis.ticks.y = element_blank(),
        legend.position = "bottom")
#ggsave("~/Desktop/tdp43_levels_proteomics.pdf")



##add pvalue manually to the normed count table

main_list <- c("TARDBP", "CD3E", "TRBC2;TRBC1;TRB")

ibm_raw_long <- ibm_raw %>% 
  filter(#PG.Genes %in% main_list) %>% # 
           #grepl("PSMB8", PG.Genes) | grepl("PSMB9", PG.Genes) | grepl("PTPRC", PG.Genes) |
           #grepl("TRIM", PG.Genes) |
           #grepl("GSK3", PG.Genes) | 
           #grepl("GBP2", PG.Genes) | 
           #grepl("CYLD", PG.Genes) | 
           #grepl("CCR5", PG.Genes) |
           #grepl("IRF8", PG.Genes) |
           #grepl("THOP", PG.Genes) | grepl("MARCH", PG.Genes) | PG.Genes %in% c("TAP1", "TAP2", "TAPBP") |
           #PG.Genes %in% c("TPP1", "TPP2") | PG.Genes == "PDIA3" | PG.Genes == "CALR" | PG.Genes == "B2M" |
           #grepl("ERAP", PG.Genes) | PG.Genes == "DERL1" | PG.Genes == "VCP") %>% 
           #grepl("METTL", PG.Genes) |
           #grepl("CD74", PG.Genes) | ###MHC II filter(grepl("CTP", PG.Genes)) %>% 
           #grepl("FAS", PG.Genes)) %>% 
           #grepl("HLA", PG.Genes)) %>% 
           #grepl("ATXN", PG.Genes)) %>% 
           #grepl("TIA", PG.Genes)) %>% 
  grepl("GZM", PG.Genes)) %>% 
  #PG.Genes == "ACE") %>% 
  pivot_longer(cols = c(4:24)) %>% 
  mutate(color_code = ifelse(grepl("CTRL", name), "control", "IBM")) %>% 
  mutate(color_code_level = ifelse(name %in% tardbp_data$level, "high", "low")) %>% 
  mutate(#Genes_f = factor(PG.Genes, levels = c('TARDBP', 'CD3E', 'TRBC2;TRBC1;TRB')),
    color_code_f = factor(color_code, levels = c("control", "IBM")))

ibm_raw_long  %>% distinct(PG.Genes) %>% pull()
ibm_raw_long$value[is.nan(ibm_raw_long$value)]<-NA

tardbp_data <- ibm_raw_long %>% 
  filter(PG.Genes == "TARDBP") %>% 
  arrange(desc(value)) %>% 
  mutate(level2 = c(rep("high", 10), rep("low",11))) %>% 
  mutate(level = ifelse(level2 == "high", name, NA))


stat.test.test <- ibm_de %>% 
  filter(Genes %in% ibm_raw_long$PG.Genes) %>% 
  #filter(PG.Genes == "ACE") %>% 
  mutate(group1 = "IBM", group2 = "control", PG.Genes = Genes, term = "PG.Genes", Genes_f = Genes,
         null.value = 0) %>% 
  arrange(PG.Genes)

stat.test.test = add_significance(stat.test.test)
stat.test.test

cel <- ibm_raw_long %>% 
  filter(PG.Genes %in% stat.test.test$Genes) %>% 
  group_by(PG.Genes) %>% 
  summarise(max_val = max(value, na.rm = TRUE))
cel

#stat.test.test <- stat.test.test %>% 
#  add_xy_position(x = "color_code", fun = "mean_se", scales = "free")
ggboxplot(
  ibm_raw_long, x = "color_code_f", y = "value",
  add = "mean_se", facet.by = "PG.Genes"
) +
  geom_point(aes(color = color_code_level), position = position_dodge2(width = 0.2), show.legend = T) +
  facet_wrap(~as.factor(PG.Genes), scales = "free") +
  #stat_pvalue_manual(stat.test.test, label = "p.adj.signif", 
  #                   tip.length = 0, y.position = 1.2*cel$max_val) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_color_viridis_d() +
  labs(x = "Condition", y = "Normalised protein levels", color = "TDP-43 expression") +
  theme_classic() +
  theme(legend.position = "top")
#ggsave("~/Desktop/suppl_5b.png", width = 9)


