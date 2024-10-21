####new muscles

###deseq2
salmon_deseq2(sample_dir =  here::here('data', 'salmon'),
              metadata_dir =  here::here('data','metadata_muscle_salmon.csv'),
              experiment_dir = "nhnn")

normed_counts_long <- normed_counts %>% 
  #filter(symbol %in% of_interest) %>% 
  filter(symbol %in% c("TRAC", "TRBC1", "TRBC2", "CD3E") | 
           symbol %in% c("HLA-A", "HLA-B", "HLA-C") |
           symbol %in% c(#"TAP1", "TAP2", "ERAP1", "ERAP2", "CALR", 
             "B2M")) %>%  #grepl("HLA", symbol)) %>% #symbol == "TARDBP") %>% 
  pivot_longer(cols = c(2:17)) %>% 
  mutate(condition = ifelse(grepl("IBM", name) | grepl("PT", name), "IBM", 
                            ifelse(grepl("CTR", name) | grepl("Ctrl", name) | grepl("Control", name), "Control",
                                   ifelse(grepl("C9", name) | grepl("FUS", name) | grepl("SOD1", name) | grepl("TDP43", name), "ALS", "Other muscle diseases")))) #%>% 
#mutate(condition = ifelse(grepl("IBM", name) | grepl("PT", name), "IBM", "Control"))
#write.table(normed_counts_long, "~/Desktop/normed_long.csv", sep = ",", row.names = F, quote = F)
#normed_clean <- read.csv("~/Desktop/normed_long.csv")

#normed_clean %>% 
#  filter(symbol %in% c("TRAC", "TRBC1", "TRBC2", "CD3E") | 
#           symbol %in% c("HLA-A", "HLA-B", "HLA-C")) %>% 
#    ggplot(aes(x = name_clean, y = value, fill = condition)) +
#  geom_col() +
#  facet_wrap(~symbol, nrow = 1, scales = "free_y") +
#  scale_fill_brewer(palette = "Set2") +
#  theme_classic() +
#  labs(x = "", y = "Normalised read counts", fill = "Condition") +
#  theme(axis.text.x = element_text(angle = 90))

#RColorBrewer::brewer.pal(n = 4, "Set2")
normed_counts_long %>% 
  group_by(condition,symbol) %>% 
  summarise(mean_val = mean(value)) %>% 
  #mutate(symbol_f = factor(symbol, levels=c('HLA-A','HLA-B','HLA-C', "B2M",'TRAC', 'TRBC1', 'TRBC2', 'CD3E' #"TAP1", "TAP2", "ERAP1", "ERAP2", "CALR", 
  #))) %>% 
  ggplot(aes(x = condition,#reorder(name,value), 
               y = mean_val, fill = condition)) +
  geom_col() +
  facet_wrap(~symbol, nrow = 2, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  labs(x = "", y = "Normalised read counts", fill = "Condition") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")


normed_counts_long %>% 
  #group_by(condition,symbol) %>% 
  #summarise(mean_val = mean(value)) %>% 
  mutate(symbol_f = factor(symbol, levels=c('HLA-A','HLA-B','HLA-C', "B2M",'TRAC', 'TRBC1', 'TRBC2', 'CD3E' #"TAP1", "TAP2", "ERAP1", "ERAP2", "CALR", 
  ))) %>% 
  ggplot(aes(x = reorder(name,value), 
             y = value, fill = condition)) +
  geom_col() +
  facet_wrap(~symbol_f, nrow = 2, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  labs(x = "", y = "Normalised read counts", fill = "Condition") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")

of_interest <- c("AARS1", "EPB41L4A", "HDGFL2", "PHF2", "PXDN", "SLC24A3", "STMN2", #"ACSF2", 
                 "ACTL6B", "ARHGAP22", "CELSR3", "DNM1", "IGLON5", "SYNJ2", "MYO18A",
                 "PTPRZ1", "NECAB2", "XPO4", "ZNF423")

splicing_bedparse_invivo <- read.table("data/significant_aggregated.clean.annotated_nhnn.bed", sep = "\t") %>%  #here::here('data', "significant_aggregated.clean.annotated.bed"), sep = "\t") %>% 
  filter(V7 %in% of_interest) %>% 
  filter(!grepl("SBMA", V4)) %>%  #filter(V7 != "MYO18A") %>% 
  mutate(V4 = gsub("\\..*","", V4)) %>% 
  mutate(condition = ifelse(grepl("IBM", V4) | grepl("PT", V4), "IBM", 
                            ifelse(grepl("CTR", V4) | grepl("Ctrl", V4) | grepl("Control", V4), "Control",
                                   ifelse(grepl("C9", V4) | grepl("FUS", V4) | grepl("SOD1", V4) | grepl("TDP43", V4), "ALS", "Other muscle diseases")))) %>% 
  group_by(V4,V7) %>% 
  arrange(desc(V5)) %>% 
  distinct(V7, .keep_all = T) %>% 
  ungroup() #%>% 
  #write.table(splicing_bedparse_invivo, "~/Desktop/bedtools.csv", sep = ",", row.names = F, quote = F)

metadata_dir <- read.csv("data/metadata_muscle_salmon.csv")

metadata_data <- metadata_dir %>% 
  dplyr::filter(experiment == "nhnn") %>% 
  mutate(V4 = sample_name,
         condition = ifelse(condition == 0, "Control", "IBM"))

splicing_bedparse_invivo_normalised <- splicing_bedparse_invivo %>% 
  left_join(metadata_data[,c(7,8)])

splicing_bedparse_invivo_normalised %>% 
  filter(condition %in% c("Control", "IBM")) %>%
  #dplyr::filter(V7 != "MYO18A" & V7 != "XPO4") %>% 
  dplyr::filter(V7 == "HDGFL2") %>% 
  ggplot(aes(x = reorder(V4,V5/paired.reads,sum), 
             y = V5/paired.reads, fill = V7)) +
  geom_col() +
  scale_fill_viridis_d(direction = -1) +
  theme_classic() +
  facet_wrap(~condition, scales = "free_x") +
  labs(x ="", y="Reads covering cryptic junction", fill = "Gene") +
  theme(axis.text.x = element_text(angle = 90))


violin_plot <- list()
for (i in c('HLA-A','HLA-B','HLA-C',"B2M",'TRAC', 'TRBC1', 'TRBC2', 'CD3E' 
            #"TAP1", "TAP2", "ERAP1", "ERAP2", "CALR", "B2M", "TARDBP"
)) {
  normed_clean_merged <- splicing_bedparse_invivo_normalised %>% 
    #dplyr::filter(#V4 != "IBM_4" & V4 != "IBM_1" &
    #
    #V4 != "IBM_7" &V4 != "IBM_5") %>% 
    filter(V7 != "MYO18A") %>% 
    mutate(name = V4) %>% 
    #filter(!is.na(count_max)) %>% 
    group_by(name) %>% 
    summarise(sum(V5)/paired.reads) %>% 
   
    left_join(normed_counts_long) %>% 
    filter(symbol == i)
  violin_plot[[i]] <- normed_clean_merged %>%
    # %>% # & name_clean != "IBM_Trieste_4") %>% 
    ggplot(aes(x = `sum(V5)/paired.reads`, 
               y = value)) +
    geom_point(aes(color = condition), alpha = 0.7, size = 3, show.legend = F) +
    geom_smooth(method = 'lm', se =F, color = "black") +
    scale_color_brewer(palette = "Set2") +
    geom_label(aes(label = name)) +
    ggpubr::stat_cor() +
    labs(x = "Reads covering CE junctions", y = paste0(i," normalised gene expression"), color = "Condition") +
    theme_bw() +
    theme(legend.position = "top")
}

patchwork::wrap_plots(violin_plot, nrow = 2) -> caio
caio


# Prepare the data (replace with your actual data)
data <- normed_counts[,c(1:17)] #%>% 
  select(-(c("IBM_1", "IBM_4", "IBM_5", "IBM_7")))
    #"Ctrl_3", "Ctrl_1", "Ctrl_4","Ctrl_2", "IBM_5", "IBM_7", "IBM_17"))) #,
    #"IBM_1", "IBM_4", "IBM_5", "IBM_7")))

# Remove the gene column and convert to numeric matrix
data_matrix <- as.matrix(data[,-1])
row.names(data_matrix) <- data$Gene

# Perform PCA
pca <- prcomp(t(data_matrix), scale = TRUE)

# Prepare PCA data for ggplot
pca_data <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Sample = colnames(data_matrix)
)

# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  ggtitle("PCA Plot of Gene Expression Data") +
  xlab(paste("PC1 (", round(100 * summary(pca)$importance[2, 1], 2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(100 * summary(pca)$importance[2, 2], 2), "%)", sep = "")) +
  theme_minimal()

####pca loadings!!


###go ontology on RNA seq



