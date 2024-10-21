####ibm data from buratti

##start from de - pca/umap/kmeans
library(tidyverse)
normed_long_both %>%
  dplyr::filter(symbol == "TARDBP" | symbol == "UPF1" | symbol == "CYFIP2" | symbol == "INSR") %>%
  ggplot(aes(x = group, y = value, color = source)) +
  geom_boxplot(alpha = 0.1, show.legend = F) +
  geom_point(position = position_dodge2(width = 1), show.legend = F) +
  #stat_pvalue_manual(stattest, hide.ns = T) + #[c(1,2),]) +
  #scale_y_log10() +
  facet_wrap(facets = vars(symbol), scales = "free_y") +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

normed_counts %>% 
  filter(symbol == "TARDBP") %>% 
  pivot_longer(cols = c(2:13)) %>% 
  ggplot(aes(x = name, y = value)) +
  geom_col()

###load splicing tables and merge


##load single files psi table and do pca / umap
library(data.table)
datalistona <- list()
samples <- c(
             "CTR_01", #"CTR_02", 
             "CTR_03", "CTR_04", "CTR_05", "CTR_06",
             #"PT_01", "PT_02", 
             "PT_03", "PT_04", "PT_05", #"PT_06", #"PT_07", 
             "PT_08", "PT_09", "PT_10") #"PT_11", #"PT_12",)
for (j in samples) {
  input_splicings <- fread(paste0("/Users/matteozanovello/Documents/phd/research_lines/muscle/muscle_splicing/data/buratti_splicing/", 
                                 j, ".Aligned.sorted.out_parsed.csv"))
  #names(input_splicing)[12] <- "base_mean_psi"
  #names(input_splicing)[13] <- "contrast_mean_psi"
  assign(paste0("single_splicing_",j), input_splicings)
  
  #a <- splicing_dots_tables_function(input_splicing) + annotate("text", x=0.9, y=0.1, label=i)
  #plot(a)
  
  datalistona[[j]] <- input_splicings
}
big_single <- rbindlist(datalistona, idcol = TRUE)
big_single$.id <- factor(big_single$.id, levels = samples)

poldippo <- big_single %>%
  filter(gene_name == "POLDIP3" &
           paste_into_igv_junction %in% c("chr22:42599793-42602770",
                                          "chr22:42599793-42601970",
                                          "chr22:42602056-42602770")) %>%
  group_by(.id) %>% 
  distinct(paste_into_igv_junction, .keep_all = T) %>% 
  mutate(relative_inclusion = mean_psi_per_lsv_junction[paste_into_igv_junction == "chr22:42599793-42602770"] /
           (mean_psi_per_lsv_junction[paste_into_igv_junction == "chr22:42599793-42601970"] + 
           mean_psi_per_lsv_junction[paste_into_igv_junction == "chr22:42602056-42602770"])) %>% 
  separate(.id, "_", into = c("group", "number"), remove = F) 

tdp_levels_rna <- data.frame(.id = samples,
           #c(CTR01,	       CTR03,	      CTR04,	      CTR05,	     CTR06,	    	
          #   PT03,     	  PT04,     	  PT05,	        PT08,	       PT09,        PT10),
            levels = c(5801.649316,	5000.666671,	5673.653342,	5550.778294, 5175.742603,
              5637.583819,	6053.868158,	5934.513623,	4717.264803, 5318.440996,	6458.295548))
tdp_levels_rna %>% 
  separate(.id, "_", into = c("group", "number"), remove = F) %>% 
  ggplot(aes(x = group, y = levels)) +
  geom_point(aes(color = .id)) +
  theme_bw()

poldippo %>%
  mutate(.id = factor(.id, levels = samples)) %>% 
  distinct(.id, .keep_all = T) %>%  
  left_join(tdp_levels_rna) %>% 
  ggplot(aes(x = levels, y = relative_inclusion)) +
  #geom_violin() +
  geom_point(aes(color = group), position = position_dodge2(width = 0.7)) +
  geom_smooth(method = "lm", se = F) +
  ggrepel::geom_text_repel(aes(label = .id)) +
  #facet_wrap(facets = vars(paste_into_igv_junction)) +
  theme_bw()


stammino <- big_single %>%
  filter(gene_name == "GPSM2") %>% 
  group_by(.id) %>% 
  distinct(paste_into_igv_junction, .keep_all = T) %>% 
  separate(.id, "_", into = c("group", "number"), remove = F) %>% 
  group_by(paste_into_igv_junction,group) %>% 
  mutate(grouped_psi = mean(mean_psi_per_lsv_junction)) #%>% 
  #ungroup() %>% 
  #filter(number == "10" | group == "CTR") 

GPSM2 chr1:108877228-108896864 <- this should be only significant
chr1:109438939-109439053 hg19

ACSF2 chr17:48539181-48539246 hg19 
chr17:50464599-50471028

HDGFRP2 chr19:4492012-4492149
ZFP91 chr11:58384466-58384527
position_lab <- position_jitter(seed = 1, height = 0, width = 0.3)

smol_delta <- big_delta %>% 
  filter(base_mean_psi < 0.05 & probability_changing > 0.8 & .id != "IBMa-IBMb")

big_single %>% 
  filter(paste_into_igv_junction == "chr1:108877228-108896864") %>% 
  separate(.id, "_", into = c("group", "number"), remove = F) %>% 
  #mutate(.id = factor(.id, levels = samples)) %>% 
  #distinct(.id, .keep_all = T) %>%  
  #left_join(tdp_levels_rna) %>% 
  ggplot(aes(x = group, y = mean_psi_per_lsv_junction)) +
  geom_violin() +
  geom_point(aes(color = group#,size = levels
                 ), position = position_lab, show.legend = F) +
  geom_smooth(method = "lm", se = F) +
  #ggrepel::geom_text_repel(aes(label = .id), position = position_lab) +
  #ggpubr::stat_pvalue_manual()
  #facet_wrap(facets = vars(paste_into_igv_junction)) +s
  scale_color_brewer(palette = "Set1") +
  xlab("") + ylab("Percent spliced in") +
  labs(size = "TARDBP levels",
       color = "Condition",
       title = "GPSM2 cryptic exon inclusion") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#chr19:7159645-7163032 INSR?

#chr12:101663560-101665057 ####


#only in myoblasts?
#xxx RANBP1 chr22:20110103-20110220
#xxx PFKP chr10:3141749-3142011


big_single_id <- big_single %>%
  select(.id, mean_psi_per_lsv_junction, lsv_junc) %>%
  pivot_wider(names_from = .id, values_from = mean_psi_per_lsv_junction) %>%
  drop_na() #%>%
 

big_single_id_pca <- big_single_id %>%
  mutate(mean_PT = rowMeans(big_single_id[2:7], na.rm = T),
         mean_CTR = rowMeans(big_single_id[8:12], na.rm = T)) %>%
  filter(mean_PT - mean_CTR > 0.1) %>%
  select(c(1:12)) %>%
  column_to_rownames("lsv_junc")


#penguins_meta <- big_single_id_pca %>%
#  rownames_to_column("lsv_junc") #%>%

#library(umap)
set.seed(420)
big_single_id_pca
umap_fit <- umap(big_single_id_pca[c(1:1000),], n_neighbors = 100, min_dist = 0.05)
#umap_fit <- dbscan(umap_fitt[,c(2:5)], eps = 0.05, minPts = 5)
#umap_fit$layout

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  dplyr::rename(UMAP1="V1",
                UMAP2="V2") %>%
  dplyr::mutate(ID = dplyr::row_number()) %>%
  dplyr::inner_join(penguins_meta, by="ID") #%>%
#dplyr::filter(color_gene_name != "none")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2)) +
             #color = category,#)) +
             #size = delta_psi)) +
  geom_point() +
  #geom_text_repel(aes(label = gene_name),
  #                max.overlaps = Inf,
  #                show.legend = F) +
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot") +
  scale_color_brewer(palette = "Set1") +
  theme_bw()





#base <- "noDox"
contrast <- c("Control-IBM", "Control-IBMa", "Control-IBMb", "IBMa-IBMb")
datalista <- list()
for (i in contrast) {
  input_splicing <- fread(paste0("/Users/matteozanovello/Documents/phd/research_lines/muscle/muscle_splicing/data/buratti/", 
                                 i, "_annotated_junctions.csv"))
  names(input_splicing)[12] <- "base_mean_psi"
  names(input_splicing)[13] <- "contrast_mean_psi"
  assign(paste0("DS_",i), input_splicing)
  
  #a <- splicing_dots_tables_function(input_splicing) + annotate("text", x=0.9, y=0.1, label=i)
  #plot(a)
  
  datalista[[i]] <- input_splicing
}
big_delta <- rbindlist(datalista, idcol = TRUE)
big_delta$.id <- factor(big_delta$.id, levels = contrast)

big_delta %>%
  mutate(color_label = ifelse(base_mean_psi < 0.05, "cryptic", "non cryptic")) %>%
  filter(probability_changing > 0.9 & mean_dpsi_per_lsv_junction > 0.1) %>%
  ggplot(aes(x = .id, y = mean_dpsi_per_lsv_junction)) +
  geom_boxplot() +
  geom_point(aes(color = color_label),
             alpha = 0.2,
             position = position_dodge2(width = 0.5)) +
  #facet_wrap(facets = vars(junc_cat)) +
  theme_classic()

big_delta_ctrl_ibm <- big_delta %>%
  filter(.id == "Control-IBM") %>%
  #filter(mean_dpsi_per_lsv_junction > 0.1) %>%
  mutate(color_label = ifelse(base_mean_psi < 0.05, "cryptic", "non cryptic")) %>%
  mutate(sig_label = ifelse(probability_changing > 0.9, "sig", "not sig"))

big_delta_ctrl_ibm %>%
  filter(mean_dpsi_per_lsv_junction > 0.1) %>%
  ggplot(aes(x = base_mean_psi, y = mean_dpsi_per_lsv_junction)) +
  #geom_boxplot() +
  geom_point(aes(color = junc_cat,
             alpha = sig_label)) +
             #position = position_dodge2(width = 0.5)) +
  geom_vline(xintercept = 0.05) +
  scale_alpha_manual(values = c(0.2,1)) +
  #facet_wrap(facets = vars(junc_cat)) +
  facet_wrap(facets = vars(junc_cat)) +
  theme_classic()

#volcano
big_delta_ctrl_ibm %>%
  ggplot(aes(x = (mean_dpsi_per_lsv_junction), y = -log10(1-probability_changing))) +
  #geom_boxplot() +
  geom_point(aes(color = color_label,
                 alpha = sig_label)) +
  #position = position_dodge2(width = 0.5)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.1)) +
  scale_alpha_manual(values = c(0.2,1)) +
  theme_classic()

###assess differences in splicing between ibm and controls


###assess differences in splicing between 2 ibm groups