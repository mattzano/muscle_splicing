#read in all splicing tables
#get the list of CE peptides
of_interest <- c("AARS1", "EPB41L4A", "HDGFL2", "PHF2", "PXDN", "SLC24A3", "STMN2", #"ACSF2", 
                 "ACTL6B", "ARHGAP22", "CELSR3", "DNM1", "IGLON5", "SYNJ2", "MYO18A","PTPRZ1", "NECAB2", "XPO4", "ZNF423")


splicing_i3muscles <- read.csv(here::here('data/splicing','i3muscles_Control-TDP43KD_annotated_junctions.csv')) %>% 
  mutate(sig = ifelse(control_mean_psi < 0.05 & mean_dpsi_per_lsv_junction > 0.1, "yes", "no"),
         gene_of_interest = ifelse(gene_name %in% of_interest, "yes", "no"),
         experiment = "muscle")

#splicing_i3MNs <- read.csv(here::here('data/splicing', "MN_Control-TDP43KD_annotated_junctions.csv")) %>% 
#  mutate(sig = ifelse(control_mean_psi < 0.05 & mean_dpsi_per_lsv_junction > 0.1, "yes", "no"),
#         gene_of_interest = ifelse(gene_name %in% of_interest, "yes", "no")) 
#sig_junctions_i3MNs <- splicing_i3MNs %>% 
#  filter(control_mean_psi < 0.05, 
#         mean_dpsi_per_lsv_junction > 0.1)

splicing_i3corticals <- read.csv(here::here('data/splicing', "Cortical_Control-TDP43KD_annotated_junctions.csv")) %>% 
  mutate(sig = ifelse(control_mean_psi < 0.05 & mean_dpsi_per_lsv_junction > 0.1, "yes", "no"),
         gene_of_interest = ifelse(gene_name %in% of_interest, "yes", "no"),
         experiment = "cortex") 
sig_junctions_i3corticals <- splicing_i3corticals %>% 
  filter(control_mean_psi < 0.05, 
         mean_dpsi_per_lsv_junction > 0.1)

#splicing_i3astros <- read.csv(here::here('data/splicing', "Astro_Control-TDP43KD_annotated_junctions.csv")) %>% 
#  mutate(sig = ifelse(control_mean_psi < 0.05 & mean_dpsi_per_lsv_junction > 0.1, "yes", "no"),
#         gene_of_interest = ifelse(gene_name %in% of_interest, "yes", "no")) 
#sig_junctions_i3astros <- splicing_i3astros %>% 
#  filter(control_mean_psi < 0.05, 
#         mean_dpsi_per_lsv_junction > 0.1)

#see where they do appear and how much overlap
splicing_i3muscles %>% 
  mutate(probability_changing = ifelse(probability_changing == 1, 0.99995,probability_changing)) %>% 
  #mutate(color_code = ifelse(paste_into_igv_junction %in% sig_junctions_i3muscles$paste_into_igv_junction, "significant","not significant")) %>%
  mutate(label_code = ifelse(sig == "yes" & gene_name %in% of_interest, gene_name, NA)) %>% 
  ggplot(aes(x = mean_dpsi_per_lsv_junction, y = -log10(1-probability_changing))) +
  geom_point(aes(color = sig)) +
  ggrepel::geom_text_repel(aes(label = label_code), max.overlaps = Inf, min.segment.length = 0.2) +
  #scale_y_continuous(trans = "log10") +
  theme_classic()
#ggsave("~/Desktop/extended_9a.png")

#splicing_i3MNs %>% 
#  mutate(probability_changing = ifelse(probability_changing == 1, 0.99995,probability_changing)) %>% 
#  #mutate(color_code = ifelse(paste_into_igv_junction %in% sig_junctions_i3muscles$paste_into_igv_junction, "significant","not significant")) %>%
#  mutate(label_code = ifelse(sig == "yes" & gene_name %in% of_interest, gene_name, NA)) %>% 
#  ggplot(aes(x = mean_dpsi_per_lsv_junction, y = -log10(1-probability_changing))) +
#  geom_point(aes(color = sig)) +
#  ggrepel::geom_text_repel(aes(label = label_code), max.overlaps = Inf, min.segment.length = 0.2) +
#  #scale_y_continuous(trans = "log10") +
#  theme_classic()
#ggsave("~/Desktop/extended_9b.png")

splicing_i3corticals %>% 
  mutate(probability_changing = ifelse(probability_changing == 1, 0.99995,probability_changing)) %>% 
  #mutate(color_code = ifelse(paste_into_igv_junction %in% sig_junctions_i3muscles$paste_into_igv_junction, "significant","not significant")) %>%
  mutate(label_code = ifelse(sig == "yes" & gene_name %in% of_interest, gene_name, NA)) %>% 
  ggplot(aes(x = mean_dpsi_per_lsv_junction, y = -log10(1-probability_changing))) +
  geom_point(aes(color = sig)) +
  ggrepel::geom_text_repel(aes(label = label_code), max.overlaps = Inf, min.segment.length = 0.2) +
  #scale_y_continuous(trans = "log10") +
  theme_classic()
#ggsave("~/Desktop/extended_9c.png")

#splicing_i3astros %>% 
#  mutate(probability_changing = ifelse(probability_changing == 1, 0.99995,probability_changing)) %>% 
  #mutate(color_code = ifelse(paste_into_igv_junction %in% sig_junctions_i3muscles$paste_into_igv_junction, "significant","not significant")) %>%
#  mutate(label_code = ifelse(sig == "yes" & gene_name %in% of_interest, gene_name, NA)) %>% 
#  ggplot(aes(x = mean_dpsi_per_lsv_junction, y = -log10(1-probability_changing))) +
#  geom_point(aes(color = sig)) +
#  ggrepel::geom_text_repel(aes(label = label_code), max.overlaps = Inf, min.segment.length = 0.2) +
  #scale_y_continuous(trans = "log10") +
#  theme_classic()

list_ce_i3corticals <- splicing_i3corticals %>% 
  filter(sig == "yes" & gene_of_interest == "yes") %>% 
  pull(gene_name) %>% unique() %>% sort()
#list_ce_i3lmns <- splicing_i3MNs %>% 
#  filter(sig == "yes" & gene_of_interest == "yes") %>% 
#  pull(gene_name) %>% unique() %>% sort()
list_ce_i3muscles <- splicing_i3muscles %>% 
  filter(sig == "yes" & gene_of_interest == "yes") %>% 
  pull(gene_name) %>% unique() %>% sort()

nauc_ce_labs <- nauc_ce %>% 
  mutate(gene_symbol = ifelse(gene_symbol == "AARS", "AARS1", gene_symbol)) %>% 
  group_by(name) %>% 
  mutate(splicing_i3 = case_when(name == "Frontal Cortex" & gene_symbol %in% list_ce_i3corticals ~ "yes",
                                 #name == "Spinal Cord" & gene_symbol %in% list_ce_i3lmns ~ "yes",
                                 name == "Skeletal Muscle" & gene_symbol %in% list_ce_i3muscles ~ "yes")) 

nauc_ce_labs %>% 
  ggplot(aes(y = gene_symbol, x = reorder(name,-value))) +
  geom_tile(aes(fill = value), show.legend = F) +
  geom_point(aes(color = splicing_i3), shape = 16, size = 3, show.legend = F) +
  scale_color_viridis_d(begin = 0.5) +
  #scale_fill_viridis_c() +
  scale_fill_gradient(low = "white", high = "black") +
  ylab("Gene") + xlab("Tissue") + #labs(fill = "gene expression") +
  scale_y_discrete(limit = rev) +
  theme_classic()
ggsave("~/Desktop/fig_xx.png")

prot_i3 <- read.csv(here::here("data/20230114_MCP_Report_protein_total 2.csv")) %>% 
  filter(PG.Genes %in% of_interest | PG.Genes == "TARDBP") %>%
  pivot_longer(cols = c(4:51)) %>% 
  mutate(value = as.numeric(ifelse(value == "Filtered", NA, value))) %>% 
  mutate(cell_type = case_when(grepl("Astro", name) ~ "Astro",
                               grepl("Motor", name) ~ "Motor_Neuron",
                               grepl("NgN2", name) ~ "Cortical_Neuron",
                               grepl("Muscle", name) ~ "Skeletal_Muscle")) %>% 
  mutate(color_code = ifelse(grepl("ctrl", name), "Control", "TDP43KD")) %>% 
  filter(!grepl("Astro", name))
prot_i3 %>% 
  ggplot(aes(x = color_code, y = value, group = color_code)) +
  geom_boxplot(show.legend = F) +
  geom_point(aes(color = color_code), position = position_dodge2(width = 0.2), show.legend = F) +
  facet_wrap(PG.Genes~cell_type, scales = "free", ncol = 9) +
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_color_viridis_d() +
  labs(x = "Condition", y = "Normalised protein levels", color = "TDP-43 expression") +
  theme_classic() #+
  #theme(legend.position = "top")

  #geom_point() + # (aes(fill = color_code),show.legend = T) +
  #scale_fill_brewer(palette = "Set1") +
  #scale_y_continuous(labels = c("","","","","")) +
  #xlab("") + ylab("TARDBP normalised protein levels") +
  #labs(fill = "Condition") +
  #theme_bw() +
  #legend(x = 1,y =350000) +
  #theme(#axis.text.x = element_blank(), #text(angle = 90),
        #axis.ticks.y = element_blank(),
  #      legend.position = "bottom")


prot_i3_nauc <- prot_i3 %>% 
  filter(color_code == "Control" & PG.Genes != "TARDBP") %>% 
  group_by(PG.Genes, cell_type) %>% 
  summarise(avg = mean(value)) %>% 
  group_by(cell_type) %>% 
  mutate(splicing_i3 = case_when(cell_type == "Cortical_Neuron" & PG.Genes %in% list_ce_i3corticals ~ "yes",
                                 cell_type == "Motor_Neuron" & PG.Genes %in% list_ce_i3lmns ~ "yes",
                                 cell_type == "Skeletal_Muscle" & PG.Genes %in% list_ce_i3muscles ~ "yes")) 
#prot_i3$cell_type

prot_i3_nauc %>% 
  ggplot(aes(y = PG.Genes, x = cell_type)) + #reorder(cell_type,avg))) +
  geom_tile(aes(fill = avg), show.legend = F) +
  geom_point(aes(color = splicing_i3), shape = 16, size = 3, show.legend = T) +
  scale_color_viridis_d(begin = 0.4) +
  #scale_fill_viridis_c() +
  scale_fill_gradient(low = "white", high = "black") +
  ylab("Gene") + xlab("Tissue") + #labs(fill = "gene expression") +
  scale_y_discrete(limit = rev) +
  theme_classic()




