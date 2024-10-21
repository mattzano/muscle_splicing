delta_psi_myoblasts <- read.csv("~/Desktop/controlmyoblast-TDP43KDmyoblast_annotated_junctions.csv")
colnames(delta_psi_myoblasts)[c(12,13)] <- c("control_psi", "disease_psi")
delta_psi_ibm <- read.csv("~/Desktop/controlpatient-IBMpatient_annotated_junctions.csv")
colnames(delta_psi_ibm)[c(12,13)] <- c("control_psi", "disease_psi")
delta_psi_fshd <- read.csv("~/Desktop/Control-FSHD_annotated_junctions.csv")
colnames(delta_psi_fshd)[c(12,13)] <- c("control_psi", "disease_psi")
  
#delta_psi_merged <- rbind(delta_psi_myoblasts, delta_psi_ibm, delta_psi_fshd)
delta_psi_merged <- do.call(rbind, list(myoblast = delta_psi_myoblasts, ibm = delta_psi_ibm, fshd = delta_psi_fshd))

delta_psi_merged <- delta_psi_merged %>%
  rownames_to_column() #%>%

delta_psi_merged$rowname <- gsub('\\..*', '', delta_psi_merged$rowname)
table(delta_psi_merged$rowname)

delta_psi_mut <- delta_psi_merged %>%
  dplyr::filter(gene_name != "CH507-528H12.1") %>%
  dplyr::filter(gene_name != "CH507-513H4.1") %>%
  dplyr::filter(gene_name != "AHNAK") %>%
  mutate(log10_test_stat = -log10(1 - probability_changing)) %>%
  mutate(log10_test_stat = ifelse(is.infinite(log10_test_stat), 4.5, log10_test_stat)) %>%
  mutate(graph_alpha = ifelse(probability_changing > 0.9 & de_novo_junctions == 1, 1, 0.2)) %>%
  mutate(label_junction = case_when(graph_alpha ==1  &
                                      mean_dpsi_per_lsv_junction > 0 ~ gene_name,
                                    T ~ ""))
table(delta_psi_mut$rowname)

delta_psi_mut_filter <- delta_psi_mut %>%
  group_by(paste_into_igv_junction) %>%
  #dplyr::filter(probability_changing > 0.9) %>%
  dplyr::filter(n() >= 2)

delta_psi_mut_filter %>%
  group_by(label_junction) %>%
  dplyr::filter(n() < 4 | label_junction == "") %>%
  ggplot(aes(x = mean_dpsi_per_lsv_junction, y = log10_test_stat, color = rowname, alpha = graph_alpha)) + 
  geom_point(show.legend = F) +
  geom_text_repel(aes(label = label_junction), color = "black", show.legend = F) +
  geom_hline(aes(yintercept = 1)) +
  geom_vline(aes(xintercept = 0)) +
  #scale_color_manual(values = c("#648FFF", "#fe6101")) +
  theme_classic() +
  theme(legend.position = 'top')
