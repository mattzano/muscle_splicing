patient_delta_psi <- read.table("data/controlpatient-IBMpatient_parsed_psi.tsv", sep = ",", header = T)

patient_delta_psi_mut <- patient_delta_psi %>%
  dplyr::filter(gene_name != "CH507-528H12.1") %>%
  dplyr::filter(gene_name != "CH507-513H4.1") %>%
  dplyr::filter(gene_name != "AHNAK") %>%
  dplyr::mutate(log10_test_stat = -log10(1 - probability_changing)) %>%
  dplyr::mutate(graph_alpha = as.character(ifelse(probability_changing > 0.9, 1, 0.2))) %>%
  dplyr::mutate(de_novo_junctions = as.character(de_novo_junctions)) %>%
  dplyr::mutate(label_junction = case_when(graph_alpha == 1 & controlpatient_mean_psi < 0.1 & mean_dpsi_per_lsv_junction > 0 & de_novo_junctions == 1 ~ gene_name, T ~ ""))

patient_delta_psi_mut %>%
  ggplot(aes(x = mean_dpsi_per_lsv_junction, y = log10_test_stat, alpha = graph_alpha, color = de_novo_junctions)) + 
  geom_point(show.legend = F) +
  geom_text_repel(aes(label = label_junction), color = "black", show.legend = F) +
  geom_hline(aes(yintercept = 1)) +
  geom_vline(aes(xintercept = 0)) +
  scale_color_manual(values = c("#648FFF", "#fe6101")) +
  theme_classic()


avg_read_counts <- featureCounts %>%
    mutate(ctrl = rowMeans(dplyr::select(featureCounts, contains("Control")), na.rm = TRUE)) %>%
    mutate(ibm = rowMeans(dplyr::select(featureCounts, contains("IBM")), na.rm = TRUE)) %>%
    mutate(delta = ctrl/ibm) %>%
  dplyr::filter(delta > 1) %>%#dplyr::select(c(20, 23, 21, 22)) %>%
    pivot_longer(cols = c("ctrl","ibm"), 
                 names_to = ".id", values_to = "avg_count") %>%
    group_by(gene_name, .id) %>%
    summarise(avg_count = max(avg_count))
avg_read_counts$.id <- factor(avg_read_counts$.id, levels = c("ctrl","ibm"))
  
avg_read_counts %>%
    dplyr::filter(gene_name %in% patient_junction_list) %>%
    #mutate(label_junction = case_when(.id ==" dox_0187"~ gene_name, T ~ "")) %>%
    ggplot(aes(x = .id, y = avg_count, group = gene_name)) +
    #geom_point(aes(color = gene_name), show.legend = F) +
    geom_line(aes(color = gene_name), show.legend = T) +
    #geom_text_repel(aes(label = label_junction), max.overlaps = Inf, size=4, show.legend = F) +
    scale_y_log10() +
    #xlab("TDP-43 knockdown level") +
    #ylab("Average read count") +
    #scale_x_discrete(labels = c("-", "+")) +
    theme_classic()



patient_junction_list <- patient_junction_list[patient_junction_list != ""]


patient_junction_list <- patient_delta_psi_mut %>%
  dplyr::filter(controlpatient_mean_psi < 0.05 & mean_dpsi_per_lsv_junction > 0.1) %>%
  pull(label_junction)

write.csv(patient_junction_list, "patient_junction_list.csv")

brain_junction_list <- brain_junction_list[,2]

brain_junction_list <- brain_junction_list %>%
  pull()

big_delta_patient <- patient_delta_psi %>%
  dplyr::filter(paste_into_igv_junction %in% brain_junction_list)

big_delta_patient_filtered <- big_delta_patient %>%
  dplyr::filter(controlpatient_mean_psi < 0.05)





patient_delta_fig <- patient_delta_psi %>%
  dplyr::filter(probability_changing > 0.9) %>%
  group_by(gene_name) %>%
  filter(n() < 5) %>%
  dplyr::mutate(color_gene = as.character(ifelse(controlpatient_mean_psi < 0.1 & mean_dpsi_per_lsv_junction > 0.1, 1,
                                    ifelse(controlpatient_mean_psi > 0.9 & mean_dpsi_per_lsv_junction < -0.1, 2, 3)))) %>%
  dplyr::mutate(label_gene = as.character(case_when((controlpatient_mean_psi < 0.1 & mean_dpsi_per_lsv_junction > 0.5) |
                                                    (controlpatient_mean_psi > 0.9 & mean_dpsi_per_lsv_junction < -0.5) ~ gene_name, T ~ "")))



patient_delta_fig %>% 
  ggplot(aes(x = controlpatient_mean_psi, y = mean_dpsi_per_lsv_junction, color = color_gene)) +
  geom_point(show.legend = F) +
  geom_text_repel(aes(label = label_gene), color = "black", max.overlaps = 1000) +
  geom_hline(aes(yintercept = 0.1), linetype = "dotted") +
  geom_hline(aes(yintercept = -0.1), linetype = "dotted") +
  geom_vline(aes(xintercept = 0.1), linetype = "dotted") +
  geom_vline(aes(xintercept = 0.9), linetype = "dotted") +
  theme_classic()

