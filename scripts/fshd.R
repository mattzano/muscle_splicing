input_psi <- list.files(path = "~/Desktop/fshd_single_psi/")
datalist_psi <- list()

for (i in input_psi) {
  #input_tsv <- paste(i, "tsv", sep = ".")
  input_psi_tsv <- read.csv(file.path(path = "~/Desktop/fshd_single_psi/", i), header = T, fill = TRUE)
  #input_psi_tsv$Model_Of_Inheritance <- gsub("\\,.*", "", input_psi_tsv$Model_Of_Inheritance)
  input_psi_clean <- input_psi_tsv[,c(1,4,5,8,10,16,17)]
  #input_psi_clean$gene_name <- input_psi_clean$Gene.Symbol
  datalist_psi[[i]] <- input_psi_clean
}

big_psi <- rbindlist(datalist_psi, idcol = TRUE)
big_psi$.id <- factor(big_psi$.id, levels = input_psi)

big_psi$.id <- gsub("\\..*", "", big_psi$.id)
big_psi$group = gsub('[0-9]+', '', big_psi$.id)

big_psi_delta <- big_psi %>%
  group_by(gene_name,paste_into_igv_junction, group) %>%
  summarize(mean_psi = mean(mean_psi_per_lsv_junction)) %>%
  ungroup() %>%
  pivot_wider(names_from = group, values_from = mean_psi) %>%
  dplyr::mutate(delta_psi = FSHD - Control)

big_psi_delta %>%
  dplyr::filter(delta_psi > 0.1 | delta_psi < - 0.1) %>%
  dplyr::mutate(color = ifelse(Control < 0.05, 1,2)) %>%
  dplyr::mutate(gene_label = case_when(color == 1 ~ gene_name, T ~ "")) %>%
  ggplot(aes(x = Control, y = delta_psi, color = color)) +
  geom_text_repel(aes(label = gene_label), max.overlaps = 500) +
  geom_point() + 
  theme_bw()

big_psi_delta_filter <- big_psi_delta %>%
  dplyr::filter(delta_psi > 0.1 & Control < 0.05) %>%
  separate(paste_into_igv_junction, into = c("chr", "start", "end"), remove = F) %>%
  dplyr::mutate(length = as.numeric(end) - as.numeric(start))
