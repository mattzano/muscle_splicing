###trust results
library(data.table)
library(tidyverse)
#list_files <- list.files(path = here::here('data', 'trust4/trust4_mayo/')) #"~/Desktop/trust_reports/")
list_files <- list.files(path = here::here("data/trust4/out_nhnn/"))
list_files

merged_file <- list()
for (i in list_files) {
  merged_file[[i]] <- read.table(paste0(here::here("data/trust4/out_nhnn/"), i), sep = "\t", header = F)
}

trust <- rbindlist(merged_file, id = T) #, fill = T)
names(trust) <- c(".id", "read_count", "frequency", "CDR3_dna", "CDR3_amino_acids", "V", "D", "J", "C", "consensus_id", "consensus_id_complete_vdj")

trust_2 <- trust %>% 
  mutate(matcher = gsub(".*(N\\d{1,4})(\\D).*", "\\1", .id)) %>%
  left_join(metadata_dir[,c(1,8)])
write.csv(trust_2, "~/Desktop/trust_nhnn.csv")


trust_annot <- trust %>% 
  mutate(matcher = gsub("^.*_(LP[0-9A-Za-z]+)_.*$", "\\1", .id)) %>%  #gsub("^TRUST_(LP[0-9]+)_.*$", "\\1", .id))# %>% 
  left_join(trust_meta) %>% 
  mutate(disease = ifelse(grepl("2608", matcher) | 
                            grepl("1123", matcher) |
                          grepl("1112", matcher) |
                          grepl("194", matcher), "Ctrl", "IBM"))

trust_tcr <- trust %>% filter(grepl("TR", V)) %>% filter("CDR3_amino_acids" != "out_of_frame")
trust_2 %>% 
  #filter(read_count > 10) %>% 
  mutate(tcr = ifelse(grepl("TR", V), "TCR", "IG")) %>% 
  filter(tcr == "TCR") %>% 
  #filter() %>% 
  ggplot(aes(x = reorder(V,read_count), y = read_count, fill = .id)) +
  geom_col(show.legend = F) +
  #facet_wrap(~tcr, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) #axis.text.x = element_text(angle = 90))


trust_2 %>% 
  #filter(grepl("TR", V)) %>% 
  filter("CDR3_amino_acids" != "out_of_frame") %>% 
  #filter(V1 > 10) %>% 
  mutate(tcr = ifelse(grepl("TR", V), "TCR", "IG")) %>% 
  filter(tcr == "TCR") %>% 
  #filter() %>% 
  ggplot(aes(x = sample_name, y = read_count, fill = V)) +
  geom_col(show.legend = F) +
  #facet_wrap(~tcr, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))


trust_2 %>% 
  #filter(read_count > 10) %>% 
  mutate(tcr = ifelse(grepl("TR", V), "TCR", "IG")) %>% 
  filter(tcr == "TCR") %>% 
  filter(CDR3_amino_acids != "out_of_frame") %>% 
  #filter() %>% 
  ggplot(aes(x = reorder(sample_name,read_count), y = read_count)) +
  geom_col(show.legend = F) +
  facet_wrap(~condition, scales = "free_x") +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) #axis.text.x = element_text(angle = 90))

trust_2 %>% 
  #filter(V1 > 10) %>% 
  mutate(tcr = ifelse(grepl("TR", V), "TCR", "IG")) %>% 
  #mutate(.id1 = grepl("\\_*", "", .id))
  filter(tcr == "TCR") %>% 
  filter(CDR3_amino_acids != "out_of_frame") %>% 
  #filter() %>% 
  ggplot(aes(x = reorder(.id, read_count), y = read_count, fill = CDR3_amino_acids)) +
  geom_col(show.legend = F) +
  facet_wrap(~condition, scales = "free_x") +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(axis.text.x = element_blank())







validated <- readxl::read_xlsx("~/Desktop/muscle_ibm_peptides/Share_with_Matteo/TRUST4_GLIPH_CLUSTERS_OF_INTEREST2.xlsx")
validated %>% 
  #filter(Freq >10) %>% 
  #mutate(condv = paste0(Condition, "_", V)) %>% 
  ggplot(aes(x = TCRb_Flag, y = Freq, fill = Condition)) +
  geom_col() +
  scale_fill_viridis_d() +
  facet_wrap(~Condition, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

validated %>% 
  filter(`HLA-B...26` == "B*08:01" | `HLA-B...27` == "B*08:01") %>% 
  mutate(tcr = ifelse(V == "TRBV7-9", "ok", "no")) %>% 
  ggplot(aes(x = Sample_ID, y = Freq, fill = J)) +
  geom_col() +
  scale_fill_viridis_d() +
  facet_wrap(tcr~Condition, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

trust_shortlist <- trust_annot %>% 
  right_join(validated)

trust_shortlist %>% 
  ggplot(aes(x = reorder(V,read_count), y = read_count)) +
  geom_col() +
  facet_wrap(~disease, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))



##### add analysis on tcrb-V types?