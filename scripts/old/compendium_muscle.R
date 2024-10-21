myoblasts <- read.csv("~/Desktop/controlmyoblast-TDP43KDmyoblast_annotated_junctions.csv") %>% 
  mutate(exp = "myoblasts")
names(myoblasts)[c(12,13)] <- c("control_mean_psi", "tdp43kd_mean_psi")

i3 <- read.csv("~/Desktop/Control-TDP43KD_annotated_junctions.csv") %>% 
  mutate(exp = "i3muscle")

muscle_invitro <- rbind(myoblasts,i3)



myoblast_sig <- myoblasts %>% 
  dplyr::filter(mean_dpsi_per_lsv_junction > 0.1 & control_mean_psi < 0.05) %>% 
  dplyr::filter(!gene_name %in% c("AHNAK", "CH507-528H12.1", "CH507-513H4.1"))

i3_sig <- i3 %>% 
  dplyr::filter(mean_dpsi_per_lsv_junction > 0.1 & control_mean_psi < 0.05)
