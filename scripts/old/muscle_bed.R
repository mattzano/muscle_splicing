###muscle junction pickup
library(data.table)
n <- c(1:6)
sample_set <- c(paste0("Control_",n),paste0("KD_",n))
mus <- list()
for (i in sample_set) {
  mus[[i]] <- read.table(paste0("~/Documents/phd/research_lines/muscle/muscle_i3_bed/iMuscle", i, ".bed"), sep = "\t")
}
muscle <- rbindlist(mus, id = T) %>% 
  filter(V4 != "QTRT2") %>% 
  mutate(paste_into_igv = paste0(V1,":",V2,"-",V3)) %>% 
  pivot_wider(names_from = .id, values_from = V6) %>% 
  mutate(controls = rowMeans(select(.,7:12)),
         kd = rowMeans(select(.,13:18)),
         delta = kd-controls,
         sum = kd+controls,
         hop = delta/sum,
         symbol = V4) %>% 
  mutate(length = V3-V2,
         aa_length = length/3) %>%
  mutate(aa_good = ifelse(length%%3 == 0, "yes", "no")) %>% 
  filter(delta > 0)
  #pivot_longer(cols = c(controls,kd))


muscle %>%
  pivot_longer(cols = c(controls,kd)) %>% 
  ggplot(aes(x = name, y = value)) + 
  geom_point() +
  geom_boxplot() +
  #scale_y_log10() +
  theme_bw()


majiq <- read.csv("~/Desktop/Control-TDP43KD_annotated_junctions_filter.csv")

muscle2 <- muscle %>% left_join(majiq)

##you should add myoblasts and splicing data from patients

##match ce reads with proteomics patients (exclude genes not expressed)
chr15:77640362-77641025







