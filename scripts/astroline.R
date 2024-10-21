###astroline

###de
salmon_deseq2(sample_dir =  here::here('data', 'salmon'),
              metadata_dir =  here::here('data','metadata_muscle_astroline.csv'),
              experiment_dir = "astroline")
#write.csv(results_table, "~/Desktop/de_astroline.csv")

results_table %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point()

normed_counts %>% 
  filter(symbol == "TARDBP") %>% 
  pivot_longer(cols = c(2:5)) %>% 
  ggplot(aes(x = name, y = value)) +
  geom_point()

##ds
splicing_astroline <- read.csv(here::here('data', 'splicing', 'Control-TDP43KD_annotated_junctions.csv'))

splicing_astroline %>% 
  ggplot(aes(x = mean_dpsi_per_lsv_junction, y = -log10(1-probability_changing))) +
  geom_point()
