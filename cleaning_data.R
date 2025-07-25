library(data.table)
library(tidyverse)

setwd("~/Desktop/GCRF/GCRF_Pub/")
morpho <- read.csv("data/SNGCRF_master.csv") %>% filter(recap=="N") %>% distinct(band_number, .keep_all = TRUE) %>% select(-BGP_ID, -X)
band_ID <- read.csv("data/BGP_Band.csv", header=F, col.names = c("BGP_ID", "band_number"))
morpho_ID <- left_join(morpho, band_ID, by="band_number") %>% select(BGP_ID, everything(), -extracted, -seq_platform)
 
sample_list <- read.delim("data/sample_list.txt", header=F, col.names = "BGP_ID")

morpho_ID <- morpho_ID %>%
  mutate(sequenced = ifelse(BGP_ID %in% sample_list$BGP_ID, "Y", "N"))

morpho_ID[morpho_ID == "" | morpho_ID == "-"] <- NA

write.csv(morpho_ID, "data/raw_GCRF_data_master.csv", row.names = F, quote = F)

#==========================================#

morpho.full <- read.csv("data/raw_GCRF_data_master.csv")

morpho <- morpho.full %>% select(BGP_ID, band_number, site_code, sex, age, wing_chord, tarsus, tail_length, nare_length, bristle_length, beak_depth, beak_width, culmen_end_length) 

feather <- read.csv("data/feather_measurements.csv")

feather_avg <- feather %>% group_by(BGP_ID) %>% summarize(pen_length=mean(penn_length), pen_num=mean(penn_num), plum_length=mean(plum_length), plum_number=mean(plum_num), nodes=mean(plum_nodes))

morpho.bf <- left_join(morpho, feather_avg, by="BGP_ID")
nrow(morpho.bf)
#removing BGP_ID 23N00557 because it doesn't have band data or any morpho and just causes problems...
morpho.bf <- morpho.bf[-112,]
nrow(morpho.bf)
write.csv(morpho.bf, "data/raw_GCRF_data_all_morpho.csv")
