library(tidyverse)
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
sps_unique<-PERDIVER %>%
  select(Clase, Order_final, Family_final, Genus, Species) %>%
  distinct()
  

morfo_a_sp<-sps_unique %>%
  filter(!grepl('sp.', Species))

morfo_a_genero<-sps_unique %>%
  filter(grepl('sp.', Species) &
           Genus != Family_final &
           Genus != Order_final)

morfo_a_familia<-sps_unique %>%
  filter(grepl('sp.', Species) &
           Genus == Family_final)

morfo_a_orden<-sps_unique %>%
  filter(grepl('sp.', Species) &
           Genus == Order_final)

a<-rbind(morfo_a_sp, morfo_a_genero, morfo_a_familia, morfo_a_orden)

morfo_a_otros<-sps_unique %>%
  filter(!Species %in% a$Species)

inters<-data.frame(Level = c('Species', 'Genus', 'Family', 'Order', 'Otros'),
           N = c(dim(morfo_a_sp)[1], dim(morfo_a_genero)[1],
                 dim(morfo_a_familia)[1], dim(morfo_a_orden)[1],
                 dim(morfo_a_otros)[1]))
inters$Prop<-100*inters$N/sum(inters$N)




krasche<-PERDIVER %>%
  filter(Host == 'kracer')

