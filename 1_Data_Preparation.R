library(readxl)
library(tidyverse)
library(taxize)
library(rgbif)
Sys.setenv(ENTREZ_KEY='7c76967d93d4938a1c715afc62fe297a8308')

PERDIVER<-read_excel("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Datos/Arthropods/PERDIVER_unificado.xlsx", sheet = 'Tabla_PERDIVER')
PERDIVER<-read_excel("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Datos/Arthropods/PERDIVER_unificado.xlsx", sheet = 'Tabla_PERDIVER')


PERDIVER<-PERDIVER %>%
  mutate(Genus = unlist(lapply(strsplit(Species, '_'), function(x) x[1]))) %>%
  filter(!Species %in% c('Sin_observaciones', '-')) %>%
  filter(!Genus %in% c('Desconocido', 'Unknown'))

PERDIVER_Stand_NCBI<-as.data.frame(sapply(unique(PERDIVER$Genus), 
          function (x) tax_name(sci = x, 
                                get = c('order', 'family', 'subfamily'), 
                                db = 'ncbi')))

PERDIVER<-PERDIVER_Stand_NCBI %>%
  unlist() %>%
  matrix(ncol = 5, byrow = T) %>%
  as.data.frame() %>%
  dplyr::select(-V1) %>%
  rename(Genus = V2, Order_NCBI = V3, 
         Family_NCBI = V4, Subfamily_NCBI = V5) %>%
  left_join(PERDIVER, 'Genus') %>%
  mutate(Order_final = ifelse(is.na(Order_NCBI), Order, Order_NCBI)) %>%
  mutate(Family_final = ifelse(is.na(Family_NCBI), Family, Family_NCBI))

PERDIVER[PERDIVER$Species == 'Gisinurus_sp.1',]$Order_final<-'Symphypleona'

save(PERDIVER, file = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData')
