library(readxl)
library(tidyverse)
library(taxize)
library(taxizedb)
library(rgbif)
Sys.setenv(ENTREZ_KEY='7c76967d93d4938a1c715afc62fe297a8308')
rm(list = ls())

PERDIVER<-read_excel("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Datos/Arthropods/PERDIVER_unificado.xlsx", sheet = 'Tabla_PERDIVER')
PERDIVER<-read_excel("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Datos/Arthropods/PERDIVER_unificado.xlsx", sheet = 'Tabla_PERDIVER')


PERDIVER<-PERDIVER %>%
  mutate(Genus = unlist(lapply(strsplit(Species, '_'), function(x) x[1]))) %>%
  filter(!Species %in% c('Sin_observaciones', '-')) %>%
  filter(!Genus %in% c('Desconocido', 'Unknown'))

#### con gbif ####
aaa<-name_backbone_checklist(data.frame(name = unique(PERDIVER$Species), phylum = 'Arthropoda'))

bbb<-aaa %>%
  dplyr::select(c(verbatim_name, species, genus, family, order, phylum)) %>%
  mutate(species = str_replace(species, ' ', '_')) %>%
  rename(Species = verbatim_name, gbif_species = species, gbif_genus = genus, gbif_family = family, gbif_order = order, gbif_phylum = phylum) %>%
  full_join(PERDIVER, 'Species')

PERDIVER<-bbb %>%
  mutate(species_final = ifelse(is.na(gbif_species) & !is.na(gbif_genus),
                               paste(gbif_genus, 
                                     unlist(lapply(strsplit(Species, '_'), function(x) x[2])), sep = '_'),
                               ifelse(is.na(gbif_species) & is.na(gbif_genus) & !is.na(gbif_family),
                                      paste(gbif_family, 
                                            unlist(lapply(strsplit(Species, '_'), function(x) x[2])), sep = '_'),
                                      ifelse(is.na(gbif_species) & is.na(gbif_genus) & is.na(gbif_family) & !is.na(gbif_order),
                                             paste(gbif_order, 
                                                   unlist(lapply(strsplit(Species, '_'), function(x) x[2])), sep = '_'),
                                             gbif_species)))) %>%
  mutate(species_final = ifelse(is.na(species_final), Species, species_final)) %>%
  mutate(gbif_family = ifelse(is.na(species_final) & 
                                is.na(gbif_family) & 
                                !is.na(Family),
                              Family, 
                              gbif_family),
         gbif_order = ifelse(is.na(gbif_order) & 
                               !is.na(Order),
                             Order, 
                             gbif_order)) %>%
  mutate(gbif_family = ifelse(!is.na(species_final) & 
                                is.na(gbif_family) & 
                                !is.na(Family),
                              Family, 
                              gbif_family),
         gbif_order = ifelse(is.na(gbif_order) & 
                               !is.na(Order),
                             Order, 
                             gbif_order)) %>%
  dplyr::select(c(Cod_censo, Host, UM, 
                  Species_original, Species, species_final,
                  gbif_family, gbif_order,
                  Clase))
PERDIVER$genus<-unlist(lapply(strsplit(PERDIVER$species_final, '_'), function(x) x[1]))

PERDIVER$gbif_order<-ifelse(PERDIVER$gbif_order == 'Acari', NA, PERDIVER$gbif_order)
PERDIVER$gbif_order<-ifelse(startsWith(PERDIVER$species_final, 'Oribatida'), 'Sarcoptiformes', PERDIVER$gbif_order)
PERDIVER$gbif_order<-ifelse(PERDIVER$gbif_family %in% c('Tomoceridae', 'Orchesellidae', 'Entomobryidae'), 
                                  'Entomobryomorpha', PERDIVER$gbif_order)
PERDIVER$gbif_order<-ifelse(is.na(PERDIVER$gbif_order) & PERDIVER$gbif_family == 'Neoliolidae',
                                  'Sarcoptiformes', PERDIVER$gbif_order)  
PERDIVER$gbif_order<-ifelse(PERDIVER$genus == 'Coccoidea', 'Hemiptera', PERDIVER$gbif_order)
PERDIVER$gbif_family<-ifelse(PERDIVER$genus == 'Coccoidea', NA, PERDIVER$gbif_family)

PERDIVER$Species<-PERDIVER$species_final
PERDIVER<-PERDIVER %>% 
  dplyr::select(-species_final) %>%
  rename(family = gbif_family, order = gbif_order)


save(PERDIVER, file = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData')
save(PERDIVER, file = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData')
