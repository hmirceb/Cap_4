library(readxl)
library(tidyverse)

PERDIVER_bichos<-read_excel("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Datos/PERDIVER_bichos.xlsx")
PERDIVER_bichos<-read_excel("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Datos/PERDIVER_bichos.xlsx")

PERDIVER_jens<-read_excel("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Datos/PERDIVER_Jens.xlsx")
PERDIVER_jens<-read_excel("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Datos/PERDIVER_Jens.xlsx")

Maria<-PERDIVER_bichos %>%
  mutate(Genus = unlist(lapply(strsplit(Morfoespecie, '_'), function(x) x[1]))) %>%
  filter(!Morfoespecie %in% c('Sin_observaciones', '-')) %>%
  filter(!Genus %in% c('Desconocido', 'Unknown')) %>%
  dplyr::select(c(Host, Morfoespecie)) %>%
  rename(Species = Morfoespecie) %>%
  mutate(Person = 'Maria') %>%
  distinct()

Jens<-PERDIVER_jens %>%
  filter(Plant != 'Krashninikovia') %>%
  dplyr::select(c(Plant, Species)) %>%
  rename(Host = Plant) %>%
  mutate(Person = 'Jens') %>%
  distinct()

# juntamos las dos listas
Maria_Jens<-Maria %>%
  full_join(Jens) %>%
  mutate(Genus = unlist(lapply(strsplit(Species, '-'), function(x) x[1]))) %>% 
  mutate(Plant = unlist(lapply(strsplit(Species, '-'), function(x) x[2])))

# Vemos cuales tienen
Maria_Jens %>%
  group_by(Species) %>%
  filter(n() == 1) %>%
  filter(!is.na(Plant)) %>%
  View()

Maria_Jens %>%
  group_by(Species) %>%
  filter(n() == 1) %>%
  filter(!is.na(Plant)) %>%
  View()

a<-Maria_Jens %>%
  group_by(Species) %>%
  filter(n() == 1) %>%
  filter(!is.na(Plant)) 
a[a$Species %in% a$Genus,]

Maria_Jens %>%
  group_by(Species) %>%
  filter(n() == 1) %>%
  filter(is.na(Plant)) %>%
  View()

