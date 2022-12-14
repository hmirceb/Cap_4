library(vegan)
library(gridExtra)
library(picante)
library(usedist)
library(tidyverse)
library(popbio)
library(PDcalc)
library(readxl)
library(ggpubr)
library(hillR)
library(spaa)
rm(list = ls())
setwd('~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions')
setwd('C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions')


load("Resultados/PERDIVER_Comm_Data.RData")
load("Resultados/Filo_PERDIVER_ALT.Rdata")
plantas_PERDIVER<-read_excel("Datos/Plantas/Habitat_PERDIVER.xlsx")
PERDIVER_abiotic<-read_excel("Datos/PERDIVER_abiotic.xlsx")
load("Resultados/Joan/zOTUs_Alfadiversidad.Rdata")
load("Resultados/Joan/ses.mpd.output.all_bySp.RData")
load("/home/hector/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/Joan/PhylogeneticDist_phyloresembl_14000Reads.Rdata")
load("/home/hector/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/Joan/TaxonomicDist_14000Reads.Rdata")

# Distancia filo entre especies media (de todos los arboles)
cophenetic_mean<-mean.list(lapply(filo_perdiver, cophenetic))

# Matriz de sitios por especies

PERDIVER_sitesXsps<-PERDIVER  %>%
  dplyr::select(c(Host, UM, Species)) %>%
  mutate(UM_s = ifelse(startsWith(UM, 'S'), 'S', 'L')) %>%
  group_by(Host, UM, Species, UM_s) %>%
  dplyr::summarize(Host,
            UM,
            UM_s,
            Species,
            n = n()) %>%
  distinct() %>%
  pivot_wider(names_from = Species,
              id_cols = c(Host, UM, UM_s),
              values_from = n,
              values_fill = 0) %>%
  ungroup()
PERDIVER_sitesXsps<-as.data.frame(PERDIVER_sitesXsps)
rownames(PERDIVER_sitesXsps)<-paste(PERDIVER_sitesXsps$Host, PERDIVER_sitesXsps$UM, sep = '_')

PERDIVER_sitesXsps_filo<-PERDIVER %>%
  filter(Clase == 'Insecta') %>%
  dplyr::select(c(Host, UM, Species)) %>%
  mutate(UM_s = ifelse(startsWith(UM, 'S'), 'S', 'L')) %>%
  group_by(Host, UM, Species, UM_s) %>%
  dplyr::summarize(Host = Host,
                   UM = UM,
                   UM_s = UM_s,
                   Species,
                   n = n()) %>%
  distinct() %>%
  pivot_wider(names_from = Species,
              id_cols = c(Host, UM, UM_s),
              values_from = n,
              values_fill = 0) %>%
  ungroup()
PERDIVER_sitesXsps_filo<-as.data.frame(PERDIVER_sitesXsps_filo)
rownames(PERDIVER_sitesXsps_filo)<-paste(PERDIVER_sitesXsps_filo$Host, PERDIVER_sitesXsps_filo$UM, sep = '_')


# Matriz de sitios por especies para las plantas coocurrentes
PLANTS_sitesXsps<-plantas_PERDIVER %>%
  dplyr::select(c(Host, UM, Especie, Abund_perc)) %>%
  mutate(UM_s = ifelse(startsWith(UM, 'S'), 'S', 'L')) %>%
  mutate(Abund_perc = as.numeric(Abund_perc)) %>%
  pivot_wider(names_from = Especie,
              id_cols = c(Host, UM, UM_s),
              values_from = Abund_perc,
              values_fill = 0) %>%
  ungroup()

PLANTS_alpha_taxo<-PLANTS_sitesXsps %>%
  dplyr::select(c(Host, UM, UM_s)) %>%
  mutate(plant_riq = specnumber(PLANTS_sitesXsps[,-c(1:3)]))

# Alpha diversity with Hill numbers

# Taxo interactuantes
qs<-seq(0, 2, by = 0.1)
a<-c()
for (i in 1:length(qs)) {
  a[[i]]<-hill_taxa(PERDIVER_sitesXsps_filo[,-c(1:3)], q = qs[[i]])
}
PERDIVER_alpha_taxo<-as.data.frame(cbind(PERDIVER_sitesXsps_filo[,c(1:3)], 
                                         do.call('cbind', a)))
names(PERDIVER_alpha_taxo)[-c(1:3)]<-as.character(qs)

# Taxo plantas
qs<-seq(0, 2, by = 0.1)
a<-c()
for (i in 1:length(qs)) {
  a[[i]]<-hill_taxa(PLANTS_sitesXsps[,-c(1:3)], q = qs[[i]])
}
PLANTS_alpha_taxo<-as.data.frame(cbind(PLANTS_sitesXsps[,c(1:3)], 
                                       do.call('cbind', a)))
names(PLANTS_alpha_taxo)[-c(1:3)]<-paste(as.character(qs), 'plants', sep = '_')

# Juntamos en un data frame y sacamos graficos
PERDIVER_alpha_taxo<-PERDIVER_alpha_taxo %>%
  left_join(PLANTS_alpha_taxo)


### TAXO Beta diversity
# Distancia con bray curtis (considera las abundancias)
PERDIVER_beta_taxo<-vegdist(decostand(PERDIVER_sitesXsps[,-c(1:3)], 'hellinger'),
                            method = 'bray')

###############################
#### phylogenetic analysis ####
##############################

# Calculo de MPD y su SES. Se hace por separado para usar un pool de espcies diferente 
# por cada especie porque si no saldria raro

temp_mpd<-list()
for (i in 1:length(unique(PERDIVER$Host))) {
  temp<-PERDIVER  %>%
    filter(Host == unique(PERDIVER$Host)[[i]]) %>%
    dplyr::select(c(Host, UM, Species)) %>%
    filter(Species %in% filo_perdiver[[1]]$tip.label) %>%
    mutate(UM_s = ifelse(startsWith(UM, 'S'), 'S', 'L')) %>%
    group_by(Host, UM, Species, UM_s) %>%
    dplyr::summarize(Host = Host,
                     UM = UM,
                     UM_s = UM_s,
                     Species,
                     n = n()) %>%
    distinct() %>%
    pivot_wider(names_from = Species,
                id_cols = c(Host, UM, UM_s),
                values_from = n,
                values_fill = 0) %>%
    ungroup()
  
  temp_mpd[[i]]<-cbind(as.data.frame(temp[,c(1:3)]), 
                       ses.mpd(samp = as.data.frame(temp[,-c(1:3)]),
                               dis = cophenetic_mean, 
                               abundance.weighted = T))
}

PERDIVER_alpha_filo<-as.data.frame(do.call('rbind', temp_mpd))

### Beta diversity with PDcalc
PERDIVER_beta_filo_list<-lapply(filo_perdiver, function(y) phyloresembl(x = decostand(PERDIVER_sitesXsps_filo[,-c(1:3)], 'hellinger'),
                                                                        phy = y, incidence = F, 'sorensen'))

PERDIVER_beta_filo<-as.dist(mean.list(lapply(PERDIVER_beta_filo_list, function(x) as.matrix(x))))

#### Tabla resumen ####

PERDIVER_taxo_all<-PERDIVER %>%
  dplyr::select(c(Host, UM, Species)) %>%
  mutate(UM_s = ifelse(startsWith(UM, 'S'), 'S', 'L')) %>%
  group_by(Host, UM, Species, UM_s) %>%
  dplyr::summarize(Host = Host,
                   UM = UM,
                   UM_s = UM_s,
                   Species,
                   n = n()) %>%
  distinct() %>%
  pivot_wider(names_from = Species,
              id_cols = c(Host, UM, UM_s),
              values_from = n,
              values_fill = 0) %>%
  ungroup()

qs<-seq(0, 2, by = 0.1)
a<-c()
for (i in 1:length(qs)) {
  a[[i]]<-hill_taxa(PERDIVER_taxo_all[,-c(1:3)], q = qs[[i]])
}
PERDIVER_alpha_taxo_all<-as.data.frame(cbind(PERDIVER_taxo_all[,c(1:3)], 
                                             do.call('cbind', a)))
names(PERDIVER_alpha_taxo_all)[-c(1:3)]<-as.character(qs)

PERDIVER_summary<-PERDIVER_alpha_taxo_all %>%
  dplyr::select(c(Host, UM, UM_s, `0`, `1`, `2`)) %>%
  mutate(Above_richness = `0`) %>%
  mutate(Above_shannon = log(`1`)) %>%
  mutate(Above_simpson = 1-(1/`2`)) %>%
  dplyr::select(-c(`0`, `1`, `2`)) %>%
  full_join(PERDIVER_alpha_taxo, c('Host', 'UM', 'UM_s')) %>%
  left_join(PERDIVER_abiotic, c('Host', 'UM')) %>%
  mutate(Insect_richness = `0`) %>%
  mutate(Insect_shannon = log(`1`)) %>%
  mutate(Insect_simpson = 1-(1/`2`)) %>%
  mutate(Plant_richness = `0_plants`) %>%
  mutate(Plant_shannon = log(`1_plants`)) %>%
  mutate(Plant_simpson = 1-(1/`2_plants`)) %>%
  dplyr::select(c('Host':'UM_s',
                  'Above_richness',
                  'Above_shannon',
                  'Above_simpson',
                  'Insect_richness',
                  'Insect_shannon',
                  'Insect_simpson',
                  'Plant_richness',
                  'Plant_shannon',
                  'Plant_simpson',
                  '%_cobertura_vegetal_(suelo)':'Perturbaciones')) %>%
  left_join(PERDIVER_alpha_filo, c('Host', 'UM', 'UM_s')) %>%
  dplyr::select(-c(ntaxa, runs)) %>% 
  dplyr::rename(Insect_mpd.obs = mpd.obs) %>%
  dplyr::rename(Insect_mpd.obs.z = mpd.obs.z) %>%
  dplyr::rename(Insect_mpd.obs.p = mpd.obs.p) %>%
  left_join(PERDIVER %>%
              group_by(Host, UM) %>%
              dplyr::summarize(Host = Host,
                        UM = UM,
                        Above_n_inter = n()) %>%
              distinct()) %>%
  left_join(PERDIVER %>%
              filter(Clase == 'Insecta') %>%
              group_by(Host, UM) %>%
              dplyr::summarize(Host = Host,
                        UM = UM,
                        Insect_n_inter = n()) %>%
              distinct())

###########################
#### Unimos data below ####
###########################
ses.mpd.output.all$sample<-rownames(ses.mpd.output.all)
metadata.alfa$sample<-rownames(metadata.alfa)

# Corregimos un pequeño error con las pobalciones de kracer
metadata.alfa[rownames(metadata.alfa) %in% c('PER08', 'PER05'),]$Poblacio_indiv<-'Large3'
metadata.alfa[rownames(metadata.alfa) %in% c('PER08', 'PER05'),]$Poblacio_indiv.1<-'Large3'
metadata.alfa[rownames(metadata.alfa) %in% c('PER08', 'PER05'),]$UM<-'L3'


dat<-metadata.alfa %>%
  full_join(ses.mpd.output.all, by = 'sample')  %>%
  dplyr::select(c(Host, UM, UM_s, 
                  richness.14000,
                  shannon.14000,
                  simpson.14000,
                  mpd.obs,
                  mpd.obs.z)) %>% 
  group_by(Host, UM) %>%
  dplyr::summarize(Host = Host,
            UM = UM,
            UM_s = UM_s,
            Below_richness = mean(richness.14000, na.rm = T),
            Below_shannon = mean(shannon.14000, na.rm = T),
            Below_simpson = mean(simpson.14000, na.rm = T),
            Below_mpd.obs = mean(mpd.obs, na.rm = T),
            Below_mpd.obs.z = mean(mpd.obs.z, na.rm = T)) %>%
  distinct()

PERDIVER_summary<-PERDIVER_summary %>%
  left_join(dat, c('Host', 'UM', 'UM_s'))

# Calculamos la distancia media entre sitios segun la comunidad de bacterias.
# Hay que hacerlo porque hay varias muestras de cada especie en cada poblacion

id_plots<-inner_join(metadata.alfa, ses.mpd.output.all, 'sample') %>%
  dplyr::select(c(Host, UM, UM_s, sample)) %>%
  mutate(plot = paste(Host, UM, sep = '_'))

a<-dist2list(otutab.dist)
a$col<-as.character(a$col)
a$row<-as.character(a$row)

z<-a %>%
  left_join((id_plots %>%
               dplyr::rename(col = sample) %>%
               dplyr::select(c(col, plot))),
            'col') %>%
  left_join((id_plots %>%
               rename(row = sample) %>%
               dplyr::select(c(row, plot))),
            'row') %>%
  mutate(group = paste(plot.x, plot.y, sep = '-')) %>%
  group_by(group) %>%
  dplyr::summarize(group = group,
            value = mean(value)) %>%
  distinct() %>%
  ungroup() %>%
  mutate(row = unlist(lapply(strsplit(group, '-'), function(x) x[1])),
         col = unlist(lapply(strsplit(group, '-'), function(x) x[2]))) %>%
  mutate(value = ifelse(row == col, 0, value)) %>%
  dplyr::select(c(col, row, value)) %>%
  mutate(col = as.factor(col),
         row = as.factor(row))  %>%
  as.data.frame()

PERDIVER_beta_taxo_below<-list2dist(z)

# Igual pero para la filo
b<-dist2list(phyl.dist)
b$col<-as.character(b$col)
b$row<-as.character(b$row)

q<-b %>%
  left_join((id_plots %>%
               rename(col = sample) %>%
               dplyr::select(c(col, plot))),
            'col') %>%
  left_join((id_plots %>%
               rename(row = sample) %>%
               dplyr::select(c(row, plot))),
            'row') %>%
  mutate(group = paste(plot.x, plot.y, sep = '-')) %>%
  group_by(group) %>%
  dplyr::summarize(group = group,
            value = mean(value)) %>%
  distinct() %>%
  ungroup() %>%
  mutate(row = unlist(lapply(strsplit(group, '-'), function(x) x[1])),
         col = unlist(lapply(strsplit(group, '-'), function(x) x[2]))) %>%
  mutate(value = ifelse(row == col, 0, value)) %>%
  dplyr::select(c(col, row, value)) %>%
  mutate(col = as.factor(col),
         row = as.factor(row))  %>%
  as.data.frame()

PERDIVER_beta_filo_below<-list2dist(q)

### Beta diversidad de plantas coocurrentes
PERDIVER_beta_plants<-vegdist(decostand(PLANTS_sitesXsps[,-c(1:3)], 'hellinger'),
                            method = 'bray')

save(PERDIVER_summary, 
     PERDIVER_beta_taxo,
     PERDIVER_beta_filo,
     PERDIVER_beta_taxo_below,
     PERDIVER_beta_filo_below,
     PERDIVER_beta_plants,
     file = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/Diversity_results.RData')

save(PERDIVER_summary, 
     PERDIVER_beta_taxo,
     PERDIVER_beta_filo, 
     PERDIVER_beta_taxo_below,
     PERDIVER_beta_filo_below,
     PERDIVER_beta_plants,
     file = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/Diversity_results.RData')

