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
library(iNEXT)

load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")

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

dat<-as.data.frame(t(PERDIVER_taxo_all[,-c(1:3)]))
names(dat)<-paste(PERDIVER_taxo_all$Host, PERDIVER_taxo_all$UM, sep = '_')

z<-iNEXT(x = dat, 
         datatype = 'abundance',
         q = c(0,1,2))
p<-as.data.frame(do.call('cbind', z$iNextEst$size_based))
p$Host<-unlist(lapply(strsplit(p$Assemblage, '_'), function(x) x[1]))
p$UM<-unlist(lapply(strsplit(p$Assemblage, '_'), function(x) x[2]))
p$UM<-unlist(lapply(strsplit(p$UM, '\\.'), function(x) x[1]))

p %>% 
  filter(Order.q == 0 & Host != 'kracer') %>%
  mutate(across(c(m,qD:SC.LCL), as.numeric)) %>%
  ggplot(aes(x = m, y = qD, color = Host, fill = Host))+
  geom_line()+
  geom_ribbon(aes(x = m, ymax = qD.UCL, ymin = qD.LCL), alpha = 0.5)+
  facet_wrap(~Host+UM, ncol = 4)+
  theme_classic()

p %>% 
  filter(Order.q == 0 & Method == 'Observed'& Host != 'kracer') %>%
  pull(SC) %>%
  as.numeric(.) %>%
  mean()

zotu <- read.delim("C:/Users/18172844S/Documents/zotutab_raw_TaxFiltered_wtax_PERDIVER.txt")
rownames(zotu)<-zotu$ZOTUId
zotu<-zotu[!colnames(zotu) %in% c('ZOTUId', 'Taxonomy')]
#out2 <- iNEXT(zotu, q=0, datatype="abundance")
#save(out2, file = 'Completeness_zOTU.RData')
sum(zotu)

load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/Completeness_zOTU.RData")
Equivalencias_muestras_suelo <- read_excel("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Datos/Equivalencias_muestras_suelo.xlsx")

Equivalencias_muestras_suelo <- Equivalencias_muestras_suelo %>%
  mutate(Host = tolower(unlist(lapply(strsplit(Population, '_'), function(x) x[1]))),
         UM = ifelse(startsWith(UM, 'S'), 'S', UM),
         Sample = str_remove(Sample, ' ')) %>%
  filter(Sample %in% unique(zotu$Assemblage))

zotu<-as.data.frame(do.call('cbind', out2$iNextEst$size_based))
zotu %>% 
  filter(Order.q == 0) %>%
  mutate(across(c(m,qD:SC.UCL), as.numeric)) %>%
  rename(Sample = Assemblage) %>%
  left_join((Equivalencias_muestras_suelo %>% dplyr::select(-Population)), 'Sample') %>%
  filter(Method == 'Observed') %>%
  group_by(UM, Host) %>%
  summarize(qD = mean(qD),
            qD.LCL = mean(qD.LCL),
            qD.UCL = mean(qD.UCL),
            SC = mean(SC),
            SC.LCL = mean(SC.LCL),
            SC.UCL = mean(SC.UCL)) %>%
  distinct() %>%
  ggplot(aes(x = m, y = qD, color = Host, fill = Host))+
  geom_line()+
  geom_ribbon(aes(x = m, ymax = qD.UCL, ymin = qD.LCL), alpha = 0.5)+
  facet_wrap(~Host+UM, ncol = 4)+
  theme_classic()

  
out2$DataInfo %>%
  rename(Sample = Assemblage) %>%
  left_join((Equivalencias_muestras_suelo %>% dplyr::select(-Population)), 'Sample') %>%
  distinct() %>%
  group_by(UM, Host) %>%
  summarize(s = mean(S.obs),
            s_sd = sd(S.obs),
            SC_ = mean(SC),
            sd = sd(SC)) %>% View()

sd(100*out2$DataInfo$SC)

ChaoRichness(zotu) %>%
  mutate(Sample = rownames(.)) %>%
  left_join((Equivalencias_muestras_suelo %>% dplyr::select(-Population)), 'Sample') %>%
  distinct() %>%
  group_by(UM, Host) %>%
  summarize(s = mean(Observed)) %>% View()
  