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
p<-do.call('rbind', z$iNextEst)
p$Host<-unlist(lapply(strsplit(rownames(p), '_'), function(x) x[1]))
p$UM<-unlist(lapply(strsplit(rownames(p), '_'), function(x) x[2]))
p$UM<-unlist(lapply(strsplit(p$UM, '\\.'), function(x) x[1]))

ggiNEXT(z, type = 2,)+facet_wrap(~site, ncol = 4)

ggplot(p %>% filter(order == 0 & Host != 'kracer'), aes(x = m, y = qD, color = Host, fill = Host))+
  geom_line()+
  geom_ribbon(aes(x = m, ymax = qD.UCL, ymin = qD.LCL), alpha = 0.5)+
  facet_wrap(~Host+UM, ncol = 4)+
  theme_classic()

ggplot(p %>% filter(order == 0 & Host != 'kracer' &
                      method != 'extrapolated'), aes(x = m, y = SC, color = Host, fill = Host))+
  geom_line(aes(linetype = method))+
  geom_point(data = p %>% 
               filter(order == 0 & Host != 'kracer' & method == 'observed'),
             aes(x = m, y = SC),
             inherit.aes = F)+
  geom_ribbon(aes(x = m, ymax = SC.UCL, ymin = SC.LCL), alpha = 0.3)+
  facet_wrap(~Host+UM, ncol = 4)+
  theme_classic()

p %>%
  filter(method == 'observed' & 
           Host != 'kracer') %>%
  pull(SC) %>%
  sd()
