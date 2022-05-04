library(lme4)
library(visreg)
library(PhyloMeasures)
library(picante)
library(tidyverse)
library(popbio)
library(corrplot)

load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER.Rdata")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER_Family.Rdata")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER_Order.Rdata")

load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER.Rdata")
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER_Family.Rdata")
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER_Order.Rdata")

PERDIVER_sitesXsps_filo<-PERDIVER  %>%
  dplyr::select(c(Host, UM, Species)) %>%
  filter(Species %in% filo_perdiver[[1]]$tip.label) %>%
  mutate(UM_s = ifelse(startsWith(UM, 'S'), 'S', 'L')) %>%
  group_by(Host, UM, Species, UM_s) %>%
  summarize(Host = Host,
            UM = UM,
            UM_s = UM_s,
            Species,
            n = n(),
            pres = 1) %>%
  distinct() %>%
  pivot_wider(names_from = Species,
              id_cols = c(Host, UM, UM_s),
              values_from = pres,
              values_fill = 0) %>%
  ungroup()


obs_host<-list()
for (q in 1:length(unique(PERDIVER_sitesXsps_filo$Host))) {
  temp1<-PERDIVER_sitesXsps_filo %>%
    filter(Host == unique(PERDIVER_sitesXsps_filo$Host)[q]) %>%
    dplyr::select(-c(Host, UM, UM_s))
  
  obs_host[[q]]<-cbind(popbio::mean.list(lapply(filo_perdiver, function(x) matrix(
    pd.query(
      tree = x, matrix = temp1),
    ncol = 1))),
    popbio::mean.list(lapply(filo_perdiver, function(x) matrix(
      pd.query(
        tree = x, matrix = temp1,
        standardize = T),
      ncol = 1))))
}
obs_host<-do.call('rbind', obs_host)

obs_host_fam<-list()
for (q in 1:length(unique(PERDIVER_sitesXsps_filo$Host))) {
  temp1<-PERDIVER_sitesXsps_filo %>%
    filter(Host == unique(PERDIVER_sitesXsps_filo$Host)[q]) %>%
    dplyr::select(-c(Host, UM, UM_s))
  
  obs_host_fam[[q]]<-cbind(popbio::mean.list(lapply(filo_perdiver_fam, function(x) matrix(
    pd.query(
      tree = x, matrix = temp1),
    ncol = 1))),
    popbio::mean.list(lapply(filo_perdiver_fam, function(x) matrix(
      pd.query(
        tree = x, matrix = temp1,
        standardize = T),
      ncol = 1))))
}
obs_host_fam<-do.call('rbind', obs_host_fam)


obs_host_ord<-list()
for (q in 1:length(unique(PERDIVER_sitesXsps_filo$Host))) {
  temp1<-PERDIVER_sitesXsps_filo %>%
    filter(Host == unique(PERDIVER_sitesXsps_filo$Host)[q]) %>%
    dplyr::select(-c(Host, UM, UM_s))
  
  obs_host_ord[[q]]<-cbind(popbio::mean.list(lapply(filo_perdiver_ord, function(x) matrix(
    pd.query(
      tree = x, matrix = temp1),
    ncol = 1))),
    popbio::mean.list(lapply(filo_perdiver_ord, function(x) matrix(
      pd.query(
        tree = x, matrix = temp1,
        standardize = T),
      ncol = 1))))
}
obs_host_ord<-do.call('rbind', obs_host_ord)

b<-cbind(PERDIVER_sitesXsps_filo[,c(1:3)], obs_host, obs_host_fam, obs_host_ord)
names(b)[4:9]<-c('obs_host', 'st_obs_host', 
                 'obs_host_fam', 'st_obs_host_fam',
                 'obs_host_ord', 'st_obs_host_ord')

rm(filo_perdiver, filo_perdiver_fam, filo_perdiver_ord)

corrplot(cor(b[,-c(1:3)]), method = 'color',
         addCoef.col = 'black',
         tl.col = 'black',
         tl.srt = 45, diag = F,
         na.label = ' ',
         number.digits = 3,
         type = 'lower')

b %>%
  group_by(Host, UM_s) %>%
  pivot_longer(cols = obs_host:st_obs_host_ord) %>%
  filter(!startsWith(name, 'st_')) %>%
  ggplot(aes(x = Host, y = value, color = name, shape = UM_s))+
  geom_point()

m<-b %>%
  group_by(Host, UM_s) %>%
  pivot_longer(cols = obs_host:st_obs_host_ord) %>%
  filter(startsWith(name, 'st_')) %>%
  lm(value~name, data = .)
summary(m)
visreg(m)
car::Anova(m)
