library(tidyverse)
library(BAT)
library(betapart)
library(vegan)
library(lme4)
library(rstanarm)
library(visreg)
library(gridExtra)
library(PhyloMeasures)
library(picante)

load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER.Rdata")

load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER.Rdata")

PERDIVER_sitesXsps<-PERDIVER  %>%
  dplyr::select(c(Host, UM, Species)) %>%
  mutate(UM_s = ifelse(startsWith(UM, 'S'), 'S', 'L')) %>%
  group_by(Host, UM, Species) %>%
  summarize(Host = Host,
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
  

# Alpha diveristy with Hill numbers
PERDIVER_alpha_taxo<-data.frame(PERDIVER_sitesXsps[,c(1:3)], 
           hill(comm = PERDIVER_sitesXsps[,-c(1:3)], q = 0), 
           hill(comm = PERDIVER_sitesXsps[,-c(1:3)], q = 1), 
           hill(comm = PERDIVER_sitesXsps[,-c(1:3)], q = 2))

PERDIVER_alpha_taxo %>%
  pivot_longer(cols = c(Hill.0:Hill.2), 
               names_to = 'q') %>%
  mutate(q = ifelse(q == 'Hill.0', 0, ifelse(q == 'Hill.1', 1, 2))) %>%
  ggplot(aes(x = q, y = value))+
  geom_line(stat = 'smooth', method = 'loess', aes(group = interaction(Host, UM), color = Host), 
              se = F, alpha = 0.5, size = 0.7)+
  geom_point(aes(shape = UM_s, fill = Host), size = 3, alpha = 1)+
  scale_shape_manual(values = c('S' = 22, 'L' = 21))+
  theme_classic()+
  scale_x_continuous(breaks = c(0, 1, 2))+
  ylab('Diversity')
  

# Beta diversity with betapart
PERDIVER_beta_taxo<-beta.pair.abund(PERDIVER_sitesXsps[,-c(1:3)])

PERDIVER_disper_taxo_total<-betadisper(PERDIVER_beta_taxo[[3]], group = PERDIVER_sitesXsps$Host)
PERDIVER_disper_taxo_nest<-betadisper(PERDIVER_beta_taxo[[1]], group = PERDIVER_sitesXsps$Host)
PERDIVER_disper_taxo_turn<-betadisper(PERDIVER_beta_taxo[[2]], group = PERDIVER_sitesXsps$Host)

Beta_Total_Data<-data.frame('PCoA1' = PERDIVER_disper_taxo_total$vectors[,1], 
                'PCoA2' = PERDIVER_disper_taxo_total$vectors[,2],
                'Host' = PERDIVER_sitesXsps$Host, 
                'UM_s' = PERDIVER_sitesXsps$UM_s,
                'Dist_cent' = PERDIVER_disper_taxo_total$distances,
                'Type' = 'Total') %>%
  bind_rows(data.frame('PCoA1' = PERDIVER_disper_taxo_nest$vectors[,1], 
                       'PCoA2' = PERDIVER_disper_taxo_nest$vectors[,2],
                       'Host' = PERDIVER_sitesXsps$Host, 
                       'UM_s' = PERDIVER_sitesXsps$UM_s,
                       'Dist_cent' = PERDIVER_disper_taxo_nest$distances,
                       'Type' = 'Nestedness')) %>%
  bind_rows(data.frame('PCoA1' = PERDIVER_disper_taxo_turn$vectors[,1], 
                       'PCoA2' = PERDIVER_disper_taxo_turn$vectors[,2],
                       'Host' = PERDIVER_sitesXsps$Host, 
                       'UM_s' = PERDIVER_sitesXsps$UM_s,
                       'Dist_cent' = PERDIVER_disper_taxo_turn$distances,
                       'Type' = 'Turnover'))

Beta_centroids<-as.data.frame(PERDIVER_disper_taxo_total$centroids[,1:2]) %>%
  bind_rows(as.data.frame(PERDIVER_disper_taxo_nest$centroids[,1:2])) %>%
  bind_rows(as.data.frame(PERDIVER_disper_taxo_turn$centroids[,1:2])) %>%
  mutate(Host = rownames(.)) %>%
  mutate(Host = substr(Host, 1, 6)) %>%
  mutate(Type = rep(c('Total', 'Nestedness', 'Turnover'), each = 6))


ggplot(Beta_Total_Data, aes(x = PCoA1, y = PCoA2))+
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_point(aes(fill = Host, shape = UM_s), size = 3, color = 'black')+
  geom_polygon(aes(group = interaction(Host, UM_s), color = Host), alpha = 0)+
  coord_fixed()+
  theme_bw()+
  scale_shape_manual(values = c('S' = 22, 'L' = 21))+
  facet_wrap(~Type)

mod<-glmer(Dist_cent~UM_s+(1|Host),
                family = 'gaussian', 
                data = Beta_Total_Data)
mod0<-glmer(Hill.0~UM_s+(1|Host),
           family = 'poisson', 
           data = PERDIVER_alpha_taxo)
mod1<-lmer(Hill.1~UM_s+(1|Host), 
            data = PERDIVER_alpha_taxo)
mod2<-lmer(Hill.2~UM_s+(1|Host), 
            data = PERDIVER_alpha_taxo)

summary(mod); visreg(mod)
summary(mod0); visreg(mod0)
car::Anova(mod0)
summary(mod1); visreg(mod1)
summary(mod2); visreg(mod2)
car::Anova(mod1)
car::Anova(mod2)

###############################
#### phylogenetic analysis ####
##############################

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

# Calculo de PD y su SES. Se hace por separado para usar un pool de espcies diferente 
# por cada especie porque si no saldria raro

pd_host<-list()
runs<-500
for (q in 1:length(unique(PERDIVER_sitesXsps_filo$Host))) {
  temp1<-PERDIVER_sitesXsps_filo %>%
    filter(Host == unique(PERDIVER_sitesXsps_filo$Host)[q]) %>%
    dplyr::select(-c(Host, UM, UM_s))
  
  obs_host<-popbio::mean.list(lapply(filo_perdiver, function(x) matrix(
    pd.query(
      tree = x, matrix = temp1),
    ncol = 1)))
  
  dist_host<-list()
  for (i in 1:runs) {
    temp2<-PERDIVER_sitesXsps_filo %>% 
      filter(Host == unique(PERDIVER_sitesXsps_filo$Host)[q]) %>% 
      dplyr::select(-c(Host, UM, UM_s)) %>% 
      randomizeMatrix(samp = ., null.model = 'richness')
    dist_host[[i]]<-popbio::mean.list(lapply(filo_perdiver, function(x) matrix(
      pd.query(
        tree = x, matrix = temp2),
      ncol = 1)))
  }
  pd_host[[q]]<-cbind(obs_host, do.call('cbind', dist_host))
  rm(dist_host, temp1, temp2, obs_host)
}


PERDIVER_pd<-lapply(pd_host, function(x) data.frame('obs.pd' = x[,1],
  'ses.pd' = (x[,1]-apply(x[,-1], 1, mean))/apply(x[,-1], 1, sd)))

PERDIVER_pd<-cbind(PERDIVER_sitesXsps_filo[,c(1:3)], do.call('rbind', PERDIVER_pd))

mod_p<-lm(Dist_cent~UM_s+(1|Host),
           family = 'gaussian', 
           data = PERDIVER_pd)
mod_p0<-lmer(ses.pd~UM_s+(1|Host), 
            data = PERDIVER_pd)

summary(mod_p0); visreg(mod_p0)
car::Anova(mod_p0)


PERDIVER_beta_filo<-phylo.beta.pair(PERDIVER_sitesXsps_filo[,-c(1:3)], tree = filo_perdiver[[1]])

PERDIVER_disper_filo_total<-betadisper(PERDIVER_beta_filo[[3]], group = PERDIVER_sitesXsps_filo$Host)
permutest(PERDIVER_disper_filo_total)

PERDIVER_disper_filo_nest<-betadisper(PERDIVER_beta_filo[[1]], group = PERDIVER_sitesXsps_filo$Host)
permutest(PERDIVER_disper_filo_nest)

PERDIVER_disper_filo_turn<-betadisper(PERDIVER_beta_filo[[2]], group = PERDIVER_sitesXsps_filo$Host)
permutest(PERDIVER_disper_filo_turn)


Beta_Total_Data_filo<-data.frame('PCoA1' = PERDIVER_disper_filo_total$vectors[,1], 
                            'PCoA2' = PERDIVER_disper_filo_total$vectors[,2],
                            'Host' = PERDIVER_sitesXsps$Host, 
                            'UM_s' = PERDIVER_sitesXsps$UM_s,
                            'Dist_cent' = PERDIVER_disper_filo_total$distances,
                            'Type' = 'Total') %>%
  bind_rows(data.frame('PCoA1' = PERDIVER_disper_filo_nest$vectors[,1], 
                       'PCoA2' = PERDIVER_disper_filo_nest$vectors[,2],
                       'Host' = PERDIVER_sitesXsps$Host, 
                       'UM_s' = PERDIVER_sitesXsps$UM_s,
                       'Dist_cent' = PERDIVER_disper_filo_nest$distances,
                       'Type' = 'Nestedness')) %>%
  bind_rows(data.frame('PCoA1' = PERDIVER_disper_filo_turn$vectors[,1], 
                       'PCoA2' = PERDIVER_disper_filo_turn$vectors[,2],
                       'Host' = PERDIVER_sitesXsps$Host, 
                       'UM_s' = PERDIVER_sitesXsps$UM_s,
                       'Dist_cent' = PERDIVER_disper_filo_turn$distances,
                       'Type' = 'Turnover'))

Beta_centroids_filo<-as.data.frame(PERDIVER_disper_filo_total$centroids[,1:2]) %>%
  bind_rows(as.data.frame(PERDIVER_disper_filo_nest$centroids[,1:2])) %>%
  bind_rows(as.data.frame(PERDIVER_disper_filo_turn$centroids[,1:2])) %>%
  mutate(Host = rownames(.)) %>%
  mutate(Host = substr(Host, 1, 6)) %>%
  mutate(Type = rep(c('Total', 'Nestedness', 'Turnover'), each = 6))


ggplot(Beta_Total_Data_filo %>% filter(Type == 'Total'), aes(x = PCoA1, y = PCoA2))+
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_point(aes(fill = Host, shape = UM_s), size = 3, color = 'black')+
  geom_polygon(aes(group = interaction(Host, UM_s), color = Host), alpha = 0)+
  coord_fixed()+
  theme_bw()+
  scale_shape_manual(values = c('S' = 22, 'L' = 21))


grid.arrange(ggplot(Beta_Total_Data %>% filter(Type == 'Total'), aes(x = PCoA1, y = PCoA2))+
            geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.3)+
            geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.3)+
            geom_point(aes(fill = Host, shape = UM_s), size = 3, color = 'black')+
            geom_polygon(aes(group = interaction(Host, UM_s), color = Host), alpha = 0)+
            coord_fixed()+
            theme_bw()+
            scale_shape_manual(values = c('S' = 22, 'L' = 21))+ggtitle('Taxo'),
          ggplot(Beta_Total_Data_filo %>% filter(Type == 'Total'), aes(x = PCoA1, y = PCoA2))+
            geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.3)+
            geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.3)+
            geom_point(aes(fill = Host, shape = UM_s), size = 3, color = 'black')+
            geom_polygon(aes(group = interaction(Host, UM_s), color = Host), alpha = 0)+
            coord_fixed()+
            theme_bw()+
            scale_shape_manual(values = c('S' = 22, 'L' = 21))+ggtitle('Filo'),
          nrow = 1)
