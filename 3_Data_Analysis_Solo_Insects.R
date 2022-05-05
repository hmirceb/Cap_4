library(betapart)
library(vegan)
library(lme4)
library(rstanarm)
library(visreg)
library(gridExtra)
library(PhyloMeasures)
library(picante)
library(usedist)
library(tidyverse)
library(popbio)
library(PDcalc)
library(readxl)
library(ggpubr)
library(hillR)
library(emmeans)

load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER.Rdata")
Habitat_PERDIVER<-read_excel("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Datos/Plantas/Habitat_PERDIVER.xlsx")

load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER.Rdata")
Habitat_PERDIVER<-read_excel("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Datos/Plantas/Habitat_PERDIVER.xlsx")


# Funcion para calcular SD de una lista de matrices
sd.list <- function(lst) {n <- length(lst)
  rc <- dim(lst[[1]]) 
  ar1 <- array(unlist(lst), c(rc, n)) 	   
  apply(ar1, c(1, 2), sd)
  }

# Distancia filo entre especies media (de todos los arboles)
cophenetic_mean<-mean.list(lapply(filo_perdiver, cophenetic))

# Matriz de sitios por especies
PERDIVER<-PERDIVER %>%
  filter(Host != 'kracer')

PERDIVER_sitesXsps_filo<-PERDIVER %>%
  dplyr::select(c(Host, UM, Species)) %>%
  filter(Species %in% filo_perdiver[[1]]$tip.label) %>%
  mutate(UM_s = ifelse(startsWith(UM, 'S'), 'S', 'L')) %>%
  group_by(Host, UM, Species, UM_s) %>%
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

# Matriz de sitios por especies para las plantas coocurrentes
PLANTS_sitesXsps<-Habitat_PERDIVER %>%
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

plot_hill_insects<-PERDIVER_alpha_taxo %>%
            pivot_longer(cols = c('0':'2'), 
                         names_to = 'q') %>%
            mutate(q = as.numeric(q)) %>%
            ggplot(aes(x = q, y = value, color = Host))+
            geom_line(aes(group = interaction(Host, UM)))+
            geom_point(data = PERDIVER_alpha_taxo %>%
                         pivot_longer(cols = c('0':'2'), 
                                      names_to = 'q') %>%
                         mutate(q = as.numeric(q)) %>%
                         filter(q %in% c(0, 1, 2)),
                       aes(x = q, y = value, 
                           fill = Host, 
                           shape = UM_s), 
                       color = 'black', 
                       size = 3,
                       inherit.aes = F)+
            scale_shape_manual(guide = 'none', values = c('S' = 22, 'L' = 21))+
            theme_classic()+
            scale_x_continuous(breaks = c(0, 1, 2))+
            guides(fill = guide_legend(override.aes=list(shape=21)))+
            ylab('Taxonomic diversity')+
            theme(text = element_text(size = 20))+
  ggtitle('A')

plot_hill_plants<-PERDIVER_alpha_taxo %>%
            dplyr::select(-c('0':'2')) %>%
            pivot_longer(cols = c('0_plants':'2_plants'), 
                         names_to = 'q') %>%
            mutate(q = as.numeric(unlist(lapply(strsplit(q, '_'), function(x) x[1])))) %>%
            ggplot(aes(x = q, y = value, color = Host))+
            geom_line(aes(group = interaction(Host, UM)))+
            geom_point(data = PERDIVER_alpha_taxo %>%
                         pivot_longer(cols = c('0_plants':'2_plants'), 
                                      names_to = 'q') %>%
                         mutate(q = as.numeric(unlist(lapply(strsplit(q, '_'), function(x) x[1])))) %>%
                         filter(q %in% c(0, 1, 2)),
                       aes(x = q, y = value, 
                           fill = Host, 
                           shape = UM_s), 
                       color = 'black', 
                       size = 3,
                       inherit.aes = F)+
            scale_shape_manual(guide = 'none', values = c('S' = 22, 'L' = 21))+
            theme_classic()+
            scale_x_continuous(breaks = c(0, 1, 2))+
            guides(fill = guide_legend(override.aes=list(shape=21)))+
            ylab('Plant diversity')+
            theme(text = element_text(size = 20))

# Diferencias en alfa taxo entre poblaciones

mod0<-glmer(`0`~UM_s+(1|Host),
            family = 'poisson', 
            data = PERDIVER_alpha_taxo)
mod1<-lmer(`1`~UM_s+(1|Host), 
           data = PERDIVER_alpha_taxo)
mod2<-lmer(`2`~UM_s+(1|Host), 
           data = PERDIVER_alpha_taxo)

summary(mod0); visreg(mod0)
emmeans(mod0, specs = 'UM_s', type = 'response')
summary(mod1); visreg(mod1)
emmeans(mod1, specs = 'UM_s', type = 'response')
summary(mod2); visreg(mod2)
emmeans(mod2, specs = 'UM_s', type = 'response')

car::Anova(mod1)
car::Anova(mod2)

### TAXO Beta diversity
# Distancia con bray curtis (considera las abundancias)
PERDIVER_beta_taxo<-vegdist(PERDIVER_sitesXsps_filo[,-c(1:3)],
                            method = 'bray')

# NMDS con las distancias
taxo_points_tot<-metaMDS(PERDIVER_beta_taxo,
                         try = 100)$points

Beta_Total_taxo<-bind_rows(data.frame('NMDS1' = taxo_points_tot[,1], 
                                      'NMDS2' = taxo_points_tot[,2],
                                      'Host' = PERDIVER_sitesXsps_filo$Host, 
                                      'UM_s' = PERDIVER_sitesXsps_filo$UM_s,
                                      'Type' = 'Taxo',
                                      'Comp' = 'Total'))

ggplot(Beta_Total_taxo, aes(x = NMDS1, y = NMDS2))+
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_point(aes(fill = Host, shape = UM_s), size = 3, color = 'black')+
  geom_polygon(aes(group = interaction(Host, UM_s), color = Host), alpha = 0)+
  coord_fixed()+
  theme_bw()+
  scale_shape_manual(values = c('S' = 22, 'L' = 21))+
  ggtitle('Taxonomic')


#### Beta diferencias L-S vs L-L ####

Taxo_dist_groups<-dist_groups(PERDIVER_beta_taxo, g = PERDIVER_sitesXsps_filo$UM_s) %>% 
  rename(UM_s = Label) %>%
  dplyr::select(c(UM_s)) %>%
  bind_cols(dist_groups(PERDIVER_beta_taxo, g = PERDIVER_sitesXsps_filo$Host)) %>%
  filter(!startsWith(as.character(Label), 'Between'))

mod<-lmer(Distance~UM_s+(1|Label), Taxo_dist_groups)
summary(mod)
visreg(mod)
car::Anova(mod)

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
  
  temp_mpd[[i]]<-ses.mpd(samp = as.data.frame(temp[,-c(1:3)]),
                         dis = cophenetic_mean, 
                         abundance.weighted = T)
}

PERDIVER_alpha_filo<-as.data.frame(cbind(PERDIVER_sitesXsps_filo[,c(1:3)], 
                                         do.call('rbind', temp_mpd)))

mod_p0<-lmer(-mpd.obs.z~UM_s+(1|Host), 
             data = PERDIVER_alpha_filo)

summary(mod_p0); visreg(mod_p0)
car::Anova(mod_p0)

plot_hill_insects_filo<-PERDIVER_alpha_filo %>%
  pivot_longer(cols = c('0':'2'), 
               names_to = 'q') %>%
  mutate(q = as.numeric(q)) %>%
  ggplot(aes(x = q, y = value, color = Host))+
  geom_line(aes(group = interaction(Host, UM)))+
  geom_point(data = PERDIVER_alpha_filo %>%
               pivot_longer(cols = c('0':'2'), 
                            names_to = 'q') %>%
               mutate(q = as.numeric(q)) %>%
               filter(q %in% c(0, 1, 2)),
             aes(x = q, y = value, 
                 fill = Host, 
                 shape = UM_s), 
             color = 'black', 
             size = 3,
             inherit.aes = F)+
  scale_shape_manual(guide = 'none', values = c('S' = 22, 'L' = 21))+
  theme_classic()+
  scale_x_continuous(breaks = c(0, 1, 2))+
  guides(fill = guide_legend(override.aes=list(shape=21)))+
  ylab('Phylogenetic diversity')+
  theme(text = element_text(size = 20))+
  ggtitle('B')

ggarrange(plot_hill_insects, plot_hill_insects_filo,
          ncol = 2, legend = 'bottom',
          common.legend = T)


### Beta diveristy with PDcalc
PERDIVER_beta_filo_list<-lapply(filo_perdiver, function(y) phyloresembl(x = PERDIVER_sitesXsps_filo[,-c(1:3)],
                                                                              phy = y, incidence = F, 'sorensen'))

PERDIVER_beta_filo<-unname(as.dist(mean.list(lapply(PERDIVER_beta_filo_list, function(x) as.matrix(x)))))

filo_points_tot<-metaMDS(PERDIVER_beta_filo,
                         try = 100)$points

Beta_Total_Filo<-bind_rows(data.frame('NMDS1' = filo_points_tot[,1], 
                                      'NMDS2' = filo_points_tot[,2],
                                      'Host' = PERDIVER_sitesXsps_filo$Host, 
                                      'UM_s' = PERDIVER_sitesXsps_filo$UM_s,
                                      'Type' = 'Phylo',
                                      'Comp' = 'Total'))


Filo_dist_groups<-dist_groups(PERDIVER_beta_filo, g = PERDIVER_sitesXsps_filo$UM_s) %>% 
  rename(UM_s = Label) %>%
  dplyr::select(c(UM_s)) %>%
  bind_cols(dist_groups(PERDIVER_beta_filo, g = PERDIVER_sitesXsps_filo$Host)) %>%
  filter(!startsWith(as.character(Label), 'Between'))

mod<-lmer(Distance~UM_s+(1|Label), Filo_dist_groups)
summary(mod)
visreg(mod)
car::Anova(mod)

### Plot de Betas con NMDS
ggarrange(ggplot(Beta_Total_taxo, aes(x = NMDS1, y = NMDS2))+
            geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.3)+
            geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.3)+
            geom_point(aes(fill = Host, shape = UM_s), size = 3, color = 'black')+
            geom_polygon(aes(group = interaction(Host, UM_s), color = Host), alpha = 0)+
            theme_classic()+
            scale_shape_manual(values = c('S' = 22, 'L' = 21))+
            ggtitle('Taxonomic beta')+
            xlim(c(-0.6, 0.6))+
            ylim(c(-0.2, 0.3)),
          ggplot(Beta_Total_Filo, aes(x = NMDS1, y = NMDS2))+
            geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.3)+
            geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.3)+
            geom_point(aes(fill = Host, shape = UM_s), size = 3, color = 'black')+
            geom_polygon(aes(group = interaction(Host, UM_s), color = Host), alpha = 0)+
            theme_classic()+
            scale_shape_manual(values = c('S' = 22, 'L' = 21))+
            ggtitle('Phylogenetic beta')+
            xlim(c(-0.6, 0.6))+
            ylim(c(-0.2, 0.3)),
          common.legend = T, legend = 'bottom')

#### Son las distancias taxo y filo mayores? ####
Total_dists_groups<-Taxo_dist_groups %>%
  mutate(across(Item1:Item2, as.character)) %>%
  mutate(Type = 'Taxo') %>%
  bind_rows(Filo_dist_groups %>% 
              mutate(Type = 'Filo') %>% 
              mutate(across(Item1:Item2, as.character)))

mod1<-lmer(Distance~UM_s*Type+(1|Label), Total_dists_groups)
DHARMa::simulateResiduals(mod1, plot = T)
summary(mod1)
car::Anova(mod1)

Dist_summary<-summary(emmeans(mod1, pairwise~Type | UM_s))
dist_contrast<-data.frame(Type = Dist_summary$emmeans$Type,
                          UM_s = Dist_summary$emmeans$UM_s,
                          emmean = Dist_summary$emmeans$emmean,
                          SE = Dist_summary$emmeans$SE,
                          df = Dist_summary$emmeans$df,
                          lower.CL = Dist_summary$emmeans$lower.CL,
                          upper.CL = Dist_summary$emmeans$upper.CL)

ggplot(dist_contrast, aes(x = Type, y = emmean))+
  geom_point(aes(group = UM_s, color = UM_s),
             position = position_dodge(width = 0.5),
             size = 3)+
  geom_errorbar(aes(group = UM_s, color = UM_s,
                    ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.5),
                width = 0.07,
                size = 1.1)+
  theme_classic()+
  ylab('Dissimilarity')+
  xlab('Type')+
  theme(text = element_text(size = 20))+
  geom_jitter(inherit.aes = F,
              data = Total_dists_groups,
              aes(x = Type,
                  y = Distance, 
                  color = UM_s,
                  group = UM_s),
              position = position_dodge(width = 0.5),
              alpha = 0.3)


#### Probamos a calcular la pendiente de la curva que hace los numeros de hill
#### como medida de evenness
PERDIVER_alpha_taxo$Even_10<-PERDIVER_alpha_taxo$`2`/PERDIVER_alpha_taxo$`1`

mod2<-lmer(Even_10~UM_s+(1|Host), PERDIVER_alpha_taxo)
DHARMa::simulateResiduals(mod2, plot = T)
summary(mod2)
car::Anova(mod2)
visreg(mod2)


#### Las poblaciones pequenyas interactuan mas quelas grandes? ####
PERDIVER_alpha_taxo$N<-apply(PERDIVER_sitesXsps_filo[,-c(1:3)], 1, sum)

mod3<-glmer(N~UM_s+(1|Host), family = 'poisson', PERDIVER_alpha_taxo)
DHARMa::simulateResiduals(mod3, plot = T)
summary(mod3)
car::Anova(mod3)
visreg(mod3)

PERDIVER_alpha_taxo$Prop<-PERDIVER_alpha_taxo$N/PERDIVER_alpha_taxo$`0`
mod4<-lmer(Prop~UM_s+(1|Host), PERDIVER_alpha_taxo)
DHARMa::simulateResiduals(mod4, plot = T)
summary(mod4)
car::Anova(mod4)
visreg(mod4)
