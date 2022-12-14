library(vegan)
library(lme4)
library(visreg)
library(gridExtra)
library(picante)
library(usedist)
library(tidyverse)
library(popbio)
library(PDcalc)
library(readxl)
library(ggpubr)
library(hillR)
library(emmeans)
library(patchwork)
library(ggnewscale)
library(car)
library(ggtext)
library(ggpattern)
library(ecodist)
library(lmerTest)
rm(list = ls())

load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/Diversity_results.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/Diversity_results.RData")
PERDIVER_summary$Plot<-paste(PERDIVER_summary$Host, PERDIVER_summary$UM, sep = '_')
PERDIVER_summary<-PERDIVER_summary %>%
  mutate(`Cobertura_vegetal_(copas)` = as.numeric(ifelse(`Cobertura_vegetal_(copas)` == '0-100', 
                                              50, 
                                              `Cobertura_vegetal_(copas)`))) %>%
  mutate(across(13:22, 
                scale)) %>%
  mutate(across(24:30, 
                scale))

paleta_plantas<-c('cypcal' = '#1f78b4', 'pinlon' = '#e31a1c', 
                  'borpyr' = '#a6cee3', 'genlut' = '#33a02c',
                  'galniv' = '#b2df8a', 'rammyc' = '#fdbf6f',
                  'kracer' = '#fb9a99')


# Diferencias en alfa taxo entre poblaciones
mod0<-glmer(Above_richness~UM_s+(1|Host),
            family = 'poisson', 
            data = PERDIVER_summary)
mod1<-lmer(exp(Above_shannon)~UM_s+(1|Host), 
           data = PERDIVER_summary)
mod2<-lmer((1/(1-Above_simpson))~UM_s+(1|Host), 
           data = PERDIVER_summary)

# Filo community structure above
mod_p0<-lmer(Insect_mpd.obs.z~UM_s+(1|Host), 
             data = PERDIVER_summary)
Anova(mod_p0)

### Beta diversity
### PERMANOVA above ###
#primero indicamos el bloque (la planta) para que tenga en cuenta eso

perm <- how(nperm = 10000)
setBlocks(perm) <- with(PERDIVER_summary, Host)
PERDIVER_summary<-PERDIVER_summary[match(labels(PERDIVER_beta_taxo), PERDIVER_summary$Plot),]
perm_taxo_above<-adonis2(PERDIVER_beta_taxo~UM_s, data = PERDIVER_summary,
        permutations = perm)

# Filo beta
PERDIVER_summary<-PERDIVER_summary[match(labels(PERDIVER_beta_filo), PERDIVER_summary$Plot),]
perm_filo_above<-adonis2(PERDIVER_beta_filo~UM_s, data = PERDIVER_summary,
        permutations = perm)


# NMDS con las distancias
taxo_points_tot<-metaMDS(PERDIVER_beta_taxo,
                         try = 1000)

Beta_Total_taxo<-bind_rows(data.frame('NMDS1' = taxo_points_tot$points[,1], 
                                      'NMDS2' = taxo_points_tot$points[,2],
                                      'Host' = PERDIVER_summary$Host, 
                                      'UM' = PERDIVER_summary$UM,
                                      'UM_s' = PERDIVER_summary$UM_s,
                                      'Type' = 'Taxo',
                                      'Comp' = 'Total',
                                      'Location' = 'Above')) %>%
  left_join(., PERDIVER_summary %>% dplyr::select(Host, UM, UM_s, Insect_richness)) %>%
  dplyr::rename('Riq' = Insect_richness) %>%
  dplyr::mutate(group = as.factor(paste(Host, UM_s)))

plot.new()
eli_taxo_above<-ordiellipse(taxo_points_tot,
                            groups = Beta_Total_taxo$group,
                            display = 'sites',
                            kind = 'se',
                            conf = 0.95,
                            label = T)

df_ell_taxo_above<-do.call('rbind', lapply(eli_taxo_above, 
                                              function(x) vegan:::veganCovEllipse(x$cov, x$center, x$scale))) %>%
  as.data.frame() %>% 
  mutate(group = rep(names(eli_taxo_above), each = dim(.)[1]/length(eli_taxo_above))) %>%
  mutate(Host = unlist(lapply(strsplit(group, ' '), function(x) x[1])),
         UM_s = unlist(lapply(strsplit(group, ' '), function(x) x[2])))

Beta_taxo_above<-ggplot(Beta_Total_taxo, aes(x = NMDS1, y = NMDS2))+
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_path(data = df_ell_taxo_above, aes(group = interaction(Host, UM_s), color = Host))+
  geom_point(aes(fill = Host, shape = UM_s, size = Riq), color = 'black')+
  coord_fixed()+
  theme_classic()+
  scale_shape_manual(name = 'Isolation', values = c('S' = 22, 'L' = 21))+
  ggtitle('(a)')+
  xlim(c(-0.6, 0.6))+
  ylim(c(-0.5, 0.4))+
  scale_size_continuous(guide = 'none')+
  scale_color_manual(guide = 'none', values = paleta_plantas)+
  scale_fill_manual(values = paleta_plantas)+
  guides(fill=guide_legend(override.aes=list(shape=21, size = 5)),
         shape=guide_legend(override.aes=list(shape =c(0,1), size = 5, fill = 'black')))+
  geom_text(inherit.aes = F,
            aes(x = 0.45, y = -0.50), label = paste('F: ', round(perm_taxo_above$F[[1]], 2), ', ', 
                                                    'p: ', round(perm_taxo_above$`Pr(>F)`[[1]], 2),
                                                    ifelse(round(perm_taxo_above$`Pr(>F)`[[1]], 2) <= 0.05, '*', ''),
                                                    sep = ''),
            fontface = 'italic',
            size = 4)


### Filo above ###
filo_points_tot<-metaMDS(PERDIVER_beta_filo,
                         try = 100)

Beta_Total_Filo<-bind_rows(data.frame('NMDS1' = filo_points_tot$points[,1], 
                                      'NMDS2' = filo_points_tot$points[,2],
                                      'Host' = PERDIVER_summary$Host, 
                                      'UM' = PERDIVER_summary$UM,
                                      'UM_s' = PERDIVER_summary$UM_s,
                                      'Type' = 'Phylo',
                                      'Comp' = 'Total',
                                      'Location' = 'Above')) %>%
  left_join(., PERDIVER_summary %>% dplyr::select(Host, UM, UM_s, Insect_mpd.obs.z)) %>%
  mutate(group = as.factor(paste(Host, UM_s)))

plot.new()
eli_filo_above<-ordiellipse(filo_points_tot,
                            groups = Beta_Total_Filo$group,
                            display = 'sites',
                            kind = 'se',
                            conf = 0.95,
                            label = T)

df_ell_filo_above<-do.call('rbind', lapply(eli_filo_above, 
                                           function(x) vegan:::veganCovEllipse(x$cov, x$center, x$scale))) %>%
  as.data.frame() %>% 
  mutate(group = rep(names(eli_filo_above), each = dim(.)[1]/length(eli_filo_above))) %>%
  mutate(Host = unlist(lapply(strsplit(group, ' '), function(x) x[1])),
         UM_s = unlist(lapply(strsplit(group, ' '), function(x) x[2])))

Beta_filo_above<-ggplot(Beta_Total_Filo, aes(x = NMDS1, y = NMDS2))+
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_path(data = df_ell_filo_above, aes(group = interaction(Host, UM_s), color = Host))+
  geom_point(aes(fill = Host, shape = UM_s, size = Insect_mpd.obs.z), color = 'black')+
  theme_classic()+
  scale_shape_manual(name = 'Isolation', values = c('S' = 22, 'L' = 21))+
  ggtitle('(b)')+
  xlim(c(-0.6, 0.6))+
  ylim(c(-0.5, 0.4))+
  coord_fixed()+
  scale_size_continuous(guide = 'none')+
  scale_color_manual(guide = 'none', values = paleta_plantas)+
  scale_fill_manual(values = paleta_plantas)+
  guides(fill=guide_legend(override.aes=list(shape=21, size = 5)),
         shape=guide_legend(override.aes=list(shape =c(0,1), size = 5, fill = 'black')))+
  geom_text(inherit.aes = F,
            aes(x = 0.45, y = -0.50), label = paste('F: ', round(perm_filo_above$F[[1]], 2), ', ', 
                                                    'p: ', round(perm_filo_above$`Pr(>F)`[[1]], 2),
                                                    ifelse(round(perm_filo_above$`Pr(>F)`[[1]], 2) <= 0.05, '*', ''),
                                                    sep = ''),
            fontface = 'italic',
            size = 4)

##########################
###### Below ground ######
##########################
mod0_b<-lmer(Below_richness~UM_s+(1|Host), 
            data = PERDIVER_summary)
mod1_b<-lmer(exp(Below_shannon)~UM_s+(1|Host), 
           data = PERDIVER_summary)
mod2_b<-lmer(1/(1-Below_simpson)~UM_s+(1|Host), 
           data = PERDIVER_summary)

# Filo community structure below
mod_p0_b<-lmer(Below_mpd.obs.z~UM_s+(1|Host), 
             data = PERDIVER_summary)
Anova(mod_p0_b)

### Permanova below ####

perm <- how(nperm = 10000)
PERDIVER_summary<-PERDIVER_summary[match(labels(PERDIVER_beta_taxo_below), PERDIVER_summary$Plot),]
setBlocks(perm) <- with(PERDIVER_summary, Host)
perm_taxo_below<-adonis2(PERDIVER_beta_taxo_below~UM_s, data = PERDIVER_summary,
                         permutations = perm)

# Filo beta
PERDIVER_summary<-PERDIVER_summary[match(labels(PERDIVER_beta_filo_below), PERDIVER_summary$Plot),]
perm_filo_below<-adonis2(PERDIVER_beta_filo_below~UM_s, data = PERDIVER_summary,
                         permutations = perm)

### NMDS Taxo below ####
taxo_points_below<-metaMDS(PERDIVER_beta_taxo_below, try = 100)

Beta_Total_Taxo_Below<-data.frame(NMDS1 = taxo_points_below$points[,1],
                NMDS2 = taxo_points_below$points[,2],
                Type = 'Taxo',
                Comp = 'Total',
                Location = 'Below') %>%
  mutate(Host = unlist(lapply(strsplit(rownames(.), '_'), function(x) x[1]))) %>%
  mutate(UM = unlist(lapply(strsplit(rownames(.), '_'), function(x) x[2]))) %>%
  mutate(UM_s = substring(UM, 1, 1)) %>%
  left_join(., PERDIVER_summary %>% dplyr::select(Host, UM, UM_s, Below_richness)) %>%
  mutate(group = as.factor(paste(Host, UM_s)))

plot.new()
eli_taxo_below<-ordiellipse(taxo_points_below,
                            groups = Beta_Total_Taxo_Below$group,
                            display = 'sites',
                            kind = 'se',
                            conf = 0.95,
                            label = T)

df_ell_taxo_below<-do.call('rbind', lapply(eli_taxo_below, 
                                           function(x) vegan:::veganCovEllipse(x$cov, x$center, x$scale))) %>%
  as.data.frame() %>% 
  mutate(group = rep(names(eli_taxo_below), each = dim(.)[1]/length(eli_taxo_below))) %>%
  mutate(Host = unlist(lapply(strsplit(group, ' '), function(x) x[1])),
         UM_s = unlist(lapply(strsplit(group, ' '), function(x) x[2])))

Beta_taxo_below<-ggplot(Beta_Total_Taxo_Below, aes(x = NMDS1, y = NMDS2))+
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_path(data = df_ell_taxo_below, aes(group = interaction(Host, UM_s), color = Host))+
  geom_point(aes(fill = Host, shape = UM_s, size = Below_richness), color = 'black')+
  theme_classic()+
  scale_shape_manual(name = 'Isolation', values = c('S' = 22, 'L' = 21))+
  ggtitle('(c)')+
  xlim(c(-0.6, 0.6))+
  ylim(c(-0.5, 0.4))+
  coord_fixed()+
  scale_size_continuous(guide = 'none')+
  scale_color_manual(guide = 'none', values = paleta_plantas)+
  scale_fill_manual(values = paleta_plantas)+
  guides(fill=guide_legend(override.aes=list(shape=21, size = 5)),
         shape=guide_legend(override.aes=list(shape =c(0,1), size = 5, fill = 'black')))+
  geom_text(inherit.aes = F,
            aes(x = 0.45, y = -0.50), label = paste('F: ', round(perm_taxo_below$F[[1]], 2), ', ', 
                                                    'p: ', round(perm_taxo_below$`Pr(>F)`[[1]], 2),
                                                    ifelse(round(perm_taxo_below$`Pr(>F)`[[1]], 2) <= 0.05, '*', ''),
                                                    sep = ''),
            fontface = 'italic',
            size = 4)


# Filo below

filo_points_below<-metaMDS(PERDIVER_beta_filo_below, try = 100)

Beta_Total_Filo_Below<-data.frame(NMDS1 = filo_points_below$points[,1],
                            NMDS2 = filo_points_below$points[,2],
                            Type = 'Phylo',
                            Comp = 'Total',
                            Location = 'Below') %>%
  mutate(Host = unlist(lapply(strsplit(rownames(.), '_'), function(x) x[1]))) %>%
  mutate(UM = unlist(lapply(strsplit(rownames(.), '_'), function(x) x[2]))) %>%
  mutate(UM_s = substring(UM, 1, 1)) %>%
  left_join(., PERDIVER_summary %>% dplyr::select(Host, UM, UM_s, Below_mpd.obs.z))%>%
  mutate(group = as.factor(paste(Host, UM_s)))

plot.new()
eli_filo_below<-ordiellipse(filo_points_below,
                            groups = Beta_Total_Filo_Below$group,
                            display = 'sites',
                            kind = 'se',
                            conf = 0.95,
                            label = T)

df_ell_filo_below<-do.call('rbind', lapply(eli_filo_below, 
                                           function(x) vegan:::veganCovEllipse(x$cov, x$center, x$scale))) %>%
  as.data.frame() %>% 
  mutate(group = rep(names(eli_filo_below), each = dim(.)[1]/length(eli_filo_below))) %>%
  mutate(Host = unlist(lapply(strsplit(group, ' '), function(x) x[1])),
         UM_s = unlist(lapply(strsplit(group, ' '), function(x) x[2])))


Beta_filo_below<-ggplot(Beta_Total_Filo_Below, aes(x = NMDS1, y = NMDS2))+
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.3)+
  geom_path(data = df_ell_filo_below, aes(group = interaction(Host, UM_s), color = Host))+
  geom_point(aes(fill = Host, shape = UM_s, size = Below_mpd.obs.z), color = 'black')+
  theme_classic()+
  scale_shape_manual(name = 'Isolation', values = c('S' = 22, 'L' = 21))+
  ggtitle('(d)')+
  xlim(c(-0.6, 0.6))+
  ylim(c(-0.5, 0.4))+
  coord_fixed()+
  scale_size_continuous(guide = 'none')+
  scale_color_manual(guide = 'none', values = paleta_plantas)+
  scale_fill_manual(values = paleta_plantas)+
  guides(fill=guide_legend(override.aes=list(shape=21, size = 5)),
         shape=guide_legend(override.aes=list(shape =c(0,1), size = 5, fill = 'black')))+
  geom_text(inherit.aes = F,
            aes(x = 0.45, y = -0.50), label = paste('F: ', round(perm_filo_below$F[[1]], 2), ', ', 
                                                    'p: ', round(perm_filo_below$`Pr(>F)`[[1]], 2),
                                                    ifelse(round(perm_filo_below$`Pr(>F)`[[1]], 2) <= 0.05, '*', ''),
                                                    sep = ''),
            fontface = 'italic',
            size = 4)

NMDS_plot<-((Beta_taxo_above+Beta_filo_above)/(Beta_taxo_below+Beta_filo_below))+
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
ggsave(NMDS_plot, width = 10, height = 10, units = 'in',
       filename = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Figures and tables/Fig_2_NMDS.jpeg')
ggsave(NMDS_plot, width = 10, height = 10, units = 'in',
       filename = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Figures and tables/Fig_2_NMDS.jpeg')

#### Graficas diferencias ####

a<-summary(emmeans(mod0, specs = 'UM_s', type = 'response'))
b<-summary(emmeans(mod1, specs = 'UM_s'))
c<-summary(emmeans(mod2, specs = 'UM_s'))
d<-summary(emmeans(mod_p0, specs = 'UM_s'))

a_b<-summary(emmeans(mod0_b, specs = 'UM_s'))
b_b<-summary(emmeans(mod1_b, specs = 'UM_s'))
c_b<-summary(emmeans(mod2_b, specs = 'UM_s'))
d_b<-summary(emmeans(mod_p0_b, specs = 'UM_s'))

names(a)<-names(b)<-names(c)
names(a_b)<-names(b_b)<-names(c_b)

z<-rbind(a,b,c)
x<-rbind(a_b,b_b,c_b)
z$type<-rep(c('0', '1', '2'), each = 2)
x$type<-rep(c('0', '1', '2'), each = 2)

y<-rbind(d, d_b)
y$type<-rep(c('Above', 'Below'), each = 2)
y$type

paleta_tamanyos<-c('L' = 'forestgreen', 'S' = 'firebrick')

# Preparamos puntos para poner luego en los plots
points_above<-PERDIVER_summary %>%
  dplyr::select(c(Host, UM, UM_s, Above_richness, Above_shannon, Above_simpson)) %>%
  mutate(Above_shannon = exp(Above_shannon), Above_simpson = (1/(1-Above_simpson))) %>%
  rename('0' = Above_richness, '1' = Above_shannon, '2' = Above_simpson) %>%
  pivot_longer(cols = c('0', '1', '2'),
               names_to = 'type',
               values_to = 'emmean')

points_below<-PERDIVER_summary %>%
  dplyr::select(c(Host, UM, UM_s, Below_richness, Below_shannon, Below_simpson)) %>%
  mutate(Below_shannon = exp(Below_shannon), Below_simpson = (1/(1-Below_simpson))) %>%
  rename('0' = Below_richness, '1' = Below_shannon, '2' = Below_simpson) %>%
  pivot_longer(cols = c('0', '1', '2'),
               names_to = 'type',
               values_to = 'emmean')

points_mpd<-PERDIVER_summary %>%
  dplyr::select(c(Host, UM, UM_s, Insect_mpd.obs.z, Below_mpd.obs.z)) %>%
  rename(Above = Insect_mpd.obs.z, Below = Below_mpd.obs.z) %>%
  pivot_longer(cols = c(Above, Below),
               names_to = 'type',
               values_to = 'emmean')
# Plots
p1<-ggplot(z, aes(y = emmean, x = type, color = UM_s))+
  geom_jitter(data = points_above,
              aes(fill = Host, x = type, y = emmean, 
                  shape = UM_s, group = interaction(UM_s, type)), 
              inherit.aes = F,
              position = position_jitterdodge(jitter.width = .05, 
                                              jitter.height = .05,
                                              dodge.width = 0.2),
              size = 2, alpha = 0.65)+
  scale_fill_manual(values = paleta_plantas, guide = 'none')+
  scale_shape_manual(name = 'Isolation', values = c('S' = 22, 'L' = 21), guide = 'none')+
  geom_point(position = position_dodge(width = 0.2), size = 5)+
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL),
                 position = position_dodge(width = 0.2),
                 size = 1, width = 0.1)+
  geom_line(aes(x = type, y = emmean, group = UM_s),
            position = position_dodge(width = 0.2), linewidth = 1)+
  theme_classic()+
  ylab('Aboveground \u03B1 diversity')+
  xlab('q')+
  scale_color_manual(values = paleta_tamanyos, guide = 'none')+
  theme(text = element_text(size = 20))+
  ggtitle('(a)')

p2<-ggplot(y, aes(x = factor(type, level = c('Below', 'Above')), y = emmean, color = UM_s))+
  geom_jitter(data = points_mpd,
              aes(fill = Host, x = factor(type, level = c('Above', 'Below')), y = emmean, 
                  shape = UM_s, group = interaction(UM_s, type)), 
              inherit.aes = F,
              position = position_jitterdodge(jitter.width = .05, 
                                              jitter.height = .05,
                                              dodge.width = 0.4),
              size = 2, alpha = 0.65)+
  geom_point(size = 5, position = position_dodge(width = 0.4))+
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL), linewidth = 1,
                 position = position_dodge(width = 0.4),
                 width = 0.15)+
  scale_fill_manual(values = paleta_plantas)+
  scale_shape_manual(values = c('S' = 22, 'L' = 21), name = 'Isolation')+
  coord_fixed(ratio = .45)+
  theme_classic()+
  xlab('Location')+
  scale_color_manual(values = paleta_tamanyos, guide = 'none')+
  theme(text = element_text(size = 20),
        legend.position = 'bottom')+
  ylab(expression(MPD["SES"]))+
  ggtitle('(c)')+
  guides(fill=guide_legend(override.aes=list(shape=21, size = 5)),
         shape=guide_legend(override.aes=list(shape =c(0,1), 
                                              size = 5)))+
  scale_x_discrete(labels = c('Aboveground', 'Belowground'))


p3<-ggplot(x, aes(y = emmean, x = type, color = UM_s))+
  geom_jitter(data = points_below,
              aes(fill = Host, x = type, y = emmean, 
                  shape = UM_s, group = interaction(UM_s, type)), 
              inherit.aes = F,
              position = position_jitterdodge(jitter.width = .05, 
                                              jitter.height = .05,
                                              dodge.width = 0.2),
              size = 2, alpha = 0.65)+
  scale_fill_manual(values = paleta_plantas, guide = 'none')+
  scale_shape_manual(name = 'Isolation', values = c('S' = 22, 'L' = 21), guide = 'none')+
  geom_point(position = position_dodge(width = 0.2), size = 5)+
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL),
                 position = position_dodge(width = 0.2),
                 size = 1, width = 0.1)+
  geom_line(aes(x = type, y = emmean, group = UM_s),
            position = position_dodge(width = 0.2), size = 1)+
  theme_classic()+
  ylab('Belowground \u03B1 diversity')+
  xlab('q')+
  scale_color_manual(values = paleta_tamanyos, guide = 'none')+
  theme(text = element_text(size = 20))+
  ggtitle('(b)')

ANOVA_plot<-((p1/p3)|p2)+plot_layout(guides = 'collect')&theme(legend.position = 'bottom')
ggsave(ANOVA_plot, 
       filename = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Figures and tables/Fig_1_Differences.jpeg',
       width = 14, height = 9, units = 'in')

### Correlation between above-below
library(corrgram)
library(rstatix)
library(corrplot)

dat<-PERDIVER_summary %>%
  dplyr::select(c(Host, Above_richness, Above_shannon, Above_simpson,
                  Below_richness, Below_shannon, Below_simpson,
                  Insect_mpd.obs.z, Below_mpd.obs.z)) %>%
  mutate(Above_shannon = exp(Above_shannon),
         Above_simpson = (1/(1-Above_simpson)),
         Below_shannon = exp(Below_shannon),
         Below_simpson = (1/(1-Below_simpson)))

  
#### Correlograma simplificado ####
  b<-list(cor.test(dat$Above_richness, dat$Below_richness),
          cor.test(dat$Above_shannon, dat$Below_shannon),
          cor.test(dat$Above_simpson, dat$Below_simpson),
          cor.test(dat$Insect_mpd.obs.z, dat$Below_mpd.obs.z))
  c<-list(ecodist::mantel(lower(PERDIVER_beta_taxo)~lower(PERDIVER_beta_taxo_below)),
          ecodist::mantel(lower(PERDIVER_beta_filo)~lower(PERDIVER_beta_filo_below)))
  
  z<-data.frame(x = 1:6, y = 1, 
                stat = c(unlist(lapply(b, function(x) x$estimate)),
                         unlist(lapply(c, function(x) x[1]))),
                lower = c(unlist(lapply(b, function(x) x$conf.int[1])),
                          unlist(lapply(c, function(x) x[5]))),
                upper = c(unlist(lapply(b, function(x) x$conf.int[2])),
                          unlist(lapply(c, function(x) x[6]))),
                var1 =c('<sup>0</sup>D', '<sup>1</sup>D', '<sup>2</sup>D', 
                        'MPD<sub>SES</sub>', '\u03B2<sub>Taxo</sub>', '\u03B2<sub>Phylo</sub>'))
z$ast <- ifelse(z$lower < 0 & z$upper > 0, '', '***')
correlation_plot <- ggplot(z, aes(x = x, y = y+0.05))+
    geom_rect(aes(xmin = x-0.5, xmax = x+0.5, ymin = y-0.5, ymax = y+0.5, 
                  fill = stat, color = ast), alpha = 0.5, linewidth = 1.5)+
    geom_text(aes(label = round(stat, 3)), size = 15)+
    geom_text(aes(x = x, y = y-0.15, 
                  label = paste('[', round(lower, 3), ', ', round(upper, 3), ']', sep = '')),
              size = 7)+
    geom_richtext(aes(x = x, y = y+0.58, label = var1), size = 10,
                  fill = NA, label.color = NA)+
    geom_text(aes(x = x, y = y-0.3, label = ast), size = 10)+
    theme_void()+
    coord_fixed()+
    scale_fill_gradient(low = '#FFB000', high = '#648FFF',
                         name = 'Correlation')+
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, "cm"),
          legend.position = 'bottom',
          text = element_text(size = 20))+
    ylim(c(0.5, 1.8))+
  scale_color_manual(values = c('\\' = 'black', '***' = 'black'),
                     guide = 'none')
ggsave(correlation_plot, 
       filename = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Figures and tables/Fig_4_Correlation.jpeg',
       width = 14, height = 9, units = 'in')


### Correlograma completo
a<-list(mantel(PERDIVER_beta_taxo, PERDIVER_beta_taxo, 'pearson'),
        mantel(PERDIVER_beta_taxo, PERDIVER_beta_filo, 'pearson'),
        mantel(PERDIVER_beta_taxo, PERDIVER_beta_taxo_below, 'pearson'),
        mantel(PERDIVER_beta_taxo, PERDIVER_beta_filo_below, 'pearson'),
        
        mantel(PERDIVER_beta_filo, PERDIVER_beta_filo, 'pearson'),
        mantel(PERDIVER_beta_filo, PERDIVER_beta_taxo, 'pearson'),
        mantel(PERDIVER_beta_filo, PERDIVER_beta_taxo_below, 'pearson'),
        mantel(PERDIVER_beta_filo, PERDIVER_beta_filo_below, 'pearson'),
        
        mantel(PERDIVER_beta_taxo_below, PERDIVER_beta_taxo_below, 'pearson'),
        mantel(PERDIVER_beta_taxo_below, PERDIVER_beta_taxo, 'pearson'),
        mantel(PERDIVER_beta_taxo_below, PERDIVER_beta_filo, 'pearson'),
        mantel(PERDIVER_beta_taxo_below, PERDIVER_beta_filo_below, 'pearson'),
        
        mantel(PERDIVER_beta_filo_below, PERDIVER_beta_filo_below, 'pearson'),
        mantel(PERDIVER_beta_filo_below, PERDIVER_beta_taxo, 'pearson'),
        mantel(PERDIVER_beta_filo_below, PERDIVER_beta_taxo_below, 'pearson'),
        mantel(PERDIVER_beta_filo_below, PERDIVER_beta_filo, 'pearson'))

b<-data.frame(var1 = unlist(lapply(a, function(x) as.character(x$call$xdis))),
              var2 = unlist(lapply(a, function(x) as.character(x$call$ydis))))
b$cor<-unlist(lapply(a, function(x) as.numeric(x$statistic)))
b$statistic<-ifelse(b$var1 == b$var2, Inf, 0.00000000001)
b$p<-unlist(lapply(a, function(x) as.numeric(x$signif)))
b$conf.low<-ifelse(b$var1 == b$var2, 1, 0.00000000001)
b$conf.high<-ifelse(b$var1 == b$var2, 1, 0.00000000001)
b$method<-'pearson'

c<-rbind(dat, b)
varnames<-c('$Above*phantom(0)^0*D', '$Above*phantom(0)^1*D', '$Above*phantom(0)^2*D', 
            '$Below*phantom(0)^0*D', '$Below*phantom(0)^1*D', '$Below*phantom(0)^2*D',
            '$Above*phantom(0)[MPD]*ses', '$Below*phantom(0)[MPD]*ses',
            '$Above~beta[taxo]', '$Above~beta[phylo]',
            '$Below~beta[taxo]', '$Below~beta[phylo]')
d<-cor_spread(c)
d<-as.matrix(d[,-1]) 
colnames(d)<-varnames
rownames(d)<-colnames(d)

p<-cor_spread(c, value = 'p')
p<-as.matrix(p[,-1])
colnames(p)<- varnames
rownames(p)<-colnames(p)

trace(corrplot, edit=TRUE)
# Then replace on line 443
# place_points = function(sig.locs, point) {
# text(pos.pNew[, 1][sig.locs], (pos.pNew[, 2][sig.locs])-0.25, 
#     labels = point, col = pch.col, cex = pch.cex, 
#     lwd = 2)

jpeg(filename = 'Correlation_plot.jpeg', res = 300, width = 30, height = 30, units = 'cm')
corrplot(corr = d,
         method = 'color',
         type = 'lower',
         diag = F,
         na.label.col = 'white',
         addCoef.col = 'black',
         outline = T,
         p.mat = p,
         insig = 'label_sig',
         pch.cex = 2.5,
         number.digits = 2,
         tl.col = 'black',
         pch.col = 'black',
         number.cex = 1.4,
         tl.srt = 45,
         tl.cex = 1.5)
dev.off()


#### Prueba con distancias de centroides ####
temp<-do.call('rbind', lapply(list(PERDIVER_beta_taxo, PERDIVER_beta_filo,
                                   PERDIVER_beta_taxo_below, PERDIVER_beta_filo_below),
                              function(x) dist_groups(x, g = PERDIVER_summary$UM_s)))
Data_dist<-temp %>%
  mutate(Type = rep(c('Taxo', 'Phylo'), each = 378, times = 2),
         Location = rep(c('Above', 'Below'), each = 2*378)) %>%
  mutate(Host = unlist(lapply(strsplit(Item1, '_'), function(x) x[1])),
         UM = unlist(lapply(strsplit(Item1, '_'), function(x) x[2])),
         Host2 = unlist(lapply(strsplit(Item2, '_'), function(x) x[1])),
         UM2 = unlist(lapply(strsplit(Item2, '_'), function(x) x[2]))) %>%
  filter(Host == Host2) %>% 
  mutate(Groups = paste(UM, UM2, sep = '_')) %>% 
  select(-c(Item1, Item2, Group1, Group2, Host2, UM2)) %>%
  mutate(Type = factor(Type, levels = c('Taxo', 'Phylo'))) 


Data_dist_mean<-Data_dist %>%
  group_by(Host, Label, Type, Location) %>%
  summarise(Host = Host,
            Label = Label,
            Dist = mean(Distance),
            SD = sd(Distance), 
            SE = 1.96*SD/sqrt(n())) %>%
  distinct()

pp<-(ggplot(Data_dist_mean %>% 
              filter(Location == 'Above'), 
            aes(x = interaction(Label, Host), y = Dist))+
       geom_bar_pattern(stat = 'identity', aes(pattern_density = Label, pattern = Label, fill = Host), color = 'black')+
    geom_errorbar(aes(ymax = Dist+SE,
                      ymin = Dist-SE),
                  width = 0.3)+
    theme_classic2()+
    scale_x_discrete(position = 'bottom')+
    scale_fill_manual(values = paleta_plantas)+
    theme(axis.title = element_blank(),
          legend.position = 'none',
          axis.text.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())+
    facet_wrap(~Type, strip.position = 'bottom'))/
  (ggplot(Data_dist_mean %>% 
            filter(Location == 'Below'), 
          aes(x = interaction(Label, Host), y = Dist))+
  geom_bar_pattern(stat = 'identity', aes(pattern_density = Label, pattern = Label, fill = Host), color = 'black')+
  geom_errorbar(aes(ymax = Dist+SE,
                    ymin = Dist-SE),
                width = 0.3)+
  theme_classic2()+
  scale_y_reverse()+
  scale_x_discrete(position = 'top')+
  scale_fill_manual(values = paleta_plantas)+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank())+
  facet_wrap(~Type))
pp

dist_AT<-lmer(Distance~Label+(1|Host), data = (Data_dist %>% 
          filter(Location == 'Above' & Type == 'Taxo')))

dist_AP<-lmer(Distance~Label+(1|Host), data = (Data_dist %>% 
          filter(Location == 'Above' & Type == 'Phylo')))

dist_BT<-lmer(Distance~Label+(1|Host), data = (Data_dist %>% 
          filter(Location == 'Below' & Type == 'Taxo')))

dist_BP<-lmer(Distance~Label+(1|Host), data = (Data_dist %>% 
          filter(Location == 'Below' & Type == 'Phylo')))

dist_all <- rbind(ls_means(dist_AT), ls_means(dist_AP),
      ls_means(dist_BT), ls_means(dist_BP))
dist_all$Type <- rep(c('Taxo', 'Phylo'), each = 2)
dist_all$Location <- rep(c('Above', 'Above', 'Below', 'Below'), each = 2)
dist_all$Label <- rep(str_remove(rownames(dist_all)[1:2], 'Label'), times = 4)

q<-ggplot(dist_all, aes(y = Estimate, x = Type))+
  theme_classic()+
  geom_point(aes(group = interaction(Location, Label), color = Label),
             position = position_dodge(width = 0.2),
             size = 5)+
  geom_errorbar(aes(ymax = upper, ymin = lower,
                     group = interaction(Location, Label),
                    color = Label),
                 position = position_dodge(width = 0.2),
                 linewidth = 1,
                width = .05)+
  facet_wrap(~Location, ncol = 1)+
  scale_color_manual(values = c('Between L and S' = 'firebrick',
                                'Within L' = 'forestgreen'))+
  ylab('Average dissimilarity')+
  theme(text = element_text(size = 20))
q
confint()
