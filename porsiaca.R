obs<-list()
for (p in 1:length(unique(PERDIVER$Host))) {
  temp<-PERDIVER  %>%
    dplyr::select(c(Host, UM, Species)) %>%
    filter(Species %in% filo_perdiver[[1]]$tip.label) %>%
    mutate(UM_s = ifelse(startsWith(UM, 'S'), 'S', 'L')) %>%
    filter(Host == unique(Host)[[p]]) %>%
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
                values_from = n,
                values_fill = 0) %>%
    ungroup()
  
  qs<-seq(0, 2, by = 0.1)
  b<-c()
  for (q in 1:length(qs)) {
    b[[q]]<-apply(do.call('cbind', lapply(filo_perdiver, 
                                          function(x) hill_phylo(temp[,-c(1:3)],
                                                                 tree = x,
                                                                 q = qs[[q]]))),
                  1, 
                  mean)
  }
  obs[[p]]<-do.call('cbind', b)
} 
PERDIVER_alpha_phylo_obs<-as.data.frame(cbind(PERDIVER_sitesXsps_filo[,c(1:3)], 
                                         do.call('rbind', obs)))
names(PERDIVER_alpha_phylo_obs)[-c(1:3)]<-as.character(qs)


ran<-list()
runs<-999
for (r in 1:runs) {
  
coso<-list()
for (p in 1:length(unique(PERDIVER$Host))) {
  temp<-PERDIVER  %>%
    dplyr::select(c(Host, UM, Species)) %>%
    filter(Species %in% filo_perdiver[[1]]$tip.label) %>%
    mutate(UM_s = ifelse(startsWith(UM, 'S'), 'S', 'L')) %>%
    filter(Host == unique(Host)[[p]]) %>%
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
                values_from = n,
                values_fill = 0) %>%
    ungroup()
  
  qs<-seq(0, 2, by = 0.1)
  b<-c()
  for (q in 1:length(qs)) {
    b[[q]]<-apply(do.call('cbind', lapply(filo_perdiver, 
                                          function(x) hill_phylo(randomizeMatrix(temp[,-c(1:3)], 'richness'),
                                                                 tree = x,
                                                                 q = qs[[q]]))), 1, mean)
  }
  
  coso[[p]]<-do.call('cbind', b)
  
}
ran[[r]]<-do.call('rbind', coso)
}

PERDIVER_alpha_phylo<-cbind(PERDIVER_alpha_phylo_obs,
                            (do.call('rbind', a)-mean.list(ran))/sd.list(ran))
names(PERDIVER_alpha_phylo)[25:45]<-paste('ses', names(PERDIVER_alpha_phylo)[4:24], sep = '_')


PERDIVER_alpha_phylo %>%
  select(c(Host:'2')) %>%
  pivot_longer(cols = c('0':'2'), 
               names_to = 'q') %>%
  mutate(q = as.numeric(q)) %>%
  ggplot(aes(x = q, y = value, color = Host))+
  geom_line(aes(group = interaction(Host, UM)))+
  geom_point(data = PERDIVER_alpha_phylo %>%
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
  ylab('Diversity')+
  theme(text = element_text(size = 20))

PERDIVER_alpha_phylo %>%
  select(c(Host:UM_s, ses_0:ses_2)) %>%
  pivot_longer(cols = c('ses_0':'ses_2'), 
               names_to = 'q') %>%
  mutate(q = as.numeric(unlist(lapply(strsplit(q, '_'), function(x) x[2])))) %>%
  ggplot(aes(x = q, y = value, color = Host))+
  geom_line(aes(group = interaction(Host, UM)))+
  geom_point(data = PERDIVER_alpha_phylo %>%
               pivot_longer(cols = c('ses_0':'ses_2'), 
                            names_to = 'q') %>%
               mutate(q = as.numeric(unlist(lapply(strsplit(q, '_'), function(x) x[2])))) %>%
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
  ylab('Diversity')+
  theme(text = element_text(size = 20))
