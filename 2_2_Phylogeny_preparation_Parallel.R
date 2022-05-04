library(V.PhyloMaker)
library(picante)
library(tidyverse)
library(phangorn)
library(castor)
library(caper)
library(PDcalc)
library(doParallel)
library(parallel)
library(foreach)

load("C:/Users/18172844S/Dropbox/Insects phylogeny/Resultados/Datos_base_PERDIVER.RData")
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Nodes_info_Chester.RData")

load("~/Dropbox/Insects phylogeny/Resultados/Datos_base_PERDIVER.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Nodes_info_Chester.RData")

# Cargamos arboles con familias y subfamilias en los nodos y juntamos los nodos
insect_tree<-read.tree('C:/Users/18172844S/Dropbox/Insects phylogeny/Datos/Arboles/Chesters_2022.bestTree')
insect_tree<-read.tree('~/Dropbox/Insects phylogeny/Datos/Arboles/Chesters_2022.bestTree')
insect_tree<-chronoMPL(multi2di(insect_tree))

rm(tree.taxa.COMPLETO, families.in.db.COMPLETO, families.in.tree.COMPLETO,
   insectos.corrected)

# Solo los insectos
PERDIVER_insectos<-PERDIVER %>% 
  filter(Clase == 'Insecta')
rm(PERDIVER)

# Matriz de especies a injertar
matriz_Perdiver<-PERDIVER_insectos %>%
  dplyr::select(c(Species, Family_final, Order_final)) %>%
  mutate(genus = unlist(lapply(strsplit(Species, '_'), function(x) x[1]))) %>%
  dplyr::rename(species = Species, family = Family_final, order = Order_final) %>%
  dplyr::select(c(species, genus, family, order)) %>%
  filter(species != "Pseudisotoma_monochaeta") %>%
  distinct() 

matriz_completas<-matriz_Perdiver %>%
  filter(genus != family) %>%
  filter(genus != order)


# Montamos la filogenia provisional. Va a dar un mensaje de error, es importante porque tendremos
# que corregir esas cosas. En caso de conflicto usamos la taxonom?a de NCBI (la nuestra)

# Creamos una distribucion de arboles posibles de tamanyo r. Es importante
# establecer el parametro output.tree = T para que nos de el arbol completo
# sin podar. Asi podemos insertar el resto de taxa sin problema
n_trees<-1000
gc()

no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)

filo_perdiver<-foreach(K = 1:n_trees,
        .packages=c('V.PhyloMaker', 
                    'phangorn',
                    'tidyverse',
                    'castor',
                    'caper',
                    'PDcalc',
                    'picante')) %dopar% {
  temp_phylo<-phylo.maker(sp.list = matriz_completas[,1:3], 
                          tree = insect_tree, 
                          nodes = nodes_Chester,
                          scenarios = c('S2'),
                          output.tree = T,
                          r = 1,
                          output.sp.list = F)
  
  # Despues de obtener la distibucion de arboles, vamos a anyadir el resto de 
  # especies. Vamos a unirlas a un nodo elegido al azar dentro de su orden.
  # Elegimos el nodo en funcion del MRCA de todas las especies del orden
  # Primero sacamos el nodo donde se va a injertar cada especie y despues 
  # las anyadimos todas a la vez.
  
  arbol<-temp_phylo$tree.scenario.2$run.1
  rm(temp_phylo); gc()
  
  missing_sps_with_family<-matriz_Perdiver %>%
    filter(!species %in% arbol$tip.label) %>%
    filter(family %in% info_tips_Chester$family) %>%
    pull(species) %>%
    unique()
  
  for (i in 1:length(missing_sps_with_family)) {
    temp_family<-PERDIVER_insectos %>% 
      filter(Species == missing_sps_with_family[[i]]) %>%
      pull(Family_final) %>%
      unique(.)
    temp_tips<-unique(c(PERDIVER_insectos %>%
                          filter(Family_final == temp_family) %>%
                          filter(Species %in% arbol$tip.label) %>%
                          pull(Species),
                        info_tips_Chester %>%
                          filter(family == temp_family) %>%
                          mutate(species = as.character(species)) %>%
                          pull(species)))
    temp_MRCA<-mrca.phylo(x = arbol, node = temp_tips)
    desc_all<-unlist(Descendants(x = arbol, node = temp_MRCA, 'all'))
    desc_tips<-unlist(Descendants(x = arbol, node = temp_MRCA, 'tips'))
    nodes_tips<-which(arbol$tip.label %in% temp_tips)
    
    conflict_tips<-info_tips_Chester %>% filter(species %in% 
                                                  arbol$tip.label[desc_tips[!desc_tips %in% nodes_tips]]) %>%
      pull(species)
    
    conflict_nodes<-unlist(Descendants(x = arbol, 
                                       node = mrca.phylo(x = arbol, 
                                                         node = conflict_tips)))
    
    temp_node<-sample(desc_all[!desc_all %in% c(desc_tips, conflict_nodes)], 1)
    rm(desc_all, desc_tips, conflict_tips, conflict_nodes)
    arbol<-int.node(phylogeny = arbol,
                    location.node = temp_node,
                    tip.label = missing_sps_with_family[i],
                    position = runif(1))
    
  }
  
  missing_sps<-matriz_Perdiver %>%
    filter(!species %in% arbol$tip.label) %>%
    pull(species) %>%
    unique()
  
  for (i in 1:length(missing_sps)) {
    temp_order<-PERDIVER_insectos %>% 
      filter(Species == missing_sps[[i]]) %>%
      pull(Order_final) %>%
      unique(.)
    temp_tips<-unique(c(PERDIVER_insectos %>%
                          filter(Order_final == temp_order) %>%
                          filter(Species %in% arbol$tip.label) %>%
                          pull(Species),
                        info_tips_Chester %>%
                          filter(order == temp_order) %>%
                          mutate(species = as.character(species)) %>%
                          pull(species)))
    
    temp_node<-sample(temp_tips, 1)
    
    arbol<-ext.node(phylogeny = arbol,
                    location.tip = temp_node,
                    tip.label = missing_sps[i],
                    position = runif(1))
  }
  
  
  arbol<-multi2di(keep.tip(arbol, matriz_Perdiver$species))
  gc()
  arbol$edge.length<-arbol$edge.length+(2*abs(min(arbol$edge.length)))
  arbol<-chronos(arbol, lambda = 0, 'correlated', quiet = T)
  class(arbol)<-'phylo'
  arbol
}
stopCluster(cl)

save(filo_perdiver, file = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER.Rdata')

filo_perdiver_fam<-list()
filo_perdiver_ord<-list()
for (q in 1:length(filo_perdiver)) {
  tree_t<-data.frame('species' = filo_perdiver[[q]]$tip.label) %>% 
    left_join(., matriz_Perdiver, 'species') %>%
    mutate(family = ifelse(is.na(family) & 
                             order == 'Thysanoptera', 
                           'Thysanoptera',
                           family)) %>%
    mutate(fam_gen = paste(family, species, sep = '_'),
           ord_gen = paste(order, species, sep = '_'))
  
  arbol_fam<-filo_perdiver[[q]] 
  arbol_fam$tip.label<-tree_t$fam_gen
  
  arbol_ord<-filo_perdiver[[q]] 
  arbol_ord$tip.label<-tree_t$ord_gen
  
  y<-getGenus(arbol_fam)
  y<-y %>%
    filter(genus != 'NA' &
             mrca != 'NA' & 
             genus != 'Machilidae')
  
  z<-getGenus(arbol_ord)
  z<-z %>%
    filter(mrca > 455)
  
  filo_perdiver_fam[[q]]<-fix.poly(tree = filo_perdiver[[q]],
                                   type = 'collapse',
                                   node = y$mrca) 
  filo_perdiver_fam[[q]]$tip.label<-tree_t$species
  
  filo_perdiver_ord[[q]]<-fix.poly(tree = filo_perdiver[[q]], 
                                   type = 'collapse', 
                                   node = z$mrca)
  filo_perdiver_ord[[q]]$tip.label<-tree_t$species
  
}

save(filo_perdiver_fam, file = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER_Family.Rdata')
save(filo_perdiver_ord, file = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER_Order.Rdata')

save(filo_perdiver_fam, file = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER_Family.Rdata')
save(filo_perdiver_ord, file = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Plant population size interactions/Resultados/Filo_PERDIVER_Order.Rdata')
