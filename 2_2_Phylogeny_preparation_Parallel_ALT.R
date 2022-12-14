library(V.PhyloMaker2)
library(picante)
library(tidyverse)
library(phangorn)
library(castor)
library(caper)
library(PDcalc)
library(doParallel)
library(parallel)
library(foreach)
library(phylocanvas)
rm(list = ls())

load("C:/Users/18172844S/Dropbox/Insects phylogeny/Resultados/Datos_base_PERDIVER.RData")
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/Nodes_info_Chester2.RData")

load("~/Dropbox/Insects phylogeny/Resultados/Datos_base_PERDIVER.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/Nodes_info_Chester2.RData")

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
  dplyr::select(c(Species, family, order)) %>%
  mutate(genus = unlist(lapply(strsplit(Species, '_'), function(x) x[1]))) %>%
  dplyr::rename(species = Species) %>%
  dplyr::select(c(species, genus, family, order)) %>%
  distinct() 

## Hay que correrlo y pararlo antes de que termine pero copiando 
# en un txt el warning 
# 'Taxonomic classification not consistent between sp.list and tree.'
# Ese warning indica que familias no cuadran entre especies del arbol
# y de nuestra lista, para que podamos corregirlo
#phylo.maker(sp.list = matriz_Perdiver[,1:3], 
            tree = insect_tree, 
            nodes = nodes_Chester,
            scenarios = c('S2'),
            output.tree = T,
            output.sp.list = F,
            r = 1)

warning <- read.table('C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Scripts/warning_phylomaker.txt')
warning <- read.table('~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Scripts/warning_phylomaker.txt')

matriz_Perdiver <- matriz_Perdiver %>%
  full_join(warning, 'genus') %>%
  mutate(family = ifelse(!is.na(family_in_tree),
                         family_in_tree,
                         family)) %>% 
  dplyr::select(c(species, genus, family, order))

n_trees<-1000
gc()

no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)
filo_perdiver<-foreach(K = 1:n_trees,
                       .packages=c('V.PhyloMaker2', 
                                   'phangorn',
                                   'tidyverse',
                                   'castor',
                                   'caper',
                                   'PDcalc',
                                   'picante')) %dopar% {
                                    temp_phylo<-phylo.maker(sp.list = matriz_Perdiver[,1:3], 
                                                             tree = insect_tree, 
                                                             nodes = nodes_Chester,
                                                             scenarios = c('S2'),
                                                             output.tree = T,
                                                             output.sp.list = F)
                                     
                                     # Despues de obtener la distibucion de arboles, vamos a anyadir el resto de 
                                     # especies. Vamos a unirlas a un nodo elegido al azar dentro de su orden.
                                     # Elegimos el nodo en funcion del MRCA de todas las especies del orden
                                     # Primero sacamos el nodo donde se va a injertar cada especie y despues 
                                     # las anyadimos todas a la vez.
                                     
                                     arbol<-temp_phylo$tree.scenario.2
                                     rm(temp_phylo); gc()
                                    
                                     missing_sps<-matriz_Perdiver %>%
                                       filter(!species %in% arbol$tip.label)
                                     
                                     for (i in 1:length(missing_sps)) {
                                       temp_order<-matriz_Perdiver %>% 
                                         filter(species == missing_sps$species[[i]]) %>%
                                         pull(order) %>%
                                         unique(.)
                                       temp_tips<-unique(c(matriz_Perdiver %>%
                                                             filter(order == temp_order) %>%
                                                             filter(species %in% arbol$tip.label) %>%
                                                             pull(species),
                                                           info_tips_Chester %>%
                                                             filter(order == temp_order) %>%
                                                             mutate(species = as.character(species)) %>%
                                                             pull(species)))
                                       
                                       temp_node<-sample(temp_tips, 1)
                                       
                                       arbol<-ext.node(phylogeny = arbol,
                                                       location.tip = temp_node,
                                                       tip.label = missing_sps[i],
                                                       position = runif(n = 1, min = 0.1, max = 0.9))
                                     }
                                     
                                     
                                     arbol<-multi2di(keep.tip(arbol, matriz_Perdiver$species))
                                     gc()
                                     arbol$edge.length<-arbol$edge.length+(2*abs(min(arbol$edge.length)))
                                     arbol<-chronoMPL(arbol)
                                   }
stopCluster(cl)

save(filo_perdiver, file = 'C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/Filo_PERDIVER_ALT.Rdata')
save(filo_perdiver, file = '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/Filo_PERDIVER_ALT.Rdata')


######## Alternativas al problema de las morfoespecies incompletas ######
#### 1 ####
# Asignarles una familia aleatoria de las que estan en el arbol y dentro de su orden
temp_phylo<-phylo.maker(sp.list = matriz_Perdiver[,1:3], 
                        tree = insect_tree, 
                        nodes = nodes_Chester,
                        scenarios = c('S2'),
                        output.tree = T,
                        output.sp.list = F,
                        r = 1)

arbol<-temp_phylo$tree.scenario.2

missing_sps<-matriz_Perdiver %>%
  filter(!species %in% arbol$tip.label)

for (i in 1:dim(missing_sps)[1]) {
  temp_order<-PERDIVER_insectos %>% 
    filter(Species == missing_sps$species[[i]]) %>%
    pull(order) %>%
    unique(.)
  missing_sps[i,]$family<-matriz_Perdiver %>%
    filter(order == temp_order) %>%
    filter(species %in% arbol$tip.label) %>%
    bind_rows(info_tips_Chester %>%
                filter(order == temp_order) %>%
                mutate(species = as.character(species))) %>%
    pull(family) %>%
    unique(.) %>%
    sample(., 1)
}

matriz_Perdiver<-missing_sps %>%
  bind_rows(matriz_Perdiver %>% 
              filter(!species %in% missing_sps$species))


temp_phylo2<-phylo.maker(sp.list = matriz_Perdiver[,1:3], 
                         tree = insect_tree, 
                         nodes = nodes_Chester,
                         scenarios = c('S2'),
                         output.tree = T,
                         output.sp.list = F)

arbol2<-temp_phylo2$tree.scenario.2
arbol2 <- multi2di(arbol2)
arbol2 <- keep.tip(arbol2, matriz_Perdiver$species)
phylocanvas(arbol2)

##### 2 ######
# Seleccionar un nodo al azar que descienda del MRCA del orden correspondiente
# (de momento no funciona)
morfo_<-matriz_Perdiver %>%
  filter(genus != family)

temp_phylo<-phylo.maker(sp.list = morfo_[,1:3], 
                        tree = insect_tree, 
                        nodes = nodes_Chester,
                        scenarios = c('S2'),
                        output.tree = T,
                        output.sp.list = F)

arbol<-temp_phylo$tree.scenario.2$run.1

morfo_familia<-matriz_Perdiver %>%
  filter(!species %in% arbol$tip.label) %>%
  filter(family %in% info_tips_Chester$family &
           genus == family)

all_tips<-bind_rows(matriz_Perdiver %>%
                      filter(species %in% arbol$tip.label),
                    info_tips_Chester %>%
                      mutate(species = as.character(species)))

for (i in 1:dim(morfo_familia)[1]){
  temp_family<-matriz_Perdiver %>% 
    filter(species == morfo_familia$species[[i]]) %>%
    pull(family) %>%
    unique(.)
  temp_sps <- all_tips %>% 
    filter(family == temp_family) %>%
    pull(species) %>%
    unique(.)
  desirednode <- get_mrca_of_set(arbol, temp_sps)
  desc <- unlist(Descendants(arbol, desirednode, 'all'))
  if(length(desc) > 1){
    location <- sample(desc[desc > desirednode],1)
  } else {
    location <- desc
  }
  
  if(length(temp_sps) > 1){
    arbol <-  int.node(phylogeny = arbol, location.node = location, 
                    tip.label = morfo_familia$species[[i]], 
                    node.label = NULL,
                    position = runif(1, min = 0.1, max = 0.9))} else {
    arbol <-  ext.node(phylogeny = arbol,
                           location.tip = location,
                           tip.label = morfo_familia$species[[i]],
                           node.label = NULL,
                           position = runif(1, min = 0.1, max = 0.9))
    }
}

missing_sps<-matriz_Perdiver %>%
  filter(!species %in% arbol$tip.label)
all_tips<-bind_rows(matriz_Perdiver %>%
                      filter(species %in% arbol$tip.label),
                    info_tips_Chester %>%
                      mutate(species = as.character(species)))

for (i in 1:dim(missing_sps)[1]){
  temp_order<-matriz_Perdiver %>% 
    filter(species == missing_sps$species[[i]]) %>%
    pull(order) %>%
    unique(.)
  temp_sps <- all_tips %>% 
    filter(order == temp_order) %>%
    pull(species) %>%
    unique(.)
  desirednode <- get_mrca_of_set(arbol, temp_sps)
  desc <- unlist(Descendants(arbol, desirednode, 'all'))
  if(length(desc) > 1){
    location <- sample(desc[desc > desirednode],1)
  } else {
    location <- desc
  }
  if(length(temp_sps) > 1){
    arbol <-  int.node(phylogeny = arbol, location.node = location, 
                       tip.label = missing_sps$species[[i]], 
                       node.label = NULL,
                       position = runif(1, min = 0.1, max = 0.9))} else {
                         arbol <-  ext.node(phylogeny = arbol,
                                            location.tip = location,
                                            tip.label = missing_sps$species[[i]],
                                            node.label = NULL,
                                            position = runif(1, min = 0.1, max = 0.9))
                       }
}
arbol<-multi2di(arbol)
arbol$edge.length<-arbol$edge.length+(2*abs(min(arbol$edge.length)))
arbol<-chronoMPL(arbol)
phylocanvas(ladderize(keep.tip(arbol, matriz_Perdiver$species)))
