library(ape)
library(tidyverse)
library(picante)
library(caper)
library(taxize)
library(taxizedb)
library(rgbif)
library(V.PhyloMaker2)
rm(list = ls())

Sys.setenv(ENTREZ_KEY='7c76967d93d4938a1c715afc62fe297a8308')

# Cargamos el arbol actualizado
insect_tree<-read.tree('C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Datos/Arboles/Chesters_2022.bestTree')
insect_tree<-read.tree('~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Datos/Arboles/Chesters_2022.bestTree')

# Obtenemos los generos en el arbol
genera_Chesters<-data.frame(species = insect_tree$tip.label, chester_genus = unlist(lapply(strsplit(insect_tree$tip.label, '_'), function(x) x[1])))




### con gbif 
info_gbif_chesters<-name_backbone_checklist(
  data.frame(name = unique(genera_Chesters$chester_genus), class = 'Insecta'))

info_tips_Chester<-info_gbif_chesters %>%
  dplyr::select(c(verbatim_name, order, family, genus)) %>%
  rename(chester_genus = verbatim_name) %>%
  full_join(genera_Chesters, 'chester_genus') %>%
  dplyr::select(c(species, chester_genus, family, order)) %>%
  rename(genus = chester_genus)

missing<-info_tips_Chester %>%
  filter(!complete.cases(.)) 

Taxonomy_taxizedb<-sapply(unique(missing$genus), function(x) name2taxid(x = x, db = 'ncbi', out_type = 'summary', verbose = T))
Taxonomy_taxizedb<-as.data.frame(apply(Taxonomy_taxizedb, 1, unlist))
Taxonomy_taxizedb <- Taxonomy_taxizedb %>%
  filter(id != 451717) # Hay que quitar manualmente esta porque es una planta
orders_tree<-taxa_at(x = Taxonomy_taxizedb$id, rank = c('order'))
families_tree<-taxa_at(x = Taxonomy_taxizedb$id, rank = 'family')
orders_tree_df<-do.call('rbind', orders_tree)
families_tree_df<-do.call('rbind', families_tree)
taxonomy_df<-merge(orders_tree_df, families_tree_df, by = 0)
rownames(taxonomy_df)<-taxonomy_df$Row.names
names(taxonomy_df)[c(1, 2, 5)]<-c('genus', 'order', 'family')
taxonomy_df<-taxonomy_df %>%
  rename(id = genus) %>%
  dplyr::select(c(id, order, family)) %>%
  full_join(Taxonomy_taxizedb, 'id') %>%
  rename(genus = name) %>%
  dplyr::select(-id)

info_tips_Chester<-info_tips_Chester %>% 
  filter(complete.cases(.)) %>%
  bind_rows(missing %>%
              dplyr::select(species, genus) %>%
              full_join(taxonomy_df, 'genus')) %>%
  filter(species %in% insect_tree$tip.label) %>%
  distinct()

nodes_Chester<-build.nodes.1(tree = keep.tip(insect_tree, info_tips_Chester$species), tips = info_tips_Chester[,c(1:3)])
save(info_tips_Chester, nodes_Chester, file = 'Nodes_info_Chester2.RData')


############## Solo con NCBI
# Como son muchos y tarda demasiado, vamos a usar taxizedb, que es mas rapido. 
# Primero obtenemos los IDs de NCBI, porque algunos nombres pueden tener conflictos
Taxonomy_taxizedb<-sapply(genera_Chesters, function(x) name2taxid(x = x, db = 'ncbi', out_type = 'summary'))

Taxonomy_taxizedb<-as.data.frame(apply(Taxonomy_taxizedb, 1, unlist))

# Sacamos los que no tengan conflictos (un solo ID) y los que si (mas de 1)
no_conflicts<-Taxonomy_taxizedb %>%
  group_by(name) %>%
  filter(n() == 1) %>%
  pull(name)
  
conflicts<-Taxonomy_taxizedb %>%
  group_by(name) %>%
  filter(n() > 1) %>%
  pull(name)

# Obtenemos el genero y subfamilia de los que no tienen conflictos
orders_tree<-taxa_at(no_conflicts, rank = 'order')
families_tree<-taxa_at(no_conflicts, rank = 'family')
subfamilies_tree<-taxa_at(no_conflicts, rank = 'subfamily')

orders_tree_df<-do.call('rbind', orders_tree)
families_tree_df<-do.call('rbind', families_tree)
subfamilies_tree_df<-do.call('rbind', subfamilies_tree)

# Juntamos orden, familia y subfamilia
taxonomy_df<-merge(orders_tree_df, families_tree_df, by = 0)
rownames(taxonomy_df)<-taxonomy_df$Row.names
taxonomy_df<-taxonomy_df[,-1]
taxonomy_df<-merge(taxonomy_df, subfamilies_tree_df, by = 0)
names(taxonomy_df)[c(1, 2, 5, 8)]<-c('genus', 'order', 'family', 'subfamily')
taxonomy_df$subfamily<-ifelse(taxonomy_df$rank == 'subfamily', taxonomy_df$subfamily, NA)
taxonomy_df<-taxonomy_df %>%
  dplyr::select(c(genus, order, family, subfamily))

# Obtenemos manualmente los datos de las especies con conflictos
cosa_conflicts<-as.data.frame(sapply(unique(conflicts), function (x) tax_name(sci = x, get = c('order', 'family', 'subfamily'), 
                                                                  db = 'ncbi')))

conflicts_solved<-cosa_conflicts %>%
  unlist() %>%
  matrix(ncol = 5, byrow = T) %>%
  as.data.frame() %>%
  dplyr::select(-V1) %>%
  rename(genus = V2, order = V3, 
         family = V4, subfamily = V5)

# Juntamos todo
info_tips_Chester<-taxonomy_df %>%
  bind_rows(conflicts_solved) %>%
  inner_join(data.frame('species' = insect_tree$tip.label, 
                        'genus' = unlist(lapply(strsplit(insect_tree$tip.label, '_'), function(x) x[1]))),
             'genus')

# En algun momento se pierden algunas especies, lo corregimos
faltan<-info_tips_Chester %>%
  filter(is.na(order)) %>%
  dplyr::select(c(genus, species))

faltan_ncbi<-as.data.frame(sapply(unique(faltan$genus), function (x) tax_name(sci = x, get = c('order', 'family', 'subfamily'), 
                                                                              db = 'ncbi')))
faltan_solved<-faltan_ncbi %>%
  unlist() %>%
  matrix(ncol = 5, byrow = T) %>%
  as.data.frame() %>%
  dplyr::select(-V1) %>%
  rename(genus = V2, order = V3, 
         family = V4, subfamily = V5) %>%
  full_join(faltan, 'genus')

# Hay una especie que no sale en NCBI, la buscamos manualmente en GBIF y listo
faltan_solved[faltan_solved$genus == 'Harveyope',]$order<-'Lepidoptera'
faltan_solved[faltan_solved$genus == 'Harveyope',]$family<-'Riodinidae'

# Juntamos todo
info_tips_Chester<-info_tips_Chester %>%
  filter(!is.na(order)) %>%
  bind_rows(faltan_solved) %>%
  dplyr::select(c(species, genus, family, subfamily, order))

nodes_Chester<-build.nodes.1(tree = keep.tip(insect_tree, info_tips_Chester$species), tips = info_tips_Chester[,c(1:3)])

save(info_tips_Chester, nodes_Chester, file = 'Nodes_info_Chester.RData')
