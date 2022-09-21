library(ape)
library(tidyverse)
library(picante)
library(caper)
library(taxize)
library(taxizedb)

Sys.setenv(ENTREZ_KEY='7c76967d93d4938a1c715afc62fe297a8308')

# Cargamos el arbol actualizado
insect_tree<-read.tree('C:/Users/18172844S/Dropbox/Insects phylogeny/Datos/Arboles/Chesters_2022.bestTree')
insect_tree<-read.tree('~/Dropbox/Insects phylogeny/Datos/Arboles/Chesters_2022.bestTree')

# Obtenemos los generos en el arbol
genera_Chesters<-unique(unlist(lapply(strsplit(insect_tree$tip.label, '_'), function(x) x[1])))

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
