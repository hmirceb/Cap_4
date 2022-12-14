library(tidyverse)
library(sf)
library(rgdal)
library(raster)
library(readxl)
library(usedist)

PERDIVER_abiotic <- read_excel("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Datos/PERDIVER_abiotic.xlsx")
sp<-PERDIVER_abiotic %>%
  dplyr::select(Host, UM, X, Y, Sist_Coord) %>%
  mutate('T' = 30)

a<-SpatialPoints(sp %>% dplyr::select(X, Y), proj4string=CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs +type=crs"))
d<-as.dist(pointDistance(a))

cosa<-vegan::metaMDS(d)

sp<-cbind(sp, cosa$points)

ggplot(sp, aes(x = MDS1, y = MDS2))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_point(aes(color = Host, shape = UM))+
  coord_fixed()+
  theme_classic()

z<-dist_groups(d = d, g = paste(sp$Host, substr(sp$UM, 1, 1)))

z %>%
  filter(substr(Group1, 1, 6) == substr(Group2, 1, 6)) %>%
  group_by(Label) %>%
  summarise(d = mean(Distance)) %>%
  View()

z %>%
  filter(substr(Group1, 1, 6) == substr(Group2, 1, 6)) %>%
  mutate(grupo = unlist(lapply(strsplit(as.character(Label), ' '), function(x) x[1]))) %>%
  group_by(grupo) %>%
  summarise(d = mean(Distance)) %>%
  View()
