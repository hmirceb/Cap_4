library(tidyverse)
library(circlize) # para los roscos
library(scales) # para la escala de colores

load("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")
load("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 4 - Plant population size interactions/Resultados/PERDIVER_Comm_Data.RData")

# Creamos la matriz de adjacencia. El nombre de las columnas de igual porque la funcion
# considera la primera columna como el origen, la segunda el destino y la tercera el 
# ancho para la flecha. En este caso son From (El Orden del interactuante), UM (La poblacion del Host),
# y value (la abundancia del interactuante en esa poblacion). Host y to estan ahi por si acaso
adjacencyData_all<-PERDIVER %>%
  dplyr::select(c(Host, UM, Species, order)) %>%
  group_by(Host, UM, order) %>%
  summarize(value = n()) %>%
  ungroup() %>%
  distinct() %>%
  rename(from = order) %>%
  mutate(to = paste(Host, UM)) %>% # Esto en realidad no hace falta
  mutate(Host = recode(Host, 'rammyc'='Ramonda myconi', 
                       'cypcal'='Cypripedium calceolus',
                       'genlut' = 'Gentiana lutea',
                       'galniv' = 'Galanthus nivalis',
                       'borpyr' = 'Borderea pyrenaica',
                       'kracer' = 'Krascheninnikovia ceratoides',
                       'pinlon' = 'Pinguicula longifolia')) %>%
  dplyr::select(c(from, UM, value, Host, to)) %>%
  filter(complete.cases(.))

# Make the circular plot
paleta<-hue_pal()(4) # Paleta de 4 colores fija para las poblaciones
grid.col<-c(S = paleta[1], 'L1' = paleta[2], 'L2' = paleta[3], 'L3' = paleta[4])

# Elegimos colores aleatorios para los Ordenes
col_mat = rand_color(length(unique(adjacencyData_all$from)), transparency = 0.5)
names(col_mat)<-unique(adjacencyData_all$from) # Hay que ponerle los nombres

circos.par(start.degree = 0, # para que los plots salgan alineados
           canvas.xlim = c(-2, 2), # aumentamos el tamanyo alrededor del plot para que las etiquetas no se corten
           canvas.ylim = c(-2, 2))
chordDiagram(adjacencyData_all, annotationTrack = "grid", # de momento solo ponemos los cuadrados para las etiquetas
             grid.col = c(grid.col, col_mat), # colores de los cuadrados
             col = col_mat, # colores de los enlaces
             small.gap = 3.5, # espacio entre grupos de origen, para que no solapen las etiquetas
             preAllocateTracks = 1) # hueco para poner las etiquetas
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim") # coordenadas x e y para poner cada etiqueta
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index") # los "sectores" son los cuadrados para etiquetas
  circos.text(mean(xlim), # pone la etiqueta en el medio del sector
              ylim[1] + .1, # la alejamos un poco
              sector.name, #etiqueta
              facing = "clockwise", # alineacion en sentido horario
              niceFacing = TRUE, # alinea para que se lea bien
              adj = c(0, 0.5)) # ajustes de la etiqueta, mejor no tocar
}, bg.border = NA) # esto anyade un borde al cuadrado (en este caso es NA porque no lo queremos)
circos.clear() # ESTA LINEA ES IMPORTANTE SI VAMOS A HACER VARIOS PLOTS O SI VAMOS A CAMBIAR ALGUNA DE LAS OPCIONES
              # PARA HACER PRUEBAS. SI NO SE MANTIENEN Y SALEN COSAS RARAS
title('Titulo', line = -4, adj = 0.47) # Podemos poner un titulo con los comandos de R base

############################
# Con varios roscos a la vez
############################

pdf(file = 'Roscos.pdf')
mus<-unique(adjacencyData_all$Host)
paleta<-hue_pal()(4)
grid.col<-c(S = paleta[1], 'L1' = paleta[2], 'L2' = paleta[3], 'L3' = paleta[4])
col_mat = rand_color(length(unique(adjacencyData_all$from)), transparency = 0.5)
names(col_mat)<-unique(adjacencyData_all$from)

for (ii in 1:length(mus)) {
  adjacencyData<-adjacencyData_all %>%
    filter(Host == mus[[ii]])
  
  circos.par(start.degree = 0,
             canvas.xlim = c(-2, 2), 
             canvas.ylim = c(-2, 2))
  chordDiagram(adjacencyData, annotationTrack = "grid", 
               grid.col = c(grid.col, col_mat),
               col = col_mat,
               small.gap = 3.5,
               preAllocateTracks = 1)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA)
  circos.clear()
  title(mus[ii], line = -4, adj = 0.47)
}
dev.off()


