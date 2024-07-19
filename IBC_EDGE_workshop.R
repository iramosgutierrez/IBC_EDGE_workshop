setwd("~")

# 0. Loading our functions####
library(ape)
source("workshop_functions.R")




# 1. Expanding phylogenies - the randtip package####
library(devtools)
install_github("iramosgutierrez/randtip")
library(randtip)

ericatree <- read.tree("files/Erica_incoplete.tre") # A sampled phylogeny
erica.table <- read.csv("files/ericatable.csv") #list of species & info

erica.info <- build_info(erica.table$taxa, ericatree, mode="list",prior.info = plants.info)


cape2include <- erica.table$taxa[erica.table$region == "CAPE" & erica.table$incl==0]
capeincluded <- erica.table$taxa[erica.table$region == "CAPE" & erica.table$incl==1]


erica.info <- edit_info(info = erica.info, taxa = cape2include, column = "taxon1", edit =  "Erica_pauciovulata")
erica.info <- edit_info(info = erica.info, taxa = cape2include, column = "taxon2", edit =  "Erica_dregei")

erica.input <- info2input(erica.info, ericatree, F)
ericatree.complete <- rand_tip(input = erica.input, ericatree, verbose = F)

plot.phylo(ericatree.complete, cex=0.6,tip.color = put_tip_col(ericatree.complete,ericatree, ) )

# 2. EDGE score calculations####
library(phylobase)
 

  GE2vals <- GE.2.calc(pext.vals)
 
  erica.pext <- data.frame("species" = erica.table$taxa ,
                           "RL"= erica.table$RL.cat,
                           "pext" = NA)
  for(sp in erica.pext$species){
    progressbar(which(erica.pext$species == sp), length(erica.pext$species))
    cat.i <- erica.pext$RL[erica.pext$species==sp]
    if(cat.i =="NE"){pext.i <- runif(1,0.0001, 0.9999)}else
      if(cat.i =="EX"){cat.i =="CR"}else{
    pext.i <- sample(GE2vals$pext[GE2vals$RL.cat==cat.i], size = 1)}
    erica.pext$pext[erica.pext$species== sp] <- pext.i
  }
  erica.pext2 <- erica.pext
  erica.pext2$RL<-NULL
  erica.edge <- EDGE2_mod(tree = ericatree.complete, pext = erica.pext2)

  edgevalues <- erica.edge[[1]]
  epdloss.trees <- erica.edge[[2]]
  EricaPD <- erica.edge[[3]]


# 3. Spatial applications####
library(raster)
library(sf)
library(reshape)
library(ggplot2)
library(RColorBrewer)

# 3.0 preparing the data####
erica <- read.csv(file = 'files/Erica_occ_cleaned.csv')


r25<-raster(xmn=-35, xmx=60, ymn=-37, ymx=72,
            crs="+proj=longlat +datum=WGS84 +no_defs ",
            resolution=c(0.25,0.25), vals = NA)

erica2 <- erica
coordinates(erica2) <- ~x+y
proj4string(erica2)<- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
r <- rasterize(erica2, r25, field="taxon", fun=function(x,...) length(unique(na.omit(x))))
plot(r, col = 'red')

# Extract values into data frame
erica.sf <- st_as_sf(erica2)
erica.grid <- c()
for(i in 1:length(unique(erica.sf$taxon))){ # takes a couple minutes to run
  progressbar(i, length(unique(erica.sf$taxon)))
  cells <- as.data.frame(raster::extract(r, erica.sf[which(erica.sf$taxon==unique(erica.sf$taxon)[i]),], cellnumbers=T))
  y<-as.data.frame(raster::extract(r, erica.sf[which(erica.sf$taxon==unique(erica.sf$taxon)[i]),], cellnumbers =T))
  y$taxon <- unique(erica.sf$taxon)[i]
  y <- cbind(y, xyFromCell(r,y[,1]))
  colnames(y)[4:5] <- c('X','Y')
  erica.grid <- rbind(erica.grid, y)
}

erica.grid$X_Y <- paste(erica.grid$X, erica.grid$Y, sep = '_')
erica.grid$layer<-NULL
erica.grid<- unique(erica.grid)

# Save
#save(erica.grid, file = './Erica/Gridded_erica.RData')
# load(file = './Erica/Gridded_erica.RData')


# 3.1 Mapping Taxon richness ####
# Grid to taxon richness raster
rich.rast <- erica.grid
rich.rast <- cast(melt(rich.rast[,c('taxon', 'X', 'Y')], id = c('X', 'Y')), fun.aggregate = length)
rich.rast <- as.data.frame(rich.rast)
coordinates(rich.rast) <- ~X+Y
gridded(rich.rast) <- T
rich.rast<-raster(rich.rast)
plot(rich.rast, col = 'red')

# from SPDF to DF, retaining ID and trename columns
df <- as.data.frame(erica2, xy=T)
ID.df<-cbind(df, xyFromCell(r,cellFromXY(r,df[c("x","y")][c(1:nrow(df)),])))
ID.df$xy <- NULL
ID.df$X <- NULL
colnames(ID.df) <- c('taxon', 'ID',  'x', 'y','trenames', 'X', 'Y')
ID.df$X_Y <- paste(ID.df$X, ID.df$Y, sep = '_')
ID.df$x <- NULL
ID.df$y <- NULL
# write.csv(ID.df, file = './Erica/Erica-ID-grid.csv')
# ID.df<-read.csv(file = './Erica/Erica-ID-grid.csv')


# Load worldmap 
library(rnaturalearth)
library(rnaturalearthhires)
library(dplyr)
worldmap <- ne_countries(scale = 'large', returnclass = 'sf')
worldmap <- worldmap %>% filter(admin!= 'Antarctica')
cropmap <- worldmap$geometry
sf_use_s2(FALSE)
cropmap<-st_crop(cropmap, st_bbox(rich.rast))
# colour scheme
cols <- brewer.pal(10, 'RdYlBu')
pal <- colorRampPalette(cols)


# Richness raster
erica.grid2 <- unique(ID.df[,c('taxon','X','Y','X_Y')])
rich.rast <- erica.grid2
rich.rast <- cast(melt(rich.rast[,c('taxon', 'X', 'Y')], id = c('X', 'Y')), fun.aggregate = length)
rich.rast <- as.data.frame(rich.rast)
coordinates(rich.rast) <- ~X+Y
gridded(rich.rast) <- T
rich.rast<-raster(rich.rast)
plot(rich.rast)

# Raster to polygon 
ep <- rasterToPolygons(rich.rast) # ep erica polygon
ep <- st_as_sf(ep)
ep$log <- log(ep$taxon)

# add log sequence
pos<-seq(min(ep$log), max(ep$log), length.out = 10)
lab.2 <- exp(pos)
lab.3 <- seq(1, max(ep$taxon), length.out= 11) 

customplot(cropmap, ep2, "log(taxon richness)", col="log")


# cropped version
st_crs(ep) <- st_crs(cropmap)
ep2 <- st_intersection(ep, cropmap)
vals <- as.character(cut(ep2$log, breaks = 10, labels = rev(pal(10))))
pos<-seq(min(ep2$log), max(ep2$log), length.out = 10)
lab.2 <- exp(pos)
lab.3 <- seq(1,max(ep$taxon), length.out= 11)




# Now for South Africa zoom 
cropmap.sa<-st_crop(cropmap, st_bbox(c(xmin = 15, ymin = -35, 
                                       xmax = 34, ymax = -21)))
ep2.sa <- st_intersection(ep2, cropmap.sa)
customplot(cropmap.sa, ep2.sa, "Taxon richness", col="log")

# 3.2 EDGE species richness####
edge <- read.csv(file = 'files/Erica_EDGE_spp.csv')


# Erica EDGE species grid
edge.grid <- erica.grid2[which(erica.grid2$taxon%in%edge$Species),] #135 grid cells 

# Grid to taxon richness raster
rs <- edge.grid
rs <- cast(melt(rs[,c('taxon', 'X', 'Y')], id = c('X', 'Y')), fun.aggregate = length)
rs <- as.data.frame(rs)
coordinates(rs) <- ~X+Y
gridded(rs) <- T
rs<-raster(rs)
plot(rs, col = 'red')

# Raster to polygon 
ep3 <- rasterToPolygons(rs)
ep3 <- st_as_sf(ep3)
cols <- brewer.pal(10, 'RdYlBu')
pal <- colorRampPalette(cols)
st_crs(ep3) <- st_crs(cropmap)
ep3 <- st_intersection(ep3, cropmap)
vals <- as.character(cut(ep3$taxon, breaks = 7, labels = rev(pal(7))))

customplot(cropmap, ep3, name="EDGE species richness", col="taxon")

# 3.3 Phylogenetic diversity####
trees <- read.tree(file = 'files/Erica_trees2.tre') # trees after imputation 


PD.grid <- data.frame('X_Y'=unique(erica.grid2$X_Y),'PD.med' = NA)

for(i in 1:length(unique(erica.grid2$X_Y))){
  progressbar(i, length(unique(erica.grid2$X_Y)))
  xy <- unique(erica.grid2$X_Y)[i]
  spp <- unique(erica.grid2$taxon[which(erica.grid2$X_Y==xy)])
  spp <- unique(erica$trnames[which(erica$taxon%in%spp)])
  pd <- c()
  for(j in 1:100){
    x <- sum(keep.tip(trees[[j]], spp)$edge.length)
    pd <- c(pd, x)
  }
  PD.grid$PD.med[which(PD.grid$X_Y==xy)] <- median(pd)
  
}
#save(PD.grid, file = 'data/Erica.pd.RData')
#load('data/Erica.pd.RData')

pdx <- unique(merge(PD.grid, erica.grid2[,c('X', 'Y', 'X_Y')], by = 'X_Y'))
pdx$X_Y <- NULL
coordinates(pdx) <- ~X+Y
gridded(pdx) <- T
pdx<-raster(pdx)

ep4 <- rasterToPolygons(pdx)
ep4 <- st_as_sf(ep4)
cols <- brewer.pal(10, 'RdYlBu')
pal <- colorRampPalette(cols)
st_crs(ep4) <- st_crs(cropmap)
ep4 <- st_intersection(ep4, cropmap)
customplot(cropmap, ep4, "Phylogenetic Diversity", col = "PD.med")

# 3.4 Threatened Evolutionary history####
epdtree <- read.tree(file = 'files/erica_tree_epd2.tre') # epdloss trees after imputation 


ePD.grid <- data.frame('X_Y'=unique(erica.grid2$X_Y),'ePDloss.med' = NA)
for(i in 1:length(unique(erica.grid2$X_Y))){
  progressbar(i, length(unique(erica.grid2$X_Y)))
  xy <- unique(erica.grid2$X_Y)[i]
  spp <- unique(erica.grid2$taxon[which(erica.grid2$X_Y==xy)])
  spp <- unique(erica$trnames[which(erica$taxon%in%spp)])
  pd <- c()
  for(j in 1:100){
    x <- sum(keep.tip(epdtree[[j]], spp)$edge.length)
    pd <- c(pd, x)
  }
  ePD.grid$ePDloss.med[which(ePD.grid$X_Y==xy)] <- median(pd)
}
#save(ePD.grid, file = 'data/Erica.epdloss.RData')
#load(file = 'data/Erica.epdloss.RData')

# Grid to taxon richness raster
pdx2 <- unique(merge(ePD.grid, erica.grid2[,c('X', 'Y', 'X_Y')], by = 'X_Y'))
pdx2$X_Y <- NULL
coordinates(pdx2) <- ~X+Y
gridded(pdx2) <- T
pdx2<-raster(pdx2)


# Raster to polygon 
ep5 <- rasterToPolygons(pdx2)
ep5 <- st_as_sf(ep5)

st_crs(ep5) <- st_crs(cropmap)
ep5 <- st_intersection(ep5, cropmap)
customplot(cropmap, ep5, "Expected PD loss", col= "ePDloss.med" )


cropmap.ib<-st_crop(cropmap, st_bbox(c(xmin = -10, ymin = 30, 
                                       xmax = 10, ymax =  50)))
ep5.ib <- st_intersection(ep5, cropmap.ib)

customplot(cropmap.ib, ep5.ib, "Expected PD loss", col= "ePDloss.med" )