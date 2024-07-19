# Erica mapping

# Working directory.
setwd("E:/UNI/4. DOCTORADO/99. OTROS/2024_06_14 EDGE Workshop/materials_workshop/")

# Library
library(rWCVP)
library(sf)
library(reshape)
library(oce)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(data.table)
library(reshape2)
library(raster)
library(sp)
library(rnaturalearth)
library(ape)
library(phytools)

# Data
erica <- read.csv(file = './Erica/Erica_occ_cleaned.csv')
edge <- read.csv(file = './Erica/Erica_EDGE_spp.csv')
og <- read.tree(file = './Erica_100tree_nameschanged.tre') # original trees
trees <- read.tree(file = './Erica/Erica_trees2.tre') # trees after imputation 
epdtree <- read.tree(file = './Erica/erica_tree_epd2.tre') # epdloss trees after imputation 

# Making a raster
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
  cells <- as.data.frame(raster::extract(r, erica.sf[which(erica.sf$taxon==unique(erica.sf$taxon)[i]),], cellnumbers=T))
  y<-as.data.frame(raster::extract(r, erica.sf[which(erica.sf$taxon==unique(erica.sf$taxon)[i]),], cellnumbers =T))
  y$taxon <- unique(erica.sf$taxon)[i]
  y <- cbind(y, xyFromCell(r,y[,1]))
  colnames(y)[4:5] <- c('X','Y')
  erica.grid <- rbind(erica.grid, y)
  print(i)
}

erica.grid$X_Y <- paste(erica.grid$X, erica.grid$Y, sep = '_')
erica.grid$layer<-NULL
erica.grid<- unique(erica.grid)

# Save
 #save(erica.grid, file = './Erica/Gridded_erica.RData')
# load(file = './Erica/Gridded_erica.RData')

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
ID.df<-cbind(df, xyFromCell(r,cellFromXY(r,df[c(3,4)][c(1:nrow(df)),])))
ID.df$xy <- NULL
colnames(ID.df) <- c('taxon', 'ID', 'trenames', 'x', 'y', 'X', 'Y')
ID.df$X_Y <- paste(ID.df$X, ID.df$Y, sep = '_')
ID.df$x <- NULL
ID.df$y <- NULL
# write.csv(ID.df, file = './Erica/Erica-ID-grid.csv')
# ID.df<-read.csv(file = './Erica/Erica-ID-grid.csv')



#----------------------------------------------------------

# Mapping Taxon richness. 

# Checking
erica.grid2 <- unique(ID.df[,c('taxon','X','Y','X_Y')])
rich.rast <- erica.grid2
rich.rast <- cast(melt(rich.rast[,c('taxon', 'X', 'Y')], id = c('X', 'Y')), fun.aggregate = length)
rich.rast <- as.data.frame(rich.rast)
coordinates(rich.rast) <- ~X+Y
gridded(rich.rast) <- T
rich.rast<-raster(rich.rast)
plot(rich.rast, col = 'red')

# Raster to polygon 
ep <- rasterToPolygons(rich.rast) # ep erica polygon
ep <- st_as_sf(ep)
ep$log <- log(ep$taxon)

# add log sequence
pos<-seq(min(ep$log), max(ep$log), length.out = 10)
lab.2 <- exp(pos)
lab.3 <- seq(1, max(ep$taxon), length.out= 11) 


# Load worldmap 
library(rnaturalearthdata)
library(rnaturalearthhires)
worldmap <- ne_countries(scale = 'large', returnclass = 'sf')
worldmap <- worldmap %>% filter(admin!= 'Antarctica')
cropmap <- worldmap$geometry
sf_use_s2(FALSE)
cropmap<-st_crop(cropmap, st_bbox(rich.rast))

# colour scheme
cols <- brewer.pal(10, 'RdYlBu')
pal <- colorRampPalette(cols)

# cropped version
st_crs(ep) <- st_crs(cropmap)
ep2 <- st_intersection(ep, cropmap)
vals <- as.character(cut(ep2$log, breaks = 10, labels = rev(pal(10))))
pos<-seq(min(ep2$log), max(ep2$log), length.out = 10)
lab.2 <- exp(pos)
lab.3 <- seq(1,max(ep$taxon), length.out= 11)


ggplot() + 
  geom_sf(fill = "antiquewhite1") +
  geom_sf(data = cropmap, fill = 'grey') +
  geom_sf(data = ep2, aes(fill = log, color = log), lwd = .01) + 
  scale_fill_gradientn(colours=rev(pal(10)), name = 'Taxon richness',
                       labels = as.character(c(round(lab.2))),
                       breaks = pos) +
  scale_color_gradientn(colours = rev(pal(10)), guide = 'none') +
  theme(panel.grid.major = element_line(color = gray(0.3), 
                                        linetype = "dashed", 
                                        size = 0.1), 
        panel.background = element_rect(fill = "antiquewhite1"),
        legend.position = 'right',
        legend.key.height = unit(1.3, "cm"))


# Now for South Africa zoom 
cropmap.sa<-st_crop(cropmap, st_bbox(c(xmin = 15, ymin = -35, 
                                       xmax = 34, ymax = -21)))
ep2.sa <- st_intersection(ep2, cropmap.sa)


ggplot() + 
  geom_sf(fill = "antiquewhite1") +
  geom_sf(data = cropmap.sa, fill = 'grey') +
  geom_sf(data = ep2.sa, aes(fill = log, color = log), lwd = .01) + 
  scale_fill_gradientn(colours=rev(pal(10)), name = 'Taxon richness',
                       labels = as.character(c(round(lab.2))),
                       breaks = pos) +
  scale_color_gradientn(colours = rev(pal(10)), guide = 'none') +
  theme(panel.grid.major = element_line(color = gray(0.3), 
                                        linetype = "dashed", 
                                        size = 0.1), 
        panel.background = element_rect(fill = "antiquewhite1"),
        legend.position = 'right',
        legend.key.height = unit(1.3, "cm"))







#---------------------------------------------------------------------------

# Part 2. EDGE richness
# Relies on the EDGE species computation being complete

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
#zoom(rs, col = 'red', drawExtent())

# Raster to polygon 
ep3 <- rasterToPolygons(rs)
ep3 <- st_as_sf(ep3)
cols <- brewer.pal(10, 'RdYlBu')
pal <- colorRampPalette(cols)
st_crs(ep3) <- st_crs(cropmap)
ep3 <- st_intersection(ep3, cropmap)
vals <- as.character(cut(ep3$taxon, breaks = 7, labels = rev(pal(7))))

# Load worldmap 
worldmap <- ne_countries(scale = 'large', returnclass = 'sf')
worldmap <- worldmap %>% filter(admin!= 'Antarctica')
cropmap <- worldmap$geometry
sf_use_s2(FALSE)
cropmap<-st_crop(cropmap, st_bbox(rich.rast))

ggplot() + 
  geom_sf(fill = "antiquewhite1") +
  geom_sf(data = cropmap, fill = 'grey') +
  geom_sf(data = ep3, aes(fill = taxon, color = taxon), lwd = .01) + 
  scale_fill_gradientn(colours=rev(pal(10)), name = 'EDGE species richness richness',
                       labels = seq(1, max(ep3$taxon), 1), 
                       breaks = seq(1, max(ep3$taxon), 1)) +
  scale_color_gradientn(colours = rev(pal(10)), guide = 'none') +
  theme(panel.grid.major = element_line(color = gray(0.3), 
                                        linetype = "dashed", 
                                        size = 0.1), 
        panel.background = element_rect(fill = "antiquewhite1"),
        legend.position = 'right',
        legend.key.height = unit(1.3, "cm"))




# -------------------------------------------------------------------------

# 3. PD


# Calculating PD - takes 10 - 15 minutes
# there is likely a quicker way of doing this. 
PD.grid <- data.frame('X_Y'=unique(erica.grid2$X_Y),'PD.med' = NA)
for(i in 1:length(unique(erica.grid2$X_Y))){
  xy <- unique(erica.grid2$X_Y)[i]
  spp <- unique(erica.grid2$taxon[which(erica.grid2$X_Y==xy)])
  spp <- unique(erica$trnames[which(erica$taxon%in%spp)])
  pd <- c()
  for(j in 1:100){
    x <- sum(keep.tip(trees[[j]], spp)$edge.length)
    pd <- c(pd, x)
  }
  PD.grid$PD.med[which(PD.grid$X_Y==xy)] <- median(pd)
  print(i)
}
#save(PD.grid, file = './Erica/Erica.pd.RData')
#load('./Erica/Erica.pd.RData')
# Joining to erica grid.

# Grid to taxon richness raster
pdx <- unique(merge(PD.grid, erica.grid2[,c('X', 'Y', 'X_Y')], by = 'X_Y'))
pdx$X_Y <- NULL
coordinates(pdx) <- ~X+Y
gridded(pdx) <- T
pdx<-raster(pdx)
plot(pdx)

# Raster to polygon 
ep4 <- rasterToPolygons(pdx)
ep4 <- st_as_sf(ep4)
cols <- brewer.pal(10, 'RdYlBu')
pal <- colorRampPalette(cols)
st_crs(ep4) <- st_crs(cropmap)
ep4 <- st_intersection(ep4, cropmap)

# Load worldmap 
worldmap <- ne_countries(scale = 'large', returnclass = 'sf')
worldmap <- worldmap %>% filter(admin!= 'Antarctica')
cropmap <- worldmap$geometry
sf_use_s2(FALSE)
cropmap<-st_crop(cropmap, st_bbox(pdx))

ggplot() + 
  geom_sf(fill = "antiquewhite1") +
  geom_sf(data = cropmap, fill = 'grey') +
  geom_sf(data = ep4, aes(fill = PD.med, color = PD.med), lwd = .01) + 
  scale_fill_gradientn(colours=rev(pal(10)), name = 'EDGE species richness richness',
                       labels = seq(1, maxValue(pdx), round(maxValue(pdx)/10)), 
                       breaks = seq(1, maxValue(pdx), round(maxValue(pdx)/10))) +
  scale_color_gradientn(colours = rev(pal(10)), guide = 'none') +
  theme(panel.grid.major = element_line(color = gray(0.3), 
                                        linetype = "dashed", 
                                        size = 0.1), 
        panel.background = element_rect(fill = "antiquewhite1"),
        legend.position = 'right',
        legend.key.height = unit(1.3, "cm"))


# -------------------------------------------------------------------------

# 4. ePD Loss


# Calculating PD
ePD.grid <- data.frame('X_Y'=unique(erica.grid2$X_Y),'ePDloss.med' = NA)
for(i in 1:length(unique(erica.grid2$X_Y))){
  xy <- unique(erica.grid2$X_Y)[i]
  spp <- unique(erica.grid2$taxon[which(erica.grid2$X_Y==xy)])
  spp <- unique(erica$trnames[which(erica$taxon%in%spp)])
  pd <- c()
  for(j in 1:100){
    x <- sum(keep.tip(epdtree[[j]], spp)$edge.length)
    pd <- c(pd, x)
  }
  ePD.grid$ePDloss.med[which(ePD.grid$X_Y==xy)] <- median(pd)
  print(i)
}
#save(ePD.grid, file = './Erica/Erica.epdloss.RData')
#load(file = './Erica/Erica.epdloss.RData')
# Joining to erica grid.


# Grid to taxon richness raster
pdx2 <- unique(merge(ePD.grid, erica.grid2[,c('X', 'Y', 'X_Y')], by = 'X_Y'))
pdx2$X_Y <- NULL
coordinates(pdx2) <- ~X+Y
gridded(pdx2) <- T
pdx2<-raster(pdx2)
plot(pdx2)
#zoom(pdx2, drawExtent())

# Raster to polygon 
ep5 <- rasterToPolygons(pdx2)
ep5 <- st_as_sf(ep5)
cols <- brewer.pal(10, 'RdYlBu')
pal <- colorRampPalette(cols)
st_crs(ep5) <- st_crs(cropmap)
ep5 <- st_intersection(ep5, cropmap)

# Load worldmap 
worldmap <- ne_countries(scale = 'large', returnclass = 'sf')
worldmap <- worldmap %>% filter(admin!= 'Antarctica')
cropmap <- worldmap$geometry
sf_use_s2(FALSE)
cropmap<-st_crop(cropmap, st_bbox(pdx))


ggplot() + 
  geom_sf(fill = "antiquewhite1") +
  geom_sf(data = cropmap, fill = 'grey') +
  geom_sf(data = ep5, aes(fill = ePDloss.med, color = ePDloss.med), lwd = .01) + 
  scale_fill_gradientn(colours=rev(pal(10)), name = 'EDGE species richness richness',
                       labels = seq(1, maxValue(pdx2), round(maxValue(pdx2)/10)), 
                       breaks = seq(1, maxValue(pdx2), round(maxValue(pdx2)/10))) +
  scale_color_gradientn(colours = rev(pal(10)), guide = 'none') +
  theme(panel.grid.major = element_line(color = gray(0.3), 
                                        linetype = "dashed", 
                                        size = 0.1), 
        panel.background = element_rect(fill = "antiquewhite1"),
        legend.position = 'right',
        legend.key.height = unit(1.3, "cm"))





#----------------------------------------------------------

# Uniting all the grids
PD.grid
ePD.grid
ID.df
edge.grid <- erica.grid2[which(erica.grid2$taxon%in%edge$Species),] #135 grid cells 

tax <- cast(melt(unique(ID.df[,c('taxon', 'X_Y')]), id = c('X_Y')), fun.aggregate = length)
edge.grid <- cast(melt(unique(edge.grid[,c('taxon', 'X_Y')]), id = c('X_Y')), fun.aggregate = length)
colnames(edge.grid)[2] <- 'Edge.richness'
colnames(tax)[2] <- 'Taxon.richness'


allvars <- merge(PD.grid, ePD.grid, by = 'X_Y', all.x = T)
allvars <- merge(allvars, tax, by = 'X_Y', all.x = T)
allvars <- merge(allvars, edge.grid, by = 'X_Y', all.x = T)
allvars[is.na(allvars)] <- 0
allvars <- merge(allvars, unique(ID.df[,c('X_Y', 'X', 'Y')]), by = 'X_Y', all.x = T)
allvars <- allvars[,c(1,6,7,2,3,4,5)]
head(allvars)

# write.csv(allvars, file = './Erica/Erica_grid_all_variables.csv', row.names = F)




