regions <- x$regions 
vertices <- x$vertices
coord$tile.id <- (x$point_region)+1
coord$.id <- seq(1:nrow(coord))

#sort regions per point
regions <- regions[coord$tile.id]
names(regions) <- coord$.id 
#Idenfity points on the border
coord$tile.border <- laply(regions,function(x){
  if (any(x==-1)){
    return(1)
  } else {
    return(0)
  }
  
})

region_coord <- llply(regions[coord$tile.border==0],function(x){
  x <- x+1
  return(vertices[x,])
})
names(region_coord) <- coord$.id[coord$tile.border==0]

region_coord2 <- region_coord
names(region_coord2) <- NULL
 areas <- get_area(region_coord2)
coord$area[coord$tile.border==0] <- unlist(areas)
med_area <- median(coord$area[coord$tile.border==0])

df <- ldply(region_coord,function(x){
  data.frame("X"=x[,1],"Y"=x[,2])
})
datapoly <- merge(coord, df, by = c(".id"))

datapoly <- datapoly[datapoly$area<(2*med_area),]

p <- ggplot(datapoly, aes(x = X.y, y = Y.y)) +
  geom_polygon(aes(fill = area, group = .id))+xlim(0,35000)+ylim(0,35000)
p
ggplotly(p)

#####Analyse overlap: take RAD51 clusters and quantifiy overlap with BRCA2 localizations

#loop over clusters and measure overlap per cell

######Take BRCA2 clusters and get overlap with RAD51

######Number of clusters in vincinity

#take 500 nm

######Minimal distance

#N localizations, RAD51, BRCA2
