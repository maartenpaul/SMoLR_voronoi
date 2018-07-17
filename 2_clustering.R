folder <- "C:/Thunderstorm test/Nik/GG1"
#folder <- "D:/Maarten/GG1"
condition <- "GG1 IR"

source("0 Functions.R")

load(file = file.path(folder,"localizations_ROI.Rdata"))

localizations <- llply(localizations,function(x){
  x$area <- 0
  x$cluster <- 0
  return(x)
})

voronoi <- llply(localizations,function(x){
  dlply(x,.variables = "Channel",function(x){
    make_voronoi(x$X,x$Y)
  })
})


####analysis for BRCA2
voronoi_area <- llply(voronoi, function(x){
  
    get_voronoi_area(x$BRCA2)
  
})


area_threshold <- llply(voronoi_area,function(x){
     med <- median(x$area[x$tile.border==0])
    return(med)
})

for (i in 1:length(voronoi_area)){
  voronoi_area[[i]] <- subset(voronoi_area[[i]],area<=area_threshold[i]&area!=0)
}


clusters <- llply(voronoi_area,function(x){
   get_clusters(x,fac = 50,conn = 3)               

})

for (i in 1:length(localizations)){
  B2 <- localizations[[i]][localizations[[i]]$Channel=="BRCA2",]
  B2[clusters[[i]]$loc_id,]$cluster <- clusters[[i]]$cluster_id
  localizations[[i]][localizations[[i]]$Channel=="BRCA2",] <- B2
}



####analysis for RAD51
voronoi_area <- llply(voronoi, function(x){
  
  get_voronoi_area(x$RAD51)
  
})


area_threshold <- llply(voronoi_area,function(x){
  med <- median(x$area[x$tile.border==0])
  return(med)
})

for (i in 1:length(voronoi_area)){
  voronoi_area[[i]] <- subset(voronoi_area[[i]],area<=area_threshold[i]&area!=0)
}


clusters <- llply(voronoi_area,function(x){
  get_clusters(x,fac = 20,conn = 3,min_size = 1000)               
  
})
for (i in 1:length(localizations)){
  R51 <- localizations[[i]][localizations[[i]]$Channel=="RAD51",]
  R51[clusters[[i]]$loc_id,]$cluster <- clusters[[i]]$cluster_id
  localizations[[i]][localizations[[i]]$Channel=="RAD51",] <- R51
}
save(localizations,file=file.path(folder,"localizations_clustered.Rdata"))
