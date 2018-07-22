folder <- "F:/180704 GG1 IR 5Gy 2h/GG1 IR 2h"
#folder <- "D:/Maarten/GG1"
condition <- "GG1 IR 2h"

source("0 Functions.R")

remove(localizations)
load(file = file.path(folder,"localizations_ROI.Rdata"))

start <-paste("Clustering started:", format(Sys.time(), "%a %b %d %X %Y"))
print(start)

ptm <- proc.time()

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

#calculate area threshold
area_threshold <- llply(voronoi_area,function(x){
     med <- median(x$area[x$tile.border==0])*2
    return(med)
})

#set threshold
for (i in 1:length(voronoi_area)){
  voronoi_area[[i]] <- subset(voronoi_area[[i]],area<=area_threshold[i]&area!=0)
}

#assign cluster ID to the localizations
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

for (i in 1:length(localizations)){
 write.table(localizations[[i]],file = file.path(folder,paste0("localizations_clustered_",i,".txt")),col.names = T)
}
print(paste("Clustering finished:", format(Sys.time(), "%a %b %d %X %Y")))
x <- proc.time() - ptm
print(paste("It took", as.character(x[3]/60),"minutes for the script to run"))
