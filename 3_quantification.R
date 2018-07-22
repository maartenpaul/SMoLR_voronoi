folder <- "F:/180704 GG1 IR 5Gy 2h/GG1 IR 2h"
#folder <- "D:/Maarten/GG1"
condition <- "GG1 IR 2h"

source("0 Functions.R")

load(file=file.path(folder,"localizations_clustered.Rdata"))

features <- ldply(localizations,function(x) smolr_point_features(x,skeleton = F))
write.table(features,file.path(folder,"features.txt"))

#measure length by minimal maximal path length, takes a long time
# length <- llply(localizations[1],  function(x) {
#   daply(x, .variables = c("Channel","cluster"),.progress = "text", function(x){
#     delaunay <- make_delaunay(x$X,x$Y)
#     edge_coord <- edge_id_to_coord(edges = delaunay$edges,points = delaunay$points)
#     edge_coord$weight<- edge_coord$dist
#     edge_coord$dist <- NULL
#     routes_igraph <- graph_from_data_frame(d = edge_coord[,-c(3,4,5,6)], directed = FALSE)
#     maxdist <- max(distances(routes_igraph))
#     return(maxdist)
#   })
#   
# })

#####Analyse overlap: take RAD51 clusters and quantifiy overlap with BRCA2 localizations
library(fields)

radius <- 100
ratio <- list()
for (i in 1:length(localizations)){
R51 <- localizations[[i]][localizations[[i]]$Channel=="RAD51",]
R51 <- subset(R51,cluster>0)
B2 <- localizations[[i]][localizations[[i]]$Channel=="BRCA2",]
#B2 <- subset(B2,cluster>0) #select if only analyse B2 localization in clusters, or all BRCA2
ratio[[i]] <- daply(R51,.variables = "cluster",function(x){
  out <- fields.rdist.near(x[,2:3],B2[,2:3],delta = radius,max.points = 1E7)
  if(length(out$ind)>2){
    N_near<- table(out$ind[,1])
    N <- length(N_near[N_near>1])
    ratio <- N/nrow(x)
    return(ratio)
  } else {
    return(0)
  }
})
}
  

#loop over clusters and measure overlap per cell

######Take BRCA2 clusters and get overlap with RAD51

######Number of B2 clusters in vincinity of R51
radius <- 500 #nm

cluster_distances <- ddply(features,.variables = ".id" ,function(x){
  R51_x <- subset(x,Channel=="RAD51")$meanX
  R51_y <- subset(x,Channel=="RAD51")$meanY
  B2_x <- subset(x,Channel=="BRCA2")$meanX
  B2_y <- subset(x,Channel=="BRCA2")$meanY
  
  out <- fields.rdist.near(cbind(R51_x,R51_y),cbind(B2_x,B2_y),delta = radius,max.points = 1E7)
  if(length(out$ind)>2){
    tab_out <- table(out$ind[,1])
    return(data.frame(tab_out))
  }
  else return(data.frame(0))
})

######Minimal distance

#N localizations, RAD51, BRCA2

