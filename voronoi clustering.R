#library(deldir)
library(plyr)
library(dplyr)
library(igraph)
library(reticulate)

#configuration of python, only run once
#py_config()
# py_install("scipy")
# 
#py_install("matplotlib")
#py_install("qt")

TS_to_smolr <- function(x){
  names(x)[3] <- "Y"
  names(x)[2] <- "X"
  names(x)[9] <- "Precision"
  return(x)
}

#conversion for imported thunderstorm files
#convert to point pattern
#localizations_ppp <- SMOLR_TO_PPP(x,marks=x$Channel)

#source_python("qt.py")
make_voronoi <- function(x,y){
  
  #fac: threshold of length between points
  scipy <- import("scipy.spatial")
  #py_install("matplotlib")
  matlabplot <- import("matplotlib")
  
  os <- import("os")

  coord <- data.frame(".id"=seq(1:length(x)),"X"=x,"Y"=y)
  out <- scipy$Voronoi(as.matrix(coord[,2:3]))
  return(out)
}


get_voronoi_area <- function(x){
  source_python('get_area.py')
  
  points <- x$points
  coord <- data.frame(points)
  names(coord) <- c("X","Y")
  regions <- x$regions 
  vertices <- x$vertices
  coord$tile.id <- (x$point_region)+1
  
  #sort regions per point
  regions <- regions[coord$tile.id]
  
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
  
  
  coord$area[coord$tile.border==0] <- ldply(get_area(region_coord))
  return(coord)
}
  

plot_voronoi <- function(x){
  library(ggplot2)
  library(plotly)
  points <- x$points
  points <- data.frame("X"=points[,1],"Y"=points[,2])
  ridge_vertices <- ldply(x$ridge_vertices)+1
  vertices <- x$vertices
  ridge_vertices[ridge_vertices==0] <- NA
  ridge_vertices <- na.omit(ridge_vertices)
  lines <- as.data.frame(cbind(vertices[ridge_vertices[,1],],vertices[ridge_vertices[,2],]))
  names(lines) <- c("x1","y1","x2","y2")
  p <- ggplot(data=points, aes(x=X,y=Y)) +xlim(0,35000)+ylim(0,35000)+
    #Plot the voronoi lines
    geom_segment(
      aes(x = x1, y = y1, xend = x2, yend = y2),
      size = 0.5,
      data = lines,
      linetype = 1,
      color= "#FFB958") + 
    #Plot the points
    geom_point(
      fill=rgb(70,130,180,255,maxColorValue=255),
      pch=21,
      size = 1,
      color="#333333") 
  ggplotly(p)
}


make_delaunay <- function(x,y){

#fac: threshold of length between points
scipy <- import("scipy.spatial")
coord <- as.matrix(cbind(x,y))
out <- scipy$Delaunay(coord)

simplices <- out$simplices
edges <- rbind(simplices[,1:2],simplices[,2:3],simplices[,c(1,3)])
edges <- as.data.frame(rbind(edges,edges[,c(2,1)]))
edges <- as.data.frame(t(apply(edges,1,sort)))
edges <- unique(edges)
edges <- edges+1

points <- out$points

return(list("edges" = edges,"points"=points))
}

edge_id_to_coord <- function(points,edges){
  edges_coord <- data.frame(cbind(edges,points[edges[,1],],points[edges[,2],]))
  names(edges_coord) <- c("ID1","ID2","X1","Y1","X2","Y2")
  edges_coord$dist <- sqrt((edges_coord$X1-edges_coord$X2)^2+(edges_coord$Y1-edges_coord$Y2)^2) 
  return(edges_coord)
}

