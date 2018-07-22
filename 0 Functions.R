#library(deldir)
library(plyr)
library(dplyr)
library(igraph)
library(reticulate)
library(SMoLR)
library(ggplot2)
library(plotly)

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
  coord$area <- 0
  
 
  #this loop is implemented to avoid to large input for the Python function, which resulted in crashes
  if(length(region_coord)>100000){
    parts <- length(region_coord)%/%100000
    rest <- length(region_coord)%%100000
    areas <- get_area(region_coord[1:100000])
    i <- 1
    if(parts>1){
    for (i in 2:parts){
      areas <- c(areas,get_area(region_coord[(1+(i-1)*100000):(i*100000)]))
    }
    }
    if (rest>0){
      areas <- c(areas,get_area(region_coord[(1+(i)*100000):length(region_coord)]))
    }
  } else { 
    areas <- get_area(region_coord)
    
    
  }
  
  
  coord$area[coord$tile.border==0] <- unlist(areas)
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

get_clusters <- function(coord, fac=50,conn=3,min_size=10){
  delaunay <- make_delaunay(coord$X,coord$Y)
  edge_coord <- edge_id_to_coord(edges = delaunay$edges,points = delaunay$points)
  points <- delaunay$points
  
  cutoff = (edge_coord$dist <= fac)
  edge_out <- edge_coord[cutoff,]
  edge_out$weight<- edge_out$dist
  edge_out$dist <- NULL
  
  routes_igraph <- graph_from_data_frame(d = edge_out[,-c(3,4,5,6)], directed = FALSE)
  #calculate degree of connections
  degree <- degree(routes_igraph)
  
  #removed edges with low connectivity
  routes_igraph <- delete_edges(routes_igraph, which(as.vector(degree)<=conn))
  #vert <- igraph::as_data_frame(routes_igraph,what="vertices")
  #names(vert) <- c("X","Y")
  subgraphs <- decompose.graph(routes_igraph)
  sizes <- laply(subgraphs, function(x) gsize(x))
  subgraphs <- subgraphs[sizes>min_size]
  
  #get the edge coordinates
  # clusters <- llply(subgraphs,function(x) {
  #   unique(apply(as_edgelist(x),2,as.numeric))
  # })
  
  clusters_coord <- llply(subgraphs,function(x){
    z <- as_data_frame(x,what="vertices")
    z <- as.numeric(z$name)
    loc_id <- as.numeric(row.names(coord)[z])
    coord <- coord[z,1:2]
    
    out <- data.frame(loc_id,coord)
    names(out) <- c("loc_id","X","Y")
    return(out)
  })
  names(clusters_coord) <- seq(1,length(clusters_coord))
  
  clusters_coord <- ldply(clusters_coord)
  names(clusters_coord) <- c("cluster_id","loc_id","X","Y")
  clusters_coord$cluster_id <- as.numeric(clusters_coord$cluster_id)
  return(clusters_coord)
}

smolr_point_features <- function(x,skeleton=TRUE){
  #avoid issues with capital letters
  ind_cluster <- grep("^cluster$",names(x),ignore.case=T)
  ind_x <- grep("^x$",names(x),ignore.case=T)
  ind_y <- grep("^y$",names(x),ignore.case=T)
  ind_ch <- grep("^ch",names(x),ignore.case=T)
  ind_prec <- grep("^prec",names(x),ignore.case=T)
  
  x$cluster <- x[,ind_cluster]
  x <- x[x$cluster>0,]
  #go through the localizations and get statistics splitted by channel and cluster
  out <- ddply(x,.variables = c("Channel","cluster"),function(x) {
    
    dx <- x[,ind_x]
    y <- x[,ind_y]
    ch <- x[,ind_ch]
    prec <- x[,ind_prec]
    
    
    D <- princomp(cbind(dx,y))
    coord <- cbind(dx,y)
    angle <- atan2(D$loadings[2,1],D$loadings[1,1])
    
    chull_coord <- chull(coord)
    area <- pracma::polyarea(coord[rev(chull_coord),1], coord[rev(chull_coord),2])
    perimeter <- pracma::poly_length(coord[rev(chull_coord),1], coord[rev(chull_coord),2])
    
    if(skeleton){
      xlim <- c(min(dx)-50,max(dx)+50)
      ylim <-  c(min(y)-50,max(y)+50)   
      kde <- SMOLR_KDE(x,threshold = 0.1,px=5,bandwidth = c(10,10),xlim=xlim,ylim=ylim)$kde_binary[,,1]
      
      skeleton <- thinImage(as.matrix(kde))
      #from pixels to nanometers    
      length_skeleton <- length(skeleton[skeleton==1])*5
      
      
      return(
        data.frame("meanX" = mean(dx),"meanY" = mean(y),
                   "sd" = ((D$sdev[1] + D$sdev[2])/2) * 2.35,"width" = (max(D$scores[,1]) - min(D$scores[,1])),"area"=area,"perimeter"=perimeter,
                   "major_axis"= D$sdev[1]*2.35,"minor_axis"= D$sdev[2]*2.35 ,"ratio" = (D$sdev[1] /D$sdev[2]),"angle" = angle,"N" = nrow(x),"skeleton"=length_skeleton
        )
      )
    } else {
      return(
        data.frame("meanX" = mean(dx),"meanY" = mean(y),
                   "sd" = ((D$sdev[1] + D$sdev[2])/2) * 2.35,"width" = (max(D$scores[,1]) - min(D$scores[,1])),"area"=area,"perimeter"=perimeter,
                   "major_axis"= D$sdev[1]*2.35,"minor_axis"= D$sdev[2]*2.35 ,"ratio" = (D$sdev[1] /D$sdev[2]),"angle" = angle,"N" = nrow(x),"skeleton"=0
        )
      ) 
    }
  })
  
}