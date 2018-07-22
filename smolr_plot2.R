SMOLR_PLOT2 <- function(x,color=NULL,size=NULL,skip=FALSE){
  library(ggplot2)
  library(plotly)
  if(class(x)!="data.frame") {stop("SMOLR_PLOT2 requires a data frame as input")}
  
  ind_x <- grep("^x$",names(x),ignore.case=T)
  ind_y <- grep("^y$",names(x),ignore.case=T)
  ind_ch <- grep("^ch",names(x),ignore.case=T)
  ind_prec <- grep("^prec",names(x),ignore.case=T)
  
  if(is.null(color)){
    color <- names(x)[ind_ch]
  }
    
    if(is.null(size)){
      size <- names(x)[ind_prec]
    }
  
 # if(length(c(ind_x,ind_y,ind_ch,ind_prec))!=4){stop("Not all parameters (x,y,channel,precision) are present once in the header")}
  p <- ggplot(data=x,aes_string(x=names(x)[ind_x],y=names(x)[ind_y],color=color,size=size))+geom_point(size=0.5)+xlim(0,35000)+ylim(0,35000)
  ggplotly(p)
}

