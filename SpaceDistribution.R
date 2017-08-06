
CalcGalXYZ <- function(data)
{
  data <- mutate (data, z = (1/gPx)*sin(gb), x = (1/gPx)*cos(gb)*cos(gl), y = (1/gPx)*cos(gb)*sin(gl))
}

DrawGalaxyPlane <- function(data, plane = "XZ", title = "Star distribution", save = NULL, dscale = 5)
{
  
  #data <- mutate (data, z = (1/gPx)*sin(b), x = (1/gPx)*cos(b)*cos(l), y = (1/gPx)*cos(b)*sin(l))
  data <- CalcGalXYZ(data)
  
  g <- ggplot()
  
  if ( (plane == "XY") | (plane == "YX"))
  {
    hrdata <- data.frame(cbind(data$x, data$y))
    g <- g + xlab("X") + ylab("Y")
  } else if ((plane == "XZ") | (plane == "ZX"))
  {
    hrdata <- data.frame(cbind(data$x, data$z))
    g <- g + xlab("X") + ylab("Z")
  }
  else
  {
    hrdata <- data.frame(cbind(data$y, data$z))
    g <- g + xlab("Y") + ylab("Z")
  }
  
  #min_x <- (max(min(hrdata[,1]),-9)%/%1)*1 - 1
  #max_x <- (min(max(hrdata[,1]), 9)%/%1)*1 + 1
  #min_y <- (max(min(hrdata[,2]), -9)%/%1)*1 - 1
  #max_y <- (min(max(hrdata[,2]), 9)%/%1)*1 + 1
  
  min_x <- -dscale
  max_x <- dscale
  min_y <- -dscale
  max_y <- dscale
  
  if(nrow(hrdata)<1000)
  {
    alpha_ <- 1
    size_ <- 1
  } else if(nrow(hrdata)<10000)
  {
    alpha_ <- 0.5
    size_ <- 1
  } else if (nrow(hrdata)<100000)
  {
    alpha_ <- 0.25
    size_ <- 0.5
  } else if (nrow(hrdata)<1000000)
  {
    alpha_ <- 0.1
    size_ <- 0.1
  } else
  {
    alpha_ <- 0.05
    size_ <- 0.1
  }
  
  g <- g +
    geom_point(data=hrdata, aes(x = hrdata[,1], y = hrdata[,2]),  alpha = alpha_, na.rm = TRUE, size = size_, shape = ".") + 
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    scale_x_continuous(breaks=seq(min_x,max_x,by=1), minor_breaks=seq(min_x,max_x,by=0.5), limits = c(min_x,max_x)) +
    ggtitle(title)
  
  
  if(!is.null(save))
  {
    ggsave(paste0(save, plane, ".png"), width = 10, height = 10)
    ggsave(paste0(save, plane, ".eps"), width = 10, height = 10)
  }
  
  return(g)
}