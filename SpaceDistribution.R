
# положение звезды в декартовых координатах в галактической системе координат
# Z - вверх, 
# X - в строну центра галактики 
# Y - в наравлении движения солнца (в направлении вращения Галактики)
CalcGalXYZ <- function(data)
{
  return (data <- mutate (data, z = (1/gPx)*sin(gb), x = (1/gPx)*cos(gb)*cos(gl), y = (1/gPx)*cos(gb)*sin(gl)))
}

DrawGalaxyPlane <- function(data, plane = "XY", title = "Star distribution", save = NULL, dscale = 5)
{
  
  #data <- mutate (data, z = (1/gPx)*sin(b), x = (1/gPx)*cos(b)*cos(l), y = (1/gPx)*cos(b)*sin(l))
  #data <- CalcGalXYZ(data)
  
  min_x <- -dscale
  max_x <- dscale
  min_y <- (-1*dscale)
  max_y <- dscale
  
  if(nrow(data)<1000)
  {
    alpha_ <- 1
    size_ <- 1
  } else if(nrow(data)<10000)
  {
    alpha_ <- 0.5
    size_ <- 1
  } else if (nrow(data)<100000)
  {
    alpha_ <- 0.25
    size_ <- 0.5
  } else if (nrow(data)<1000000)
  {
    alpha_ <- 0.1
    size_ <- 0.1
  } else if (nrow(data)<10000000)
  {
    alpha_ <- 0.01
    size_ <- 0.01
   } else if (nrow(data)<100000000)
   {
     alpha_ <- 0.01
     size_ <- 0.01
  } else
   {
     alpha_ <- 0.01
     size_ <- 0.001
   }
  
  g <- ggplot(data)
  
  
  if ( (plane == "XY") | (plane == "YX"))
  {
    
    g <- g + geom_point(aes(x = x, y = y),  alpha = alpha_, na.rm = TRUE, size = size_, shape = ".")
    g <- g + xlab("X") + ylab("Y")
  } else if ((plane == "XZ") | (plane == "ZX"))
  {
    g <- g + geom_point(aes(x = x, y = z),  alpha = alpha_, na.rm = TRUE, size = size_, shape = ".")
    g <- g + xlab("X") + ylab("Z")
  }  else
  {
    g <- g + geom_point(aes(x = y, y = z),  alpha = alpha_, na.rm = TRUE, size = size_, shape = ".")
    g <- g + xlab("Y") + ylab("Z")
  }
  
  #min_x <- (max(min(hrdata[,1]),-9)%/%1)*1 - 1
  #max_x <- (min(max(hrdata[,1]), 9)%/%1)*1 + 1
  #min_y <- (max(min(hrdata[,2]), -9)%/%1)*1 - 1
  #max_y <- (min(max(hrdata[,2]), 9)%/%1)*1 + 1
  
  
  g <- g +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    scale_x_continuous(breaks=seq(min_x,max_x,by=1), minor_breaks=seq(min_x,max_x,by=0.5), limits = c(min_x,max_x)) +
    ggtitle(title)

  
  if(!is.null(save))
  {
    ggsave(paste0(save, plane, ".png"), width = 10, height = 10)
    #ggsave(paste0(save, plane, ".eps"), width = 10, height = 10)
  }
  
  return(g)
}


#ggplot(ga, aes(x = x, y = y)) + stat_bin2d(binwidth = 0.01, aes(fill = ..count..), alpha = 0.8) + theme_bw() + scale_x_continuous(breaks = seq(-5, 5, 1), limits = c(-5, 5)) + scale_y_continuous(breaks = seq(-5, 5, 1), limits = c(-5, 5)) + scale_fill_gradient(low = "white", high = "black", tran = "log")
#ggplot(ga, aes(x = y, y = z)) + stat_bin2d(binwidth = 0.01, aes(fill = ..count..), alpha = 0.8) + theme_bw() + scale_x_continuous(breaks = seq(-5, 5, 1), limits = c(-5, 5)) + scale_y_continuous(breaks = seq(-5, 5, 1), limits = c(-5, 5)) + scale_fill_gradient(low = "white", high = "black", tran = "log")
#ggplot(gaia, aes(x = x, y = y)) + stat_bin2d(binwidth = 0.01, aes(fill = ..count..), alpha = 0.8) + theme_bw() + scale_x_continuous(breaks = seq(-10, 10, 1), limits = c(-10, 10)) + scale_y_continuous(breaks = seq(-10, 10, 1), limits = c(-10, 10)) + scale_fill_gradient(low = "white", high = "black", tran = "log")
