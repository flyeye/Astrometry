


HRDiagram <- function(data, photometric = "APASS", title = "Hertzsprung-Russell", save = NULL, 
                      L5lim = TRUE, L3lim = TRUE, BV_lim = c(-1, 3), M_lim = c(10, -10))
{
  if (photometric == "TYCHO")
  {
    data <- data %>% mutate(M = NA, B_V = NA)
    s <- (data$gPx>0) & (!is.na(data$tyc_m))
    data$M[s] <- data$tyc_m[s] + 5 + 5*log10(data$gPx[s]/1000)
    hrdata <- data.frame(cbind( M = data$M[!is.na(data$M)], B_V = data$tyc_bv[!is.na(data$M)]))
  } else if (photometric == "APASS")
  {
    data <- data %>% mutate(M = NA)
    index <-(data$gPx>0)&(!is.na(data$apasm_v))
    data$M[index] <- data$apasm_v[index] + 5 + 5*log10(data$gPx[index]/1000)
    hrdata <- data.frame(cbind( M = data$M[index], B_V = (data$apasm_b[index]-data$apasm_v[index]), LC = data$LClass[!is.na(data$M)]))
  }
  else if (photometric == "TGAS")
  {
    data <- data %>% mutate(M = NA)
    data$M[data$gPx>0] <- data$Gm_mag[data$gPx>0] + 5 + 5*log10(data$gPx[data$gPx>0]/1000)
    hrdata <- data.frame(cbind( M = data$M[!is.na(data$M)], B_V = data$B_V[!is.na(data$M)], LC = data$LClass[!is.na(data$M)]))
  } else 
  {
    hrdata <- data.frame(cbind( M = data$M[!is.na(data$M)], B_V = data$B_V[!is.na(data$M)], LC = data$LClass[!is.na(data$M)]))
  }
  
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
  
  #g <- ggplot() + geom_point(data=hrdata, aes(x = hrdata$B_V, y = hrdata$M), alpha_ = 0.05, na.rm = TRUE, size_ = 0.1) + scale_y_reverse()
  
  g <- ggplot() +
    #scale_colour_manual("L Class",  breaks = colnames(c(1,2,3,4,5)),
    #           values = c("blue", "brown", "red", "yellow", "green")) +
    #scale_color_continuous(high = "red", low = "blue") + 
    geom_point(data=hrdata, aes(x = hrdata$B_V, y = hrdata$M),  alpha = alpha_, na.rm = TRUE, size = size_, shape = ".") +  
    scale_y_reverse(breaks=seq(M_lim[1],M_lim[2],by=-1), minor_breaks=seq(M_lim[1],M_lim[2],by=-0.5), limits = M_lim) +
    scale_x_continuous(breaks=seq(-1,3,by=0.25), minor_breaks=seq(-1,3,by=0.125), limits = c(-1,3)) +
    xlab("B-V") + ylab("M") + ggtitle(title)
  
  if (L5lim)
  {
    #ms_top_limit <- matrix(0, nrow = 5, ncol = 2)
    #ms_top_limit[,1] <- c(-0.1, 0.8, 1.0, 1.5, 1.7)
    #ms_top_limit[,2] <- c(-1.8, 3.4, 6.0, 8.0, 10.0 )
    #ms_bottom_limit <- matrix(0, nrow = 4, ncol = 2)
    #ms_bottom_limit[,1] <- c(-0.1, 0.5, 1.3, 1.35)
    #ms_bottom_limit[,2] <- c(1.8, 5.6, 9.0, 10.0)
    
    #g <- g + geom_line(aes(x = ms_top_limit[,1], y = ms_top_limit[,2])) +
    #        geom_line(aes(x = ms_bottom_limit[,1], y = ms_bottom_limit[,2]));
    
    a <- seq(from = -1, to = 2.0, by = 0.01)
    g <- g + geom_line(aes(x = a, y = max_M(a))) +
      geom_line(aes(x = a, y = min_M(a)));
  }
  
  if (L3lim)
  {
    rg_top_limit <- matrix(0, nrow = 4, ncol = 2)
    rg_top_limit[,1] <- c(0.8, 0.8, 2.5, 2.5)
    rg_top_limit[,2] <- c(2.5, -1.5, -1.5, 2.5)
    rg_bottom_limit <- matrix(0, nrow = 2, ncol = 2)
    rg_bottom_limit[,1] <- c(0.8, 2.5)
    rg_bottom_limit[,2] <- c(2.5, 2.5)
    
    g <- g + geom_line(aes(x = rg_top_limit[,1], y = rg_top_limit[,2])) +
      geom_line(aes(x = rg_bottom_limit[,1], y = rg_bottom_limit[,2]))
  }
  
  
  if(!is.null(save))
  {
    ggsave(paste0(save, "-HR.jpg"), width = 10, height = 10)
    ggsave(paste0(save, "-HR.eps"), width = 10, height = 10)
  }
  
  return(g)
}