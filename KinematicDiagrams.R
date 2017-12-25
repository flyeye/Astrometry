draw_OMParComp <- function(parameter, sol1, sol2, title = "", xat = "", yat = "")
{
  
  min_x <- (min(sol1$X[,parameter])%/%0.2)*0.2 - 1
  max_x <- (max(sol1$X[,parameter])%/%0.2)*0.2 + 1
  stepx <- 1
  min_y <- (min(sol2$X[,parameter])%/%0.2)*0.2 - 1
  max_y <- (max(sol2$X[,parameter])%/%0.2)*0.2 + 1
  stepy <- 1
  
  g<- ggplot()
  
  g <- g +
    #scale_y_continuous(breaks=seq(-18,18,by=3), minor_breaks=seq(-18,18,by=1), limits = c(-17,17)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=stepy), minor_breaks=seq(min_y,max_y,by=stepy/2), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0.25,4,by=0.25), minor_breaks=seq(0.25,4,by=0.125), limits = c(0.25,4))
    scale_x_continuous(breaks=seq(min_x,max_x,by=stepx), minor_breaks=seq(min_x,max_x,by=stepx/2), limits = c(min_x,max_x))
  #scale_x_continuous(breaks=seq(1,4.5,by=0.5), minor_breaks=seq(1,4.5,by=0.25), limits = c(1,4.5))
  
  g <- g + xlab(xat) + ylab(yat) +ggtitle(title)
  
  g <- g + 
    geom_line(aes(x = sol1$X[,parameter], y = sol2$X[,parameter]), size = 1) + 
    geom_errorbar(aes(x = sol1$X[,parameter], ymin = sol2$X[,parameter] - sol1$S_X[,parameter], ymax = sol2$X[,parameter] + sol1$S_X[,parameter])) +
    geom_errorbarh(aes(x = sol1$X[,parameter], xmin = sol1$X[,parameter] - sol2$S_X[,parameter], xmax = sol1$X[,parameter] + sol2$S_X[,parameter], y = sol2$X[,parameter])) + 
    geom_point(aes(x = sol1$X[,parameter], y = sol2$X[,parameter])) #+  
  
  return (g)
}



draw_OM <- function(res, title = "Ogorodnikov-Miln Model")
{
  min_x <- (min(res$Parameters[,4])%/%0.2)*0.2
  max_x <- (max(res$Parameters[,4])%/%0.2)*0.2 + 0.2
  stepx <- 0.2
  #stepx <- 10
  
  #min_y <- (min(res$X[,4:11])%/%1) - 1
  #max_y <- (max(res$X[,4:11])%/%1) + 1
  min_y <- -17
  max_y <- 18
  
  
  g <- ggplot()
  g <- g +
    #scale_y_continuous(breaks=seq(-18,18,by=3), minor_breaks=seq(-18,18,by=1), limits = c(-17,17)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0.25,4,by=0.25), minor_breaks=seq(0.25,4,by=0.125), limits = c(0.25,4))
    scale_x_continuous(breaks=seq(min_x,max_x,by=stepx), minor_breaks=seq(min_x,max_x,by=stepx/2), limits = c(min_x,max_x))
  #scale_x_continuous(breaks=seq(1,4.5,by=0.5), minor_breaks=seq(1,4.5,by=0.25), limits = c(1,4.5))
  
  g <- g + xlab("<px>, kpc") + ylab("O-M, km/s/kpc") +ggtitle(title) +
    scale_colour_manual("Parameters",  breaks = colnames(res$X)[4:11],
                        values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "green", "#0072B2", "#D55E00", "#CC79A7")) +
    #geom_line(aes(x = res$Parameters[,4], y = res$X[,1], colour = colnames(res$X)[1]), size = 1) +
    #geom_line(aes(x = res$Parameters[,4], y = res$X[,2], colour = colnames(res$X)[2]), size = 1) +
    #geom_line(aes(x = res$Parameters[,4], y = res$X[,3], colour = colnames(res$X)[3])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,4], colour = colnames(res$X)[4]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,4] - res$S_X[,4], ymax = res$X[,4] + res$S_X[,4], colour = colnames(res$X)[4])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,5], colour = colnames(res$X)[5]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,5] - res$S_X[,5], ymax = res$X[,5] + res$S_X[,5], colour = colnames(res$X)[5])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,6], colour = colnames(res$X)[6]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,6] - res$S_X[,6], ymax = res$X[,6] + res$S_X[,6], colour = colnames(res$X)[6])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,7], colour = colnames(res$X)[7]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,7] - res$S_X[,7], ymax = res$X[,7] + res$S_X[,7], colour = colnames(res$X)[7])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,8], colour = colnames(res$X)[8]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,8] - res$S_X[,8], ymax = res$X[,8] + res$S_X[,8], colour = colnames(res$X)[8])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,9], colour = colnames(res$X)[9]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,9] - res$S_X[,9], ymax = res$X[,9] + res$S_X[,9], colour = colnames(res$X)[9])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,10], colour = colnames(res$X)[10]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,10] - res$S_X[,10], ymax = res$X[,10] + res$S_X[,10], colour = colnames(res$X)[10])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,11], colour = colnames(res$X)[11]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,11] - res$S_X[,11], ymax = res$X[,11] + res$S_X[,11], colour = colnames(res$X)[11])) +
    geom_point(aes(x = res$Parameters[,4], y = rep(0,nrow(res$Parameters)))) #+  
  #geom_smooth(aes(x = res$Parameters[,4], y = res$X[,6], colour = colnames(res$X)[6])) +
  #geom_smooth(aes(x = res$Parameters[,4], y = res$X[,9], colour = colnames(res$X)[9])) +
  #geom_smooth(aes(x = res$Parameters[,4], y = res$X[,10], colour = colnames(res$X)[10])) 
  return(g)
}


draw_OM_diff <- function(res, title = "Ogorodnikov-Miln Model")
{
  min_x <- (min(res$Parameters[,4])%/%0.2)*0.2
  max_x <- (max(res$Parameters[,4])%/%0.2)*0.2 + 0.2
  
  min_y <- (min(res$X[,4:11])%/%1) - 1
  max_y <- (max(res$X[,4:11])%/%1) + 1
  
  g <- ggplot()
  
  g <- g +
    #scale_y_continuous(breaks=seq(-6,7,by=0.5), minor_breaks=seq(-6,7,by=0.25), limits = c(-6,7)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0,4,by=0.25), minor_breaks=seq(0,4,by=0.125), limits = c(0.2,4))
    #scale_x_continuous(breaks=seq(0,2.5,by=0.25), minor_breaks=seq(0,2.5,by=0.125), limits = c(0.2,2.5))
    scale_x_continuous(breaks=seq(min_x,max_x,by=0.2), minor_breaks=seq(min_x,max_x,by=0.1), limits = c(min_x,max_x))
  
  
  g <- g + xlab("<px>, kpc") + ylab("O-M, km/s/kpc") +ggtitle(title) +
    scale_colour_manual("Parameters",  breaks = colnames(res$X)[4:11],
                        values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "green", "#0072B2", "#D55E00", "#CC79A7")) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,4], colour = colnames(res$X)[4]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,5], colour = colnames(res$X)[5]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,6], colour = colnames(res$X)[6]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,7], colour = colnames(res$X)[7]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,8], colour = colnames(res$X)[8]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,9], colour = colnames(res$X)[9]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,10], colour = colnames(res$X)[10]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,11], colour = colnames(res$X)[11]), size = 1) +
    geom_point(aes(x = res$Parameters[,4], y = rep(0,nrow(res$Parameters))))
  return(g)
}

draw_OM_Solar <- function(res, title = "Solar motion")
{
  min_x <- (min(res$Parameters[,4])%/%0.2)*0.2
  max_x <- (max(res$Parameters[,4])%/%0.2)*0.2 + 0.2
  
  min_y <- (min(res$X[,1:3])%/%1) - 1
  max_y <- (max(res$X[,1:3])%/%1) + 1
  
  g <- ggplot() +
    #scale_y_continuous(breaks=seq(0,20,by=1), minor_breaks=seq(0,20,by=0.5), limits = c(-0.5,20)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0,4,by=0.25), minor_breaks=seq(0,4,by=0.125), limits = c(0.25,4)) +
    #scale_x_continuous(breaks=seq(0,2.5,by=0.25), minor_breaks=seq(0,2.5,by=0.125), limits = c(0.25,2.5)) +
    scale_x_continuous(breaks=seq(min_x,max_x,by=0.2), minor_breaks=seq(min_x,max_x,by=0.1), limits = c(min_x,max_x))+
    xlab("<px>, kpc") + ylab("km/s") +ggtitle(title) +
    scale_colour_manual("Parameters",  breaks = colnames(res$X)[1:3],
                        values = c("green", "blue", "brown")) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,1], colour = colnames(res$X)[1]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,1] - res$S_X[,1], ymax = res$X[,1] + res$S_X[,1], colour = colnames(res$X)[1])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,2], colour = colnames(res$X)[2]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,2] - res$S_X[,2], ymax = res$X[,2] + res$S_X[,2], colour = colnames(res$X)[2])) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,3], colour = colnames(res$X)[3]), size = 1) +
    geom_errorbar(aes(x = res$Parameters[,4], ymin = res$X[,3] - res$S_X[,3], ymax = res$X[,3] + res$S_X[,3], colour = colnames(res$X)[3])) +
    geom_point(aes(x = res$Parameters[,4], y = rep(0,nrow(res$Parameters))))
  return(g)
}


draw_OM_Solar_diff <- function(res, title = "Solar motion")
{
  min_x <- (min(res$Parameters[,4])%/%0.2)*0.2
  max_x <- (max(res$Parameters[,4])%/%0.2)*0.2 + 0.2
  
  min_y <- (min(res$X[,1:3])%/%1) - 1
  max_y <- (max(res$X[,1:3])%/%1) + 1
  g <- ggplot() +
    #scale_y_continuous(breaks=seq(-2,7,by=0.5), minor_breaks=seq(-2,7,by=0.25), limits = c(-2,7)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0,4,by=0.5), minor_breaks=seq(0,4,by=0.1), limits = c(0.25,4)) +
    scale_x_continuous(breaks=seq(min_x,max_x,by=0.2), minor_breaks=seq(min_x,max_x,by=0.1), limits = c(min_x,max_x))+
    xlab("<px>, kpc") + ylab("km/s") +ggtitle(title) +
    scale_colour_manual("Parameters",  breaks = colnames(res$X)[1:3],
                        values = c("green", "blue", "brown")) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,1], colour = colnames(res$X)[1]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,2], colour = colnames(res$X)[2]), size = 1) +
    geom_line(aes(x = res$Parameters[,4], y = res$X[,3], colour = colnames(res$X)[3]), size = 1) +
    geom_point(aes(x = res$Parameters[,4], y = rep(0,nrow(res$Parameters))))
  return(g)
}


draw_Oort <- function(res, title = "Oort`s parameters")
{
  min_x <- (min(res$Parameters[,4])%/%0.2)*0.2
  max_x <- (max(res$Parameters[,4])%/%0.2)*0.2 + 0.2
  
  #min_y <- (min(res$Oort)%/%1) - 1
  #max_y <- (max(res$Oort)%/%1) + 1
  min_y <- -30
  max_y <- 30
  
  
  g <- ggplot() +
    #scale_y_continuous(breaks=seq(0,20,by=1), minor_breaks=seq(0,20,by=0.5), limits = c(-0.5,20)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0,4,by=0.25), minor_breaks=seq(0,4,by=0.125), limits = c(0.25,4)) +
    #scale_x_continuous(breaks=seq(0,2.5,by=0.25), minor_breaks=seq(0,2.5,by=0.125), limits = c(0.25,2.5)) +
    scale_x_continuous(breaks=seq(min_x,max_x,by=0.2), minor_breaks=seq(min_x,max_x,by=0.1), limits = c(min_x,max_x))+
    xlab("<r>, kpc") + ylab("km/s/kpc") +ggtitle(title) +
    scale_colour_manual("Parameters",  breaks = colnames(res$Oort)[1:6],
                        values = c("green", "blue", "brown", "black", "gray", "orange"))
  
  for (i in 1:ncol(res$Oort))
  {
    g <- g + geom_line(aes_string(x = res$Parameters[,4], y = res$Oort[,i], colour = shQuote(colnames(res$Oort)[i])), size = 1) +
      geom_errorbar(aes_string(x = res$Parameters[,4], ymin = res$Oort[,i] - res$s_Oort[,i], ymax = res$Oort[,i] + res$s_Oort[,i], colour = shQuote(colnames(res$Oort)[i])))
  }
  # geom_line(aes(x = res$Parameters[,4], y = res$Oort[,1], colour = colnames(res$Oort)[1]), size = 1) +
  # geom_errorbar(aes(x = res$Parameters[,4], ymin = res$Oort[,1] - res$s_Oort[,1], ymax = res$Oort[,1] + res$s_Oort[,1], colour = colnames(res$Oort)[1])) +
  # geom_line(aes(x = res$Parameters[,4], y = res$Oort[,2], colour = colnames(res$Oort)[2]), size = 1) +
  # geom_errorbar(aes(x = res$Parameters[,4], ymin = res$Oort[,2] - res$s_Oort[,2], ymax = res$Oort[,2] + res$s_Oort[,2], colour = colnames(res$Oort)[2])) +
  # geom_line(aes(x = res$Parameters[,4], y = res$Oort[,3], colour = colnames(res$Oort)[3]), size = 1) +
  # geom_errorbar(aes(x = res$Parameters[,4], ymin = res$Oort[,3] - res$s_Oort[,3], ymax = res$Oort[,3] + res$s_Oort[,3], colour = colnames(res$Oort)[3])) +
  # geom_line(aes(x = res$Parameters[,4], y = res$Oort[,4], colour = colnames(res$Oort)[4]), size = 1) +
  # geom_errorbar(aes(x = res$Parameters[,4], ymin = res$Oort[,4] - res$s_Oort[,4], ymax = res$Oort[,4] + res$s_Oort[,4], colour = colnames(res$Oort)[4])) +
  # geom_point(aes(x = res$Parameters[,4], y = rep(0,nrow(res$Parameters))))
  return(g)
}

draw_solution_Oort <- function(solution, title = "Oort`s parameters")
{
  min_x <- 0
  max_x <- 4
  step_x <- 0.5
  
  #min_y <- (min(res$Oort)%/%1) - 1
  #max_y <- (max(res$Oort)%/%1) + 1
  min_y <- -17
  max_y <- 18
  
  g <- ggplot() +
    #scale_y_continuous(breaks=seq(0,20,by=1), minor_breaks=seq(0,20,by=0.5), limits = c(-0.5,20)) +
    scale_y_continuous(breaks=seq(min_y,max_y,by=1), minor_breaks=seq(min_y,max_y,by=0.5), limits = c(min_y,max_y)) +
    #scale_x_continuous(breaks=seq(0,4,by=0.25), minor_breaks=seq(0,4,by=0.125), limits = c(0.25,4)) +
    #scale_x_continuous(breaks=seq(0,2.5,by=0.25), minor_breaks=seq(0,2.5,by=0.125), limits = c(0.25,2.5)) +
    scale_x_continuous(breaks=seq(min_x,max_x,by=step_x), minor_breaks=seq(min_x,max_x,by=step_x/2), limits = c(0.2,3.5))+
    xlab("<r>, kpc") + ylab("km/s/kpc") +ggtitle(title) +
    #scale_colour_manual("Parameters",  breaks = colnames(solution$MS_All$Oort)[1:4],
    #                   values = c("green", "blue", "brown", "black")) +
    scale_colour_manual("Parameters",  breaks = c("A Main Sequence", "A Red Giants Disk", "A Red Giants Galo", 
                                                  "B Main Sequence", "B Red Giants Disk", "B Red Giants Galo", 
                                                  "C Main Sequence", "C Red Giants Disk", "C Red Giants Galo", 
                                                  "K Main Sequence", "K Red Giants Disk", "K Red Giants Galo"),
                        values = c("green1", "green3", "green4", "blue1", "blue3", "blue4","brown1", "brown3", "brown4", "gray0","gray20", "gray40")) +
    scale_fill_manual("Parameters",  breaks = c("A Main Sequence", "A Red Giants Disk", "A Red Giants Galo", 
                                                "B Main Sequence", "B Red Giants Disk", "B Red Giants Galo", 
                                                "C Main Sequence", "C Red Giants Disk", "C Red Giants Galo", 
                                                "K Main Sequence", "K Red Giants Disk", "K Red Giants Galo"),
                      values = c("green1", "green3", "green4", "blue1", "blue3", "blue4","brown1", "brown3", "brown4", "gray0","gray20", "gray40")) +    
    #geom_line(aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,1], colour = colnames(solution$MS_All$Oort)[1]), size = 1) +
    #geom_errorbar(aes(x = solution$MS_All$Parameters[,4], ymin = solution$MS_All$Oort[,1] - solution$MS_All$s_Oort[,1], ymax = solution$MS_All$Oort[,1] + solution$MS_All$s_Oort[,1], colour = colnames(solution$MS_All$Oort)[1])) +
    geom_line(  aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,1], colour = "A Main Sequence"), size = 1) +
    geom_point( aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,1], fill = "A Main Sequence"), shape = 21) + 
    geom_errorbar(aes(x = solution$MS_All$Parameters[,4], ymin = solution$MS_All$Oort[,1] - solution$MS_All$s_Oort[,1], ymax = solution$MS_All$Oort[,1] + solution$MS_All$s_Oort[,1], colour = "A Main Sequence")) +
    
    geom_line( aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,2], colour = "B Main Sequence"), size = 1) +
    geom_point(aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,2], fill = "B Main Sequence"), shape = 21) + 
    geom_errorbar(aes(x = solution$MS_All$Parameters[,4], ymin = solution$MS_All$Oort[,2] - solution$MS_All$s_Oort[,2], ymax = solution$MS_All$Oort[,2] + solution$MS_All$s_Oort[,2], colour = "B Main Sequence")) +
    
    
    geom_line(aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,3], colour = "C Main Sequence"), size = 1) +
    geom_errorbar(aes(x = solution$MS_All$Parameters[,4], ymin = solution$MS_All$Oort[,3] - solution$MS_All$s_Oort[,3], ymax = solution$MS_All$Oort[,3] + solution$MS_All$s_Oort[,3], colour = "C Main Sequence")) +
    geom_point(aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,3], fill = "C Main Sequence"), shape = 21) + 
    
    geom_line(aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,4], colour = "K Main Sequence"), size = 1) +
    geom_errorbar(aes(x = solution$MS_All$Parameters[,4], ymin = solution$MS_All$Oort[,4] - solution$MS_All$s_Oort[,4], ymax = solution$MS_All$Oort[,4] + solution$MS_All$s_Oort[,4], colour = "K Main Sequence")) +
    geom_point(aes(x = solution$MS_All$Parameters[,4], y = solution$MS_All$Oort[,4], fill = "K Main Sequence"), shape = 21) + 
    
    geom_point(aes(x = solution$MS_All$Parameters[,4], y = rep(0,nrow(solution$MS_All$Parameters))))
  
  
  g <- g + 
    #geom_line(aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,1], colour = colnames(solution$RG_Disk$Oort)[1]), size = 1) +
    #geom_errorbar(aes(x = solution$RG_Disk$Parameters[,4], ymin = solution$RG_Disk$Oort[,1] - solution$RG_Disk$s_Oort[,1], ymax = solution$RG_Disk$Oort[,1] + solution$RG_Disk$s_Oort[,1], colour = colnames(solution$RG_Disk$Oort)[1])) +
    geom_line(aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,1], colour = "A Red Giants Disk"), size = 1) +
    geom_point( aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,1], fill = "A Red Giants Disk"), shape = 22) + 
    geom_errorbar(aes(x = solution$RG_Disk$Parameters[,4], ymin = solution$RG_Disk$Oort[,1] - solution$RG_Disk$s_Oort[,1], ymax = solution$RG_Disk$Oort[,1] + solution$RG_Disk$s_Oort[,1], colour = "A Red Giants Disk")) +
    
    geom_line(aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,2], colour = "B Red Giants Disk"), size = 1) +
    geom_point( aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,2], fill = "B Red Giants Disk"), shape = 22) + 
    geom_errorbar(aes(x = solution$RG_Disk$Parameters[,4], ymin = solution$RG_Disk$Oort[,2] - solution$RG_Disk$s_Oort[,2], ymax = solution$RG_Disk$Oort[,2] + solution$RG_Disk$s_Oort[,2], colour = "B Red Giants Disk")) +
    
    geom_line(aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,3], colour = "C Red Giants Disk"), size = 1) +
    geom_point( aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,3], fill = "C Red Giants Disk"), shape = 22) + 
    geom_errorbar(aes(x = solution$RG_Disk$Parameters[,4], ymin = solution$RG_Disk$Oort[,3] - solution$RG_Disk$s_Oort[,3], ymax = solution$RG_Disk$Oort[,3] + solution$RG_Disk$s_Oort[,3], colour = "C Red Giants Disk")) +
    
    geom_line(aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,4], colour = "K Red Giants Disk"), size = 1) +
    geom_point( aes(x = solution$RG_Disk$Parameters[,4], y = solution$RG_Disk$Oort[,4], fill = "K Red Giants Disk"), shape = 22) + 
    geom_errorbar(aes(x = solution$RG_Disk$Parameters[,4], ymin = solution$RG_Disk$Oort[,4] - solution$RG_Disk$s_Oort[,4], ymax = solution$RG_Disk$Oort[,4] + solution$RG_Disk$s_Oort[,4], colour = "K Red Giants Disk")) +
    
    geom_point(aes(x = solution$RG_Disk$Parameters[,4], y = rep(0,nrow(solution$RG_Disk$Parameters))))
  
  g <- g + 
    #geom_line(aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,1], colour = colnames(solution$RG_Galo$Oort)[1]), size = 1) +
    #geom_errorbar(aes(x = solution$RG_Galo$Parameters[,4], ymin = solution$RG_Galo$Oort[,1] - solution$RG_Galo$s_Oort[,1], ymax = solution$RG_Galo$Oort[,1] + solution$RG_Galo$s_Oort[,1], colour = colnames(solution$RG_Galo$Oort)[1])) +
    geom_line(  aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,1], colour = "A Red Giants Galo"), size = 1) +
    geom_point( aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,1], fill = "A Red Giants Galo"), shape = 23) + 
    geom_errorbar(aes(x = solution$RG_Galo$Parameters[,4], ymin = solution$RG_Galo$Oort[,1] - solution$RG_Galo$s_Oort[,1], ymax = solution$RG_Galo$Oort[,1] + solution$RG_Galo$s_Oort[,1], colour = "A Red Giants Galo")) +
    
    geom_line(aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,2], colour = "B Red Giants Galo"), size = 1) +
    geom_point( aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,2], fill = "B Red Giants Galo"), shape = 23) + 
    geom_errorbar(aes(x = solution$RG_Galo$Parameters[,4], ymin = solution$RG_Galo$Oort[,2] - solution$RG_Galo$s_Oort[,2], ymax = solution$RG_Galo$Oort[,2] + solution$RG_Galo$s_Oort[,2], colour = "B Red Giants Galo")) +
    
    geom_line(aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,3], colour = "C Red Giants Galo"), size = 1) +
    geom_point( aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,3], fill = "C Red Giants Galo"), shape = 23) + 
    geom_errorbar(aes(x = solution$RG_Galo$Parameters[,4], ymin = solution$RG_Galo$Oort[,3] - solution$RG_Galo$s_Oort[,3], ymax = solution$RG_Galo$Oort[,3] + solution$RG_Galo$s_Oort[,3], colour = "C Red Giants Galo")) +
    
    geom_line(aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,4], colour = "K Red Giants Galo"), size = 1) +
    geom_point( aes(x = solution$RG_Galo$Parameters[,4], y = solution$RG_Galo$Oort[,4], fill = "K Red Giants Galo"), shape = 23) + 
    geom_errorbar(aes(x = solution$RG_Galo$Parameters[,4], ymin = solution$RG_Galo$Oort[,4] - solution$RG_Galo$s_Oort[,4], ymax = solution$RG_Galo$Oort[,4] + solution$RG_Galo$s_Oort[,4], colour = "K Red Giants Galo")) +
    
    geom_point(aes(x = solution$RG_Galo$Parameters[,4], y = rep(0,nrow(solution$RG_Galo$Parameters))))
  
  # g <- g + geom_line(aes(x = solution$SG_Disk$Parameters[,4], y = solution$SG_Disk$Oort[,1], colour = colnames(solution$SG_Disk$Oort)[1]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Disk$Parameters[,4], ymin = solution$SG_Disk$Oort[,1] - solution$SG_Disk$s_Oort[,1], ymax = solution$SG_Disk$Oort[,1] + solution$SG_Disk$s_Oort[,1], colour = colnames(solution$SG_Disk$Oort)[1])) +
  #   
  #   geom_line(aes(x = solution$SG_Disk$Parameters[,4], y = solution$SG_Disk$Oort[,2], colour = colnames(solution$SG_Disk$Oort)[2]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Disk$Parameters[,4], ymin = solution$SG_Disk$Oort[,2] - solution$SG_Disk$s_Oort[,2], ymax = solution$SG_Disk$Oort[,2] + solution$SG_Disk$s_Oort[,2], colour = colnames(solution$SG_Disk$Oort)[2])) +
  #   
  #   geom_line(aes(x = solution$SG_Disk$Parameters[,4], y = solution$SG_Disk$Oort[,3], colour = colnames(solution$SG_Disk$Oort)[3]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Disk$Parameters[,4], ymin = solution$SG_Disk$Oort[,3] - solution$SG_Disk$s_Oort[,3], ymax = solution$SG_Disk$Oort[,3] + solution$SG_Disk$s_Oort[,3], colour = colnames(solution$SG_Disk$Oort)[3])) +
  #   
  #   geom_line(aes(x = solution$SG_Disk$Parameters[,4], y = solution$SG_Disk$Oort[,4], colour = colnames(solution$SG_Disk$Oort)[4]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Disk$Parameters[,4], ymin = solution$SG_Disk$Oort[,4] - solution$SG_Disk$s_Oort[,4], ymax = solution$SG_Disk$Oort[,4] + solution$SG_Disk$s_Oort[,4], colour = colnames(solution$SG_Disk$Oort)[4])) +
  #   
  #   geom_point(aes(x = solution$SG_Disk$Parameters[,4], y = rep(0,nrow(solution$SG_Disk$Parameters))))
  # 
  # g <- g + geom_line(aes(x = solution$SG_Galo$Parameters[,4], y = solution$SG_Galo$Oort[,1], colour = colnames(solution$SG_Galo$Oort)[1]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Galo$Parameters[,4], ymin = solution$SG_Galo$Oort[,1] - solution$SG_Galo$s_Oort[,1], ymax = solution$SG_Galo$Oort[,1] + solution$SG_Galo$s_Oort[,1], colour = colnames(solution$SG_Galo$Oort)[1])) +
  #   
  #   geom_line(aes(x = solution$SG_Galo$Parameters[,4], y = solution$SG_Galo$Oort[,2], colour = colnames(solution$SG_Galo$Oort)[2]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Galo$Parameters[,4], ymin = solution$SG_Galo$Oort[,2] - solution$SG_Galo$s_Oort[,2], ymax = solution$SG_Galo$Oort[,2] + solution$SG_Galo$s_Oort[,2], colour = colnames(solution$SG_Galo$Oort)[2])) +
  #   
  #   geom_line(aes(x = solution$SG_Galo$Parameters[,4], y = solution$SG_Galo$Oort[,3], colour = colnames(solution$SG_Galo$Oort)[3]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Galo$Parameters[,4], ymin = solution$SG_Galo$Oort[,3] - solution$SG_Galo$s_Oort[,3], ymax = solution$SG_Galo$Oort[,3] + solution$SG_Galo$s_Oort[,3], colour = colnames(solution$SG_Galo$Oort)[3])) +
  #   
  #   geom_line(aes(x = solution$SG_Galo$Parameters[,4], y = solution$SG_Galo$Oort[,4], colour = colnames(solution$SG_Galo$Oort)[4]), size = 1) +
  #   geom_errorbar(aes(x = solution$SG_Galo$Parameters[,4], ymin = solution$SG_Galo$Oort[,4] - solution$SG_Galo$s_Oort[,4], ymax = solution$SG_Galo$Oort[,4] + solution$SG_Galo$s_Oort[,4], colour = colnames(solution$SG_Galo$Oort)[4])) +
  #   
  #   geom_point(aes(x = solution$SG_Galo$Parameters[,4], y = rep(0,nrow(solution$SG_Galo$Parameters))))
  # 
  return(g)
}


#----------------------------------------------------


draw_OMParameter <- function(solution, 
                             parameter = 1,
                             title = "", 
                             x_lim = c(0, 4, 0.5), 
                             y_lim = c(5, 40, 5), 
                             clr = c("blue", "green4", "brown", "black", "red", "orange"),
                             x_par = 4,
                             x_title = "<r>, kpc", 
                             y_title = "km/s/kpc", 
                             is_legend = TRUE)
{
  
  names <- vector("character", 0)
  for (i in 1:length(solution))
  {
    names[i] <- solution[[i]]$Name
  }
  
  if (title == "")
  {
    title <- colnames(solution[[i]]$X)[parameter]
  }
  
  g <- ggplot() +
    #scale_y_continuous(breaks=seq(y_lim[1],y_lim[2],by=y_lim[3]), minor_breaks=seq(y_lim[1],y_lim[2],by=y_lim[3]/2)) +
    #scale_x_continuous(breaks=seq(x_lim[1],x_lim[2],by=x_lim[3]), minor_breaks=seq(x_lim[1],x_lim[2],by=x_lim[3]/2)) +
    scale_y_continuous(breaks=seq(y_lim[1],y_lim[2],by=y_lim[3]), minor_breaks=seq(y_lim[1],y_lim[2],by=y_lim[3]/2), limits = c(y_lim[1],y_lim[2])) +
    scale_x_continuous(breaks=seq(x_lim[1],x_lim[2],by=x_lim[3]), minor_breaks=seq(x_lim[1],x_lim[2],by=x_lim[3]/2), limits = c(x_lim[1],x_lim[2])) +
    xlab(x_title) + ylab(y_title) +
    scale_colour_manual("Parameters",  breaks = names, values = clr) +
    scale_fill_manual("Parameters",  breaks = names, values = clr) +    
    scale_shape_manual("Parameters", breaks = names, values = c(21, 22, 23, 24, 25, 26)) +
    scale_linetype_manual("Parameters", breaks = names, values = c(1, 2, 3, 4, 5, 6)) + 
    coord_cartesian(ylim = c(y_lim[1],y_lim[2]), xlim =c(x_lim[1],x_lim[2])) 
  
  for (i in 1:length(solution))
  {
    
    g <- g + geom_line( aes_string(x = solution[[i]]$Parameters[,x_par], y = solution[[i]]$X[,parameter], 
                                   colour = shQuote(solution[[i]]$Name), 
                                   linetype = shQuote(solution[[i]]$Name)
    ), 
    size = 1)
    
    g <- g + geom_point(aes_string(x = solution[[i]]$Parameters[,x_par], y = solution[[i]]$X[,parameter],
                                   fill = shQuote(solution[[i]]$Name),
                                   shape = shQuote(solution[[i]]$Name)
    ))
    
    g <- g + geom_errorbar(aes_string(x = solution[[i]]$Parameters[,x_par],
                                      ymin = solution[[i]]$X[,parameter] - solution[[i]]$S_X[,parameter],
                                      ymax = solution[[i]]$X[,parameter] + solution[[i]]$S_X[,parameter],
                                      colour = shQuote(solution[[i]]$Name)),
                           width = 0.05)
  }
  
  g <- g + theme_bw() + theme(axis.text=element_text(size=rel(1.0)))
    
  if (is_legend == FALSE)
  {
    g <- g + guides(colour = "none", shape = "none", fill = "none", linetype = "none")
    ## theme(axis.text=element_text(size=18,face="bold"))      
  }
  else {
    g <- g + ggtitle(title)
    
  }
  
  return(g)
}

draw_OortParameter <- function(solution, 
                               parameter = 1,
                               title = "Oort`s parameter A", 
                               x_lim = c(0, 3.5, 0.5), 
                               y_lim = c(8, 18, 1), 
                               clr = c("blue", "green4", "brown", "black", "red", "orange"),
                               x_par = 4,  # 4 - фактическое среднее расстояние выборки
                               x_title = "<r>, kpc", 
                               y_title = "km/s/kpc", 
                               is_legend = TRUE)  
{
  
  names <- vector("character", 0)
  for (i in 1:length(solution))
  {
    names[i] <- solution[[i]]$Name
  }
  
  g <- ggplot() +
    #scale_y_continuous(breaks=seq(y_lim[1],y_lim[2],by=y_lim[3]), minor_breaks=seq(y_lim[1],y_lim[2],by=y_lim[3]/2)) +
    #scale_x_continuous(breaks=seq(x_lim[1],x_lim[2],by=x_lim[3]), minor_breaks=seq(x_lim[1],x_lim[2],by=x_lim[3]/2)) +
    scale_y_continuous(breaks=seq(y_lim[1],y_lim[2],by=y_lim[3]), minor_breaks=seq(y_lim[1],y_lim[2],by=y_lim[3]/2), limits = c(y_lim[1],y_lim[2])) +
    scale_x_continuous(breaks=seq(x_lim[1],x_lim[2],by=x_lim[3]), minor_breaks=seq(x_lim[1],x_lim[2],by=x_lim[3]/2), limits = c(x_lim[1],x_lim[2])) +
    xlab(x_title) + ylab(y_title) +
    scale_colour_manual("Parameters",  breaks = names, values = clr) +
    scale_fill_manual("Parameters",  breaks = names, values = clr) +    
    scale_shape_manual("Parameters", breaks = names, values = c(21, 22, 23, 24, 25, 26)) +
    scale_linetype_manual("Parameters", breaks = names, values = c(1, 2, 3, 4, 5, 6)) + 
    coord_cartesian(ylim = c(y_lim[1],y_lim[2]), xlim =c(x_lim[1],x_lim[2])) 
  
  for (i in 1:length(solution))
  {
    
    g <- g + geom_line( aes_string(x = solution[[i]]$Parameters[,x_par], y = solution[[i]]$Oort[,parameter], 
                                   colour = shQuote(solution[[i]]$Name), 
                                   linetype = shQuote(solution[[i]]$Name)
    ), 
    size = 1)
    
    g <- g + geom_point(aes_string(x = solution[[i]]$Parameters[,x_par], y = solution[[i]]$Oort[,parameter],
                                   fill = shQuote(solution[[i]]$Name),
                                   shape = shQuote(solution[[i]]$Name)
    ))
    
    g <- g + geom_errorbar(aes_string(x = solution[[i]]$Parameters[,x_par],
                                      ymin = solution[[i]]$Oort[,parameter] - solution[[i]]$s_Oort[,parameter],
                                      ymax = solution[[i]]$Oort[,parameter] + solution[[i]]$s_Oort[,parameter],
                                      colour = shQuote(solution[[i]]$Name)), 
                            width = 0.05)
  }
  
  g <- g + theme_bw() + theme(axis.text=element_text(size=rel(1.0)))
  
  if (is_legend == FALSE)
  {
    g <- g + guides(colour = "none", shape = "none", fill = "none", linetype = "none")

  }
  else {
    g <- g + ggtitle(title)
  
  }
  
  return(g)
}



draw_Physical <- function(solution,                        
                          parameter = 1,
                          title = "Linear galactic velocity at Solar distance", 
                          x_lim = c(0, 4, 0.5), y_lim = c(185, 245, 10), 
                          clr = c("blue", "green4", "brown", "black", "red"), 
                          x_par = 4,  # 4 - фактическое среднее расстояние выборки
                          x_title = "<r>, kpc", 
                          y_title = "km/s", 
                          is_legend = TRUE)
{
  
  names <- vector("character", 0)
  for (i in 1:length(solution))
  {
    names[i] <- solution[[i]]$Name
  }
  
  g <- ggplot() +
    scale_y_continuous(breaks=seq(y_lim[1],y_lim[2],by=y_lim[3]), minor_breaks=seq(y_lim[1],y_lim[2],by=y_lim[3]/2), limits = c(y_lim[1],y_lim[2])) +
    scale_x_continuous(breaks=seq(x_lim[1],x_lim[2],by=x_lim[3]), minor_breaks=seq(x_lim[1],x_lim[2],by=x_lim[3]/2), limits = c(x_lim[1],x_lim[2])) +
    xlab(x_title) + ylab(y_title) +
    scale_colour_manual("Parameters",  breaks = names, values = clr) +
    scale_fill_manual("Parameters",  breaks = names, values = clr) +    
    scale_shape_manual("Parameters", breaks = names, values = c(21, 22, 23, 24, 25)) +
    scale_linetype_manual("Parameters", breaks = names, values = c(1, 2, 3, 4, 5)) 
  
  
  for (i in 1:length(solution))
  {
    g <- g + geom_line( aes_string(x = solution[[i]]$Parameters[,x_par], y = solution[[i]]$Physical[,parameter], 
                                   colour = shQuote(solution[[i]]$Name), 
                                   linetype = shQuote(solution[[i]]$Name)), 
                        size = 1)
    
    g <- g + geom_point(aes_string(x = solution[[i]]$Parameters[,x_par], y = solution[[i]]$Physical[,parameter],
                                   fill = shQuote(solution[[i]]$Name),
                                   shape = shQuote(solution[[i]]$Name)))
    
    g <- g + geom_errorbar(aes_string(x = solution[[i]]$Parameters[,x_par],
                                      ymin = solution[[i]]$Physical[,parameter] - solution[[i]]$s_Physical[,parameter],
                                      ymax = solution[[i]]$Physical[,parameter] + solution[[i]]$s_Physical[,parameter],
                                      colour = shQuote(solution[[i]]$Name)), 
                           width = 0.05)
    
    #  g <- g + geom_point(aes(x = solution$MS_All$Parameters[,x_par], y = rep(0,nrow(solution$MS_All$Parameters))))
  }
  
  g <- g + theme_bw() + theme(axis.text=element_text(size=rel(1.0)))
  
  if (is_legend == FALSE)
  {
    g <- g + guides(colour = "none", shape = "none", fill = "none", linetype = "none")
    
  }
  else {
    g <- g + ggtitle(title)
    
  }
  
  
  return(g)
}