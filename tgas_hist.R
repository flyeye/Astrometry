

draw_tgas_hist_Z <- function(tgas_, name_){
  tgas_ <- CalcGalXYZ(tgas_)
  g<- ggplot(data = tgas_) + geom_histogram(aes(x = z), fill = "gray70", colour = "gray10") + 
    scale_x_continuous(breaks=seq(-1,1,by=0.1), minor_breaks=seq(-1,1,by=0.025), limits = c(-1,1)) #+ 
    #scale_y_continuous(breaks=seq(0,250000,by=25000), minor_breaks=seq(0,250000,by=12500)) + 
    theme_bw()
  
  ggsave(filename = paste0(name_, ".jpeg"), width = 10, height = 10)
         
  return(g);
}


tgas_make_z_hist <- function(tgas)
{
  tgas_ <- tgas %>% filter(!is.na(B_V) & !is.na(M) & (LClass_apass == 5))
  draw_tgas_hist_Z(tgas_, name_ = "TGAS_MS_Z_all_stars")
  
  tgas_ <- tgas %>% filter(!is.na(B_V) & !is.na(M) & (LClass_apass == 5) & ((B_V>1.42) & (B_V<Inf)))
  draw_tgas_hist_Z(tgas_, name_ = "TGAS_MS_Z_M-stars")
                  
  tgas_ <- tgas %>% filter(!is.na(B_V) & !is.na(M) & (LClass_apass == 5) & ((B_V>0.85) & (B_V<1.42)))
  draw_tgas_hist_Z(tgas_, name_ = "TGAS_MS_Z_K-stars")
  
  tgas_ <- tgas %>% filter(!is.na(B_V) & !is.na(M) & (LClass_apass == 5) & ((B_V>0.58) & (B_V<0.85)))
  draw_tgas_hist_Z(tgas_, name_ = "TGAS_MS_Z_G-stars")
  
  tgas_ <- tgas %>% filter(!is.na(B_V) & !is.na(M) & (LClass_apass == 5) & ((B_V>0.29) & (B_V<0.58)))
  draw_tgas_hist_Z(tgas_, name_ = "TGAS_MS_Z_F-stars")
  
  tgas_ <- tgas %>% filter(!is.na(B_V) & !is.na(M) & (LClass_apass == 5) & ((B_V>0.0) & (B_V<0.29)))
  draw_tgas_hist_Z(tgas_, name_ = "TGAS_MS_Z_A-stars")
  
  tgas_ <- tgas %>% filter(!is.na(B_V) & !is.na(M) & (LClass_apass == 5) & ((B_V>(-0.3)) & (B_V<0.0)))
  draw_tgas_hist_Z(tgas_, name_ = "TGAS_MS_Z_B-stars")
  
  tgas_ <- tgas %>% filter(!is.na(B_V) & !is.na(M) & (LClass_apass == 5)) %>% filter(B_V>(-Inf)) %>% filter(B_V<(-0.3))
  draw_tgas_hist_Z(tgas_, name_ = "TGAS_MS_Z_O-stars")
  
}


