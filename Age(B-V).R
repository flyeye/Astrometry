

g <- ggplot() + 
  scale_x_continuous(breaks = seq(0.3, 1.5, by = 0.3), limits = c(0.2, 1.6)) + 
  scale_y_continuous(breaks = seq(1, 9, by = 1), limits = c(1,9)) + 
  geom_line(aes(x = stars_ph_prop$B.V, y = stars_ph_prop$Age)) +
  xlab("B-V") + ylab("Возраст, миллиардов лет") +
  theme_bw() + 
  geom_vline(xintercept=c(0.3, 0.42, 0.58, 0.69, 0.85, 1.16, 1.42), linetype = "dashed")
ggsave(filename = "Age_B-V.jpeg", plot = g, width = 4, height = 4)
g


g <- ggplot() + 
  scale_x_reverse(breaks = seq(7500, 4000, by = -500), limits = c(7500,3500)) + 
  scale_y_continuous(breaks = seq(1, 9, by = 1), limits = c(1,9)) + 
  geom_line(aes(x = stars_ph_prop$Teff, y = stars_ph_prop$Age)) +
  xlab("Тэфф") + ylab("Возраст, миллиардов лет") +
  geom_vline(xintercept=c(7078, 6540, 5925, 5564, 5122, 4456, 4011), linetype = "dashed") + 
  theme_bw()  
 
ggsave(filename = "Age_Teff.jpeg", plot = g, width = 4, height = 4)
g