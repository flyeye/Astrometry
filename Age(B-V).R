

ggplot() + 
  scale_x_continuous(breaks = seq(0.1, 1.6, by = 0.1), limits = c(0.2, 1.6)) + 
  scale_y_continuous(breaks = seq(1, 9, by = 1), limits = c(1,9)) + 
  geom_line(aes(x = stars_ph_prop$B.V, y = stars_ph_prop$Age)) +
  xlab("B-V") + ylab("Возраст, миллиардов лет")