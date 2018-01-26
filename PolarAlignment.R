

# calculate observed ra and de with given 
#     {1} polar axis azimuth error dA and {2} altitude error dH, 
#     {3} catalugue ra on date (ac, in hours) and {4} catalogue de on date (dc, degree),
#     {5} sideral time on the moment, {6} phi - observation site latitude (degree)
#  returns marix with two columns: {1} ra_obs (in hours), {2} de_obs (in degrees)
Get_RaDeObserved <- function (dA, dH, ac, dc, Ts, phi)
{
  a_o <- ac - (dH/15) * sin((Ts-ac)*15*pi/180) * tan(dc*pi/180) - (dA/15) * (sin(phi*pi/180) - cos((Ts-ac)*15*pi/180) * cos(phi*pi/180) * tan(dc*pi/180))
  d_o <- dc + (dH)    * cos((Ts-ac)*15*pi/180)                 +  (dA)    * cos(phi*pi/180) *sin((Ts-ac)*15*pi/180)
  return(cbind(a_o, d_o))
 }


# calulate azimuth and altitude error of polar axis
# Input: 
#   observations - array of data, each row - 1 observation, 
#       columns: {1} ra_obs (hours), {2} de_obs (degree), {3} ra_cat (hours), {4} de_cat (degree), {5} sideral time (hours)
#   phi - observation site latitude (degree)
Get_dADdH <- function(observations, phi)
{
  d_a <-  (observations[,1] - observations[,3])*15
  d_d <-  (observations[,2] - observations[,4])
  
  y <- rbind(as.matrix(d_a, nrow(observations), 1), as.matrix(d_d, nrow(observations), 1))
  
  print(y)
  
  x1 <- matrix(0, nrow(observations), 2)
  x1[,1] <- -sin(phi*pi/180) + cos((observations[,5]-observations[,3])*15*pi/180)*cos(phi*pi/180)*tan(observations[,4]*pi/180)
  x1[,2] <- -sin((observations[,5]-observations[,3])*15*pi/180)*tan(observations[,4]*pi/180)
  #x1 <- x1*15
  
  x2 <- matrix(0, nrow(observations), 2)
  x2[,1] <- cos(phi*pi/180)*sin((observations[,5]-observations[,3])*15*pi/180)
  x2[,2] <- cos((observations[,5]-observations[,3])*15*pi/180)
  
  x <- rbind(x1, x2)
  
  print(x)
  
  res <- TLS_Gen(In = x, y = y, mode = 1, scaling = 0, ef = ncol(x))
  
  return(res);
}


# calculation test
#   stars - matrix with 3 columns, 
#    {1}  catalogue ra (hours) 
#    {2}  cataloge de (hours)
#    {3}  sideral time (hours)
#    each row - one observation
#   da, dh - azimuth and altitude errors of polar alignment (degrees)
#   phi - observation site latitude (degree)
#   s0 - noise level in standard deviations (degrees)
# example:
# dAdH_test(stars_ah[,3:5], da = 0.167, dh = -0.167, s0 = 0.01)
dAdH_test <- function(stars, da = 0.0167, dh = 0.0167, phi = 59.8388, s0 = 0)
{
  print("Initial error:")
  print(c(da, dh))
  observations <- matrix(0, nrow = nrow(stars), ncol = 5)
  observations[,1:2] <- Get_RaDeObserved(da, dh, stars[,1], stars[,2], stars[,3], phi)
  observations[,3:5] <- stars[,1:3]
  print("Observations:")
  print(observations)
  
  print("Noise:")
  print(s0)
  noise <- matrix(rnorm(n = nrow(stars)*2, mean = 0, sd = s0), nrow = nrow(stars), ncol = 2)
  print(noise)
  observations[,1:2] <- observations[,1:2] + noise
  
  
  res <- Get_dADdH(observations = observations, phi = phi)
  
  print("Estimated errors:")
  print(res$X)
  print(res$s_X)
  return(res)
  
}