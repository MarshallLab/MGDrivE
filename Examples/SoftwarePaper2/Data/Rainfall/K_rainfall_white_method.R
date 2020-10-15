### White et al. 2011 Rainfall Dynamics 
### Tomás M. León 2020

## Comoros Data

km_precip <- read.csv("km_precip_2010s_daily.csv", header = T) #Data from 2010s (daily rainfall)

## White Methods (see paper):
# i)	Carrying capacity proportional to mean past rainfall
# ii)	Carrying capacity proportional to linearly weighted past rainfall
# iii)	Carrying capacity proportional to exponentially weighted past rainfall


tau <- 4 # number of past days of rainfall to account for (White et al. found 4 to be best for their model)
lambda <- 1 #fitted scaling factor (differs by site; fit to data in White et al.); in our case we rescale after the fact
L_eq <- 120000 #In the situation of the Comoros, this is the maximum value that we allow carrying capacity to attain
constant_frac <- 2/3 #In the situation of the Comoros, some fraction of K is ever present and constant

K_km_out <- array(NA, dim = c(nrow(km_precip)-tau+1,3,3))
for(j in 1:dim(K_km_out)[2]){
  rain <- km_precip[,j]
  for(i in tau:length(rain)){
    K_km_out[i-tau+1,j,1] <- lambda / tau * sum(rain[(i-tau+1):i])
    K_km_out[i-tau+1,j,2] <- lambda / (tau^2) * sum(rain[(i-tau+1):i]*(1:tau))
    K_km_out[i-tau+1,j,3] <- lambda / (tau * (1 - exp(-i/tau))) * sum(rain[1:i]*exp(-(i-(1:i))/tau))
  }
}

#dimnames(K_km_out) <- list(1:dim(K_km_out)[1], colnames(km_precip), c("Method 1", "Method 2", "Method 3"))

K_km_out_method3 <- as.data.frame(K_km_out[,,3]) #White et al. found Method 3 was best for their purposes

colnames(K_km_out_method3) <- colnames(km_precip)

K_km_out_method3 <- K_km_out_method3 / max(K_km_out_method3) * L_eq * (1 - constant_frac) + L_eq * constant_frac #Rescale for reasonable values based on input parameters

K_km_out_method3 <- cbind(Day = 1:dim(K_km_out_method3)[1], K_km_out_method3)

#readr::write_csv(K_km_out_method3, path = "K_km_v2.csv")

K_km_plot <- reshape2::melt(K_km_out_method3, id.vars = "Day")
colnames(K_km_plot) <- c("Day", "Island", "K")

ggplot(K_km_plot) + geom_line(aes(x = Day, y = K, col = Island)) + theme_bw() + ggtitle("Carrying Capacity in Comoros") +
  theme(plot.title = element_text(hjust = 0.5)) + expand_limits(y = 0)
