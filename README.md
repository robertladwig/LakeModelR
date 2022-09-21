# LakeModelR
<a href="url"><img src="Figs/mendota.gif" align="right" height="220" width="220" ></a>

-----

:busts_in_silhouette: Robert Ladwig  
:email: [contact](mailto:ladwigjena@gmail.com)
:computer: [more info](https://robertladwig.github.io)

-----

R package to run a modularized 1D integral energy model for water temperature dynamics in a lake. The package can be installed through

```{r gh-installation, eval = FALSE}
#install.packages("remotes")
require(remotes)
remotes::install_github("robertladwig/LakeModelR")
```

Click [here](https://github.com/robertladwig/LakeModeling/blob/main/Manual/1D_IntegralEnergy.pdf) for an overview of the model algorithm.

To test the model code run either example.R or

```{r gh-run, eval = FALSE}
#' Example workflow to run the 1D integral energy model
#' Long-term Mendota data were obtained from North Temperate Lakes Long Term
#' Ecological Research program (#DEB-1440297)
#' @author: Robert Ladwig
#' @email: ladwigjena@gmail.com

## CLEAN WORKSPACE
rm(list = ls())

## LOAD PACKAGE(S)
library(LakeModelR)
require(tidyverse)

## GENERAL LAKE CONFIGURATION
zmax = 25 # maximum lake depth
nx = 25 # number of layers we want to have
dt = 3600  # temporal step (here, one hour because it fits boundary data)
dx = zmax/nx # spatial step

## HYPSOGRAPHY OF THE LAKE
hyps_all <- get_hypsography(hypsofile = system.file('extdata', 'bathymetry.csv',
                            package = 'LakeModelR'),
                            dx = dx, nx = nx)

## ATMOSPHERIC BOUNDARY CONDITIONS
meteo_all <- provide_meteorology(meteofile = system.file('extdata', 'meteorology.csv',
                            package = 'LakeModelR'),
                            secchifile = NULL)

### TIME INFORMATION
startingDate <- meteo_all[[1]]$datetime[1]
startTime = 1
endTime = 365 *24 * 3600 # seconds
total_runtime = endTime / 24 / 3600 # days

# INTERPOLATE ATMOSPHERIC BOUNDARY CONDITIONS
meteo = get_interp_drivers(meteo_all = meteo_all,
                           total_runtime = total_runtime,
                           dt = dt,
                           method = "integrate",
                           secchi = F)

## DEFINE INITIAL WATER TEMPERATURE FROM OBSERVED DATA
u_ini <- initial_profile(initfile = system.file('extdata', 'observedTemp.txt',
                            package = 'LakeModelR'),
                            nx = nx, dx = dx,
                            depth = hyps_all[[2]],
                            processed_meteo = meteo_all[[1]])

## RUN THE LAKE MODEL
res <-  run_thermalmodel(u = u_ini,
                          startTime = startTime,
                          endTime =  endTime,
                          ice = FALSE,
                          Hi = 0,
                          iceT = 6,
                          supercooled = 0,
                          kd_light = 0.5,
                          sw_factor = 1.0,
                          zmax = zmax,
                          nx = nx,
                          dt = dt,
                          dx = dx,
                          area = hyps_all[[1]], # area
                          depth = hyps_all[[2]], # depth
                          volume = hyps_all[[3]], # volume
                          daily_meteo = meteo,
                          Cd = 0.0013,
                          scheme = 'implicit')

## SAVE THE RESULTS
temp = res$temp
mixing = res$mixing
ice = res$icethickness
avgtemp = res$average

## POST-PROCESSING OF THE RESULTS
time =  startingDate + seq(1, ncol(temp), 1) * dt
avgtemp = as.data.frame(avgtemp)
colnames(avgtemp) = c('time', 'epi', 'hyp', 'tot', 'stratFlag', 'thermoclineDep')
avgtemp$time = time

## AVERAGE TEMPERATURES IN EPILIMNION AND HYPOLIMNION
ggplot(avgtemp) +
  geom_line(aes(time, epi, col = 'epilimnion')) +
  geom_line(aes(time, hyp, col = 'hypolimnion')) +
  geom_line(aes(time, tot, col = 'total')) +
  theme_minimal()

## CREATE DATAFRAME FOR FULL TEMPERATURE PROFILES
df <- data.frame(cbind(time, t(temp)) )
colnames(df) <- c("time", as.character(paste0(seq(1,nrow(temp)))))
m.df <- reshape2::melt(df, "time")
m.df$time <- time

## CREATE DATAFRAME FOR ICE
df.ice = data.frame('time' = time,
                    'ice_h' = ice)

## HEATMAP OF WATER TEMPERATURE WITH THERMOCLINE DEPTH AND ICE THICKNESS
ggplot(m.df, aes((time), dx*as.numeric(as.character(variable)))) +
  geom_raster(aes(fill = as.numeric(value)), interpolate = TRUE) +
  scale_fill_gradientn(limits = c(-1,30),
                       colours = rev(RColorBrewer::brewer.pal(11, 'Spectral')))+
  theme_minimal()  +xlab('Time') +
  ylab('Depth [m]') +
  labs(fill = 'Temp [degC]')+
  geom_line(data = avgtemp, aes(time, thermoclineDep, col = 'thermocline depth'), linetype = 'dashed', col = 'brown') +
  geom_line(data = df.ice, aes(time, ice_h * (-1), col = 'ice thickness'), linetype = 'solid', col = 'darkblue') +
  scale_y_reverse()

```
