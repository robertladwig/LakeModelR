## function to calculate density from temperature
calc_dens <-function(wtemp){
  dens = 999.842594 + (6.793952 * 1e-2 * wtemp) - (9.095290 * 1e-3 *wtemp**2) +
    (1.001685 * 1e-4 * wtemp**3) - (1.120083 * 1e-6* wtemp**4) +
    (6.536336 * 1e-9 * wtemp**5)
  return(dens)
}

## this is our attempt for turbulence closure, estimating eddy diffusivity
eddy_diffusivity <-function(rho, depth, g, rho_0, ice, area){
  buoy = rep(1, (nx)) * 7e-5
  buoy[1:(nx-1)] = abs(rho[2:nx] - rho[1:(nx-1)]) / (depth[2:nx] - depth[1:(nx-1)]) * g/rho_0
  # for (i in seq(1, nx-1)){#range(0, nx - 1):
  #   buoy[i] = ( abs(rho[i+1] - rho[i]) / (depth[i+1] - depth[i]) * g/rho_0 )
  # }
  buoy[nx] = ( abs(rho[nx-1] - rho[nx]) / abs(depth[nx-1] - depth[nx]) *
                 g/rho_0 )

  low_values_flags = buoy < 7e-5  # Where values are low
  buoy[low_values_flags] = 7e-5

  if (ice){
    ak <- 0.000898
  } else{
    ak <- 0.00706 *( max(area)/1E6)**(0.56)
  }

  kz = ak * (buoy)**(-0.43)
  return(kz)
}

provide_meteorology <- function(meteofile, secchifile = NULL,
                                windfactor = 1.0){
  meteo <- read_csv(meteofile)

  daily_meteo <- meteo
  daily_meteo$date = daily_meteo$datetime
  daily_meteo$Cloud_Cover <- calc_cc(date = as.POSIXct(daily_meteo$date),
                                                airt = daily_meteo$Air_Temperature_celsius,
                                                relh = daily_meteo$Relative_Humidity_percent,
                                                swr = daily_meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared,
                                                lat = 40, lon = -105,
                                                elev = 3110)
  daily_meteo$dt <- as.POSIXct(daily_meteo$date) - (as.POSIXct(daily_meteo$date)[1]) + 1
  daily_meteo$ea <- (daily_meteo$Relative_Humidity_percent * (4.596 * exp((17.27*(daily_meteo$Air_Temperature_celsius))/
                                                                            (237.3 + (daily_meteo$Air_Temperature_celsius) )))/100)
  daily_meteo$ea <- (101.325 * exp(13.3185 * (1 - (373.15 / (daily_meteo$Air_Temperature_celsius + 273.15))) -
                                     1.976 * (1 - (373.15 / (daily_meteo$Air_Temperature_celsius + 273.15)))**2 -
                                     0.6445 * (1 - (373.15 / (daily_meteo$Air_Temperature_celsius + 273.15)))**3 -
                                     0.1229 * (1 - (373.15 / (daily_meteo$Air_Temperature_celsius + 273.15)))**4)) *daily_meteo$Relative_Humidity_percent/100
  daily_meteo$ea <- (daily_meteo$Relative_Humidity_percent/100) * 10^(9.28603523 - 2322.37885/(daily_meteo$Air_Temperature_celsius + 273.15))
  startDate <- daily_meteo$datetime[1]

  ## calibration parameters
  daily_meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared <-
    daily_meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared
  daily_meteo$Ten_Meter_Elevation_Wind_Speed_meterPerSecond <-
    daily_meteo$Ten_Meter_Elevation_Wind_Speed_meterPerSecond * windfactor# wind speed multiplier

  # kd = 0.4 # 0.2# 1.0 #0.2 # light attenuation coefficient

  ## light
  # Package ID: knb-lter-ntl.31.30 Cataloging System:https://pasta.edirepository.org.
  # Data set title: North Temperate Lakes LTER: Secchi Disk Depth; Other Auxiliary Base Crew Sample Data 1981 - current.
  if (!is.null(secchifile)){
    secview <- read_csv(secchifile) %>%
      dplyr::filter(sampledate >= startDate)
    if (secview$sampledate[1] >= startDate){
      secview <- rbind(data.frame('sampledate' = startDate,
                                  'secnview' = secview$secnview[1]),
                       secview)
    }
    secview$dt <- as.POSIXct(secview$sampledate) - (as.POSIXct(secview$sampledate)[1]) + 1
    secview$kd <- 1.7 / secview$secnview
    secview$kd  <- zoo::na.approx(secview$kd)
  } else {
    secview = NULL
  }


  return(list(daily_meteo, secview))
}

# FROM GOTMTOOLS, original author: Tadhg Moore
calc_cc <- function(date, airt, relh = NULL, dewt = NULL, swr, lat, lon, elev, daily = F){
  orig_date = date
  timestep = difftime(orig_date[2], orig_date[1], units = "secs")

  # If the time step is 24 hours or more, create artificial hourly time steps
  if(timestep >= as.difftime(24, units = "hours")){
    date = seq.POSIXt(from = date[1], to = (date[length(date)] + timestep - 1 * 60 * 60), by = '1 hour')
  }

  yday <- lubridate::yday(date)
  hour <- lubridate::hour(date)
  hour[hour == 0] <- 24

  std.mer = seq(-90,90, 15)
  Lsm = std.mer[which.min(abs(lon - std.mer))] # Local standard meridian (degrees)

  Hsc = 1390 # Solar constant (W/m2)
  cd = 0.06 # Dust coefficient
  Rg = 0.045 # Reflectivity of the ground - extended mixed forest


  theta = lat*pi/180 # Latitude in radians

  r = 1 + 0.017 * cos((2*pi/365)*(186-yday)) # Relative earth-sun distance

  d = 23.45 * pi/180 * cos((2*pi/365)*(172-yday)) # Declination of the sun

  dts = (1/15) * (Lsm-lon) # Fraction of 15-degree increment Llm is east of Lsm
  value = (sin(theta)*sin(d))
  value = value/(cos(theta)*cos(d))
  tss = (12/pi) * acos(-value) + dts + 12 # Time of sunset
  tsu = -tss + (2 * dts) + 24 # Time of sunrise

  gamma = rep(0, length(tss)) # Correction factor
  dum = which(hour>tsu & hour<tss)
  gamma[dum] = 1

  #Calculate Hb and Htheta
  dum1 = which(hour <=12 )
  dum2 = which(hour > 12 )
  hb1  = pi/12*(hour-1-dts)
  hb1[dum1] = hb1[dum1]+pi
  hb1[dum2] = hb1[dum2]-pi
  hb  = hb1
  dum3 = which(hb1 > 2*pi)
  hb[dum3] = hb[dum3] - 2 * pi
  dum4 = which(hb1 < 0)
  hb[dum4] = hb[dum4] + 2 * pi
  #rm(c(dum3, dum4))
  he1  = pi/12*(hour-dts)
  he1[dum1] = he1[dum1]+pi
  he1[dum2] = he1[dum2]-pi
  he  = he1
  dum3 = which(he1 > 2*pi)
  he[dum3] = he[dum3] - 2*pi
  dum4 = which(he1 < 0)
  he[dum4] = he[dum4] + 2*pi
  #clear dum1 dum2 dum3 dum4

  Ho = Hsc/(r^2)*(sin(theta)*sin(d)+12/pi*cos(theta)*cos(d)*(sin(he)-sin(hb)))*gamma

  # Radiation scattering and absorption #####################################

  w = (he+hb)/2 # Hour angle
  alpha1 = abs(sin(theta)*sin(d)+cos(theta)*cos(d)*cos(w))
  alpha = atan(alpha1/sqrt(1-alpha1^2)) # Solar altitude

  theta_am1 = ((288-0.0065*elev)/288)^5.256
  theta_am2 = sin(alpha)+0.15*((alpha*180/pi)+3.855)^(-1.253)
  theta_am = theta_am1/theta_am2 # Optical air mass

  # Dewpoint temperature
  if(is.null(dewt)){
    dewt <- 243.04*(log(relh/100)+((17.625*airt)/(243.04+airt)))/(17.625-log(relh/100)-((17.625*airt)/(243.04+airt)))
  }
  if(timestep >= as.difftime(2, units = "hours")){
    dewt = rep(dewt, each = as.numeric(difftime(orig_date[2], orig_date[1], units = "hours")))
  }

  Pwc = 0.85*exp(0.11+0.0614*dewt) # Precipitable atmospheric water content

  a2 = exp(-(0.465+0.134*Pwc)*(0.179+0.421*exp(-0.721*theta_am))*theta_am) # Atmospheric transmission coefficient after scattering and absorption
  a1 = exp(-(0.465+0.134*Pwc)*(0.129+0.171*exp(-0.88*theta_am))*theta_am)
  at = (a2+0.5*(1-a1-cd))/(1-0.5*Rg*(1-a1-cd)) # attenuation (scattering and absorption)
  #att = mean(at)

  Ho = at*Ho
  #Ho = att*Ho

  dum5 = which(Ho<0.0)
  Ho[dum5] = 1

  df = data.frame(DateTime = date,Ho = Ho)
  if(timestep >= as.difftime(2, units = "hours")){
    df = aggregate(list(Ho = df$Ho), by = list(DateTime = cut(df[,1], paste(timestep, "s"))), mean, na.rm = T)
  }
  df$swr =swr

  df$ccsim <- NA

  for(i in 1:nrow(df)){
    if(df$Ho[i] < df$swr[i]){
      df$ccsim[i] <- NaN
    }else{
      df$ccsim[i] <- sqrt((1 - (df$swr[i]/df$Ho[i]))/0.65)
    }
  }

  ccsim = df$ccsim
  ccsim[ccsim > 1] <- 1

  # Fill gaps with the mean value of the previous and posterior.
  sta = min(which(!is.nan(ccsim)))
  stp = max(which(!is.nan(ccsim)))
  ccsim[sta:stp] <- zoo::na.approx(ccsim[sta:stp])
  if(sta != 1){
    ccsim[1:sta] <- ccsim[sta]
  }
  if(stp != length(ccsim)){
    ccsim[stp:length(ccsim)] <- ccsim[stp]
  }

  return(ccsim)
}

initial_profile <- function(initfile, nx, dx, depth, processed_meteo){
  meteo <- processed_meteo
  startDate <- meteo$datetime[1]
  obs <- read_csv(initfile)
  init.df <- obs %>%
    mutate(ditt = as.numeric(abs((datetime) - (startDate)))) %>%
    dplyr::filter(ditt == min(ditt)) %>%
    arrange(Depth_meter)
  if (max(depth) > max(init.df$Depth_meter)){
    init.df <- rbind(init.df, init.df[nrow(init.df),])
    init.df$Depth_meter[nrow(init.df)] <- max(depth)
  }
  u = approx(init.df$Depth_meter, as.numeric(unlist(init.df[,3])),
             seq(0, nx * dx, length.out= nx), rule = 2)$y
  warning(paste0('Meteorological starting date is ',as.Date(startDate),', but observed data starts ',min(init.df$ditt),' days later on ',
                 as.Date(min(init.df$datetime))))
  return(u)
}

get_hypsography <- function(hypsofile, dx, nx){
  hyps <- read_csv(hypsofile)
  area = approx(hyps$Depth_meter,hyps$Area_meterSquared,seq(1,nx*dx,
                                                            length.out= nx))$y
  # area[which.min(area)] <- 1e-2
  area[nx] <- area[nx-1] -1
  depth = seq(1,nx*dx, length.out = nx)
  volume <- c(rev(diff(pracma::cumtrapz(area, depth))*(-1)),0)
  volume[which(volume == 0)] = min(volume[-which(volume == 0)])
  volume <- rep(0, (length(depth)-1))
  for (p in 1:length(volume)){
    volume[p] <- pracma::trapz(depth[p:(p+1)],area[p:(p+1)])
  }
  # volume <- c(volume, 1000)
  volume <- c(volume, volume[length(volume)])
  return(list(area, depth, volume))
}

longwave <- function(cc, sigma, Tair, ea, emissivity, Jlw){  # longwave radiation into
  Tair = Tair + 273.15
  p <- (1.33 * ea/Tair)
  Ea <- 1.24 * (1 + 0.17 * cc**2) * p**(1/7)
  lw <- emissivity * Ea *sigma * Tair**4
  return(lw)
}
backscattering <- function(emissivity, sigma, Twater, eps){ # backscattering longwave
  # radiation from the lake
  Twater = Twater + 273.15
  back = (eps * sigma * (Twater )^4)
  return((-1) * back)
}
# sensible <- function(p2, B, Tair, Twater, Uw){ # convection / sensible heat
#   Twater = Twater + 273.15
#   Tair = Tair + 273.15
#   fu = 4.4 + 1.82 * Uw + 0.26 *(Twater - Tair)
#   sensible <- ( p2 * B * fu * (Twater - Tair))
#   return((-1) * sensible)
# }

PSIM <- function(zeta){
  # Function to compute stability functions for momentum
  if (zeta < 0.0){
    X = (1 - 16*zeta)^0.25;
    psim = 2*log((1 + X)/2) + log((1 + X*X)/2)-2*atan(X) + pi/2;
  }
  else if(zeta > 0.0){
    if (zeta > 0.5){
      if (zeta > 10.0){
        psim = log(zeta) - 0.76*zeta - 12.093;
      }
      else{
        psim = 0.5/(zeta*zeta) - 4.25/zeta - 7.0*log(zeta) - 0.852;
      }
    }
    else {
      psim = -5*zeta ;
    }

  }   # Stable case
  else {
    psim = 0.0;
  }

  return(psim)
}



PSITE = function(zeta){
  # Function to compute stability functions for sensible and latent heat
  if (zeta < 0.0){
    X = (1 - 16*zeta)^0.25;
    psite = 2*log((1 + X*X)/2);
  }
  else if (zeta > 0.0)  { # Stable case
    if (zeta > 0.5)   {
      if (zeta > 10.0) {
        psite = log(zeta) - 0.76*zeta - 12.093;
      }
      else{
        psite = 0.5/(zeta*zeta) - 4.25/zeta - 7.0*log(zeta) - 0.852;
      }
    }
    else {
      psite = -5*zeta ;
    }
  }
  else {
    psite = 0.0;
  }
  return(psite)
}
sensible <- function(Tair, Twater, Uw, p2, pa, ea, RH, A, Cd = 0.013){ # evaporation / latent heat
  # https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2009JD012839

  # Tair =0
  # Twater = 0
  # Uw = 0.01
  # pa = 98393
  # ea = 6.079572
  # A = 31861
  # Cd = 0.0037
  #


  const_SpecificHeatAir = 1005;           # Units : J kg-1 K-1
  const_vonKarman = 0.41;                 # Units : none
  const_Gravity = 9.81;                   # Units : m s-2
  const_Charnock = Cd;

  U_Z = Uw
  if (Uw <= 0){
    U_Z = 1e-3
  }
  T = Tair
  if (Tair == 0){
    T = runif(1, min = 1e-7, max = 1e-5)
  }
  T0 = Twater
  if (Twater == 0){
    T0 = runif(1, min = 1e-7, max = 1e-5)
  }
  Rh=RH
  p=pa/100
  z=2

  # Step 2c - Compute saturated vapour pressure at air temperature
  e_s = 6.11*exp(17.27*T/(237.3+T)); # Units : mb ##REF##
  # Step 2d - Compute vapour pressure
  e_a = Rh*e_s/100; # Units : mb
  ### End step 2

  ### Step 3 - Compute other values used in flux calculations
  # Step 3a - Compute specific humidity
  q_z = 0.622*e_a/p; # Units: kg kg-1
  # Step 3b - Compute saturated vapour pressure at water temperature
  e_sat = 6.11*exp(17.27*T0/(237.3+T0)); # Units : mb ##REF##
  # Step 3c - Compute humidity at saturation (Henderson-Sellers 1986 eqn 36)
  q_s = 0.622*e_sat/p; # Units: kg kg-1
  # Step 3d - Compute latent heat of vaporisation
  L_v = 2.501e6-2370*T0; # Units : J kg-1 ** EQUATION FROM PIET ##REF##
  # Step 3e - Compute gas constant for moist air
  R_a = 287*(1+0.608*q_z); # Units : J kg-1 K-1
  # Step 3f - Compute air density
  rho_a = 100*p/(R_a*(T+273.16)); # Units : kg m-3
  # Step 3g - Compute kinematic viscosity of air
  v = (1./rho_a)*(4.94e-8*T + 1.7184e-5); # Units : m2 s-1
  # Step 3h - Compute virtual air temperature and virtual air-water temperature difference
  T_v = (T+273.16)*(1+0.61*q_z); # Units - K
  T_ov = (T0+273.16)*(1+0.61*q_s); # Units - K
  del_theta = T_ov - T_v;
  # Step 3h - Compute water density
  rho_w = 1000*(1-1.9549*0.00001*abs(T0-3.84)^1.68);
  ### End step 3

  # step 4
  u_star = U_Z *sqrt(0.00104+0.0015/(1+exp((-U_Z+12.5)/1.56))); # Amorocho and DeVries, initialise ustar using U_Z

  if (u_star == 0){
    u_star = 1e-6
  }

  z_0 = (const_Charnock*u_star^2./const_Gravity) + (0.11*v/u_star);
  z_0_prev=z_0*1.1; # To initiate the iteration
  for (i1 in 1:length(U_Z)){
    while (abs((z_0[i1] - z_0_prev[i1]))/abs(z_0_prev[i1]) > 0.000001){ # Converge when z_0 within 0.0001# of previous value
      u_star[i1]=const_vonKarman*U_Z[i1]/(log(z/z_0[i1]));  # Compute u_star
      dummy = z_0[i1]; # Used to control while loop
      z_0[i1]=(const_Charnock*u_star[i1]^2./const_Gravity) + (0.11*v[i1]/u_star[i1]); # Compute new roughness length
      z_0_prev[i1] = dummy; # Used to control while loop
    }
  }

  # Step 4d - Compute initial neutral drag coefficient
  C_DN = (u_star^2)/(U_Z^2); # Units - none
  # Step 4e - Compute roughness Reynolds number
  Re_star = u_star*z_0/v; # Units - none
  # Step 4f - Compute initial roughness length for temperature
  z_T = z_0*exp(-2.67*(Re_star)^(1/4) + 2.57); # Units - m
  z_T = Re(z_T); # Get real components, and NaN can create imag component despite no data
  # Step 4g - Compute initial roughness length for vapour
  z_E = z_0*exp(-2.67*(Re_star)^(1/4) + 2.57); # Units - m
  z_E = Re(z_E); # Get real components, and NaN can create imag component despite no data
  # Step 4h - Compute initial neutral sensible heat transfer coefficient
  C_HN = const_vonKarman*sqrt(C_DN)/(log(z/z_T));
  # Step 4i - Compute initial neutral latent heat transfer coefficient
  C_EN = const_vonKarman*sqrt(C_DN)/(log(z/z_E));
  ### End step 4

  ### Step 5 - Start iteration to compute corrections for atmospheric stability
  for (i1 in 1:length(U_Z)){
    # Step 5a - Compute initial sensible heat flux based on neutral coefficients
    H_initial = rho_a[i1]*const_SpecificHeatAir*C_HN[i1]*U_Z[i1]*(T0[i1]-T[i1]); # Units : W m-2
    # Step 5b - Compute initial latent heat flux based on neutral coefficients
    E_initial = rho_a[i1]*L_v[i1]*C_EN[i1]*U_Z[i1]*(q_s[i1]-q_z[i1]); # Units : W m-2
    # Step 5c - Compute initial Monin-Obukhov length
    L_initial = (-rho_a[i1]*u_star[i1]^3*T_v[i1])/(const_vonKarman*const_Gravity*(H_initial[i1]/const_SpecificHeatAir + 0.61*E_initial[i1]*(T[i1]+273.16)/L_v[i1])); # Units - m
    # Step 5d - Compute initial stability parameter
    zeta_initial = z/L_initial[i1];
    # Step 5e - Compute initial stability function
    psim=PSIM(zeta_initial[i1]); # Momentum stability function
    psit=PSITE(zeta_initial[i1]); # Sensible heat stability function
    psie=PSITE(zeta_initial[i1]); # Latent heat stability function
    # Step 5f - Compute corrected coefficients
    C_D=const_vonKarman*const_vonKarman/(log(z/z_0[i1])-psim[i1])^2;
    C_H=const_vonKarman*sqrt(C_D[i1])/(log(z/z_T[i1])-psit[i1]);
    C_E=const_vonKarman*sqrt(C_D[i1])/(log(z/z_E[i1])-psie[i1]);
    # Step 5g - Start iteration
    L_prev = L_initial[i1];
    L = L_prev*1.1; # Initialise while loop
    count=0;
    while (abs((L[i1] - L_prev))/abs(L_prev) > 0.000001){
      # Iteration counter
      count[i1]=count[i1]+1;
      if (count[i1] > 20){
        break;
      }
      # Step 5i - Compute new z_O, roughness length for momentum
      z_0= (const_Charnock*u_star[i1]^2./const_Gravity) + (0.11*v[i1]/u_star[i1]);
      # Step 5j - Compute new Re_star
      Re_star = u_star[i1]*z_0[i1]/v[i1];
      # Step 5k - Compute new z_T, roughness length for temperature
      z_T = z_0[i1]*exp(-2.67*(Re_star[i1])^(1/4) + 2.57);
      # Step 5l - Compute new z_E, roughness length for vapour
      z_E = z_0[i1]*exp(-2.67*(Re_star[i1])^(1/4) + 2.57);
      # Step 5p - Compute new stability parameter
      zeta = z/L[i1];
      #fprintf('zeta #g\n',zeta[i1]);
      # Step 5q - Check and enforce bounds on zeta
      if (zeta[i1] > 15){
        zeta[i1] = 15}
      else if (zeta[i1] < -15) {
        zeta[i1] = -15}
      # Step 5r - Compute new stability functions
      psim=PSIM(zeta[i1]); # Momentum stability function
      psit=PSITE(zeta[i1]); # Sensible heat stability function
      psie=PSITE(zeta[i1]); # Latent heat stability function
      # Step 5s - Compute corrected coefficients
      C_D=const_vonKarman*const_vonKarman/(log(z/z_0[i1])-psim[i1])^2;
      C_H=const_vonKarman*sqrt(C_D[i1])/(log(z/z_T[i1])-psit[i1]);
      C_E=const_vonKarman*sqrt(C_D[i1])/(log(z/z_E[i1])-psie[i1]);
      # Step 5m - Compute new H (now using corrected coefficients)
      H = rho_a[i1]*const_SpecificHeatAir*C_H[i1]*U_Z[i1]*(T0[i1]-T[i1]);
      # Step 5n - Compute new E (now using corrected coefficients)
      E = rho_a[i1]*L_v[i1]*C_E[i1]*U_Z[i1]*(q_s[i1]-q_z[i1]);
      # Step 5h - Compute new u_star
      u_star=sqrt(C_D[i1]*U_Z[i1]^2);
      # Step 5o - Compute new Monin-Obukhov length
      dummy = L[i1]; # Used to control while loop
      L = (-rho_a[i1]*u_star[i1]^3*T_v[i1])/(const_vonKarman*const_Gravity*(H[i1]/const_SpecificHeatAir + 0.61*E[i1]*(T[i1]+273.16)/L_v[i1]));
      L_prev = dummy; # Used to control while loop
    } # Converge when L within 0.0001# or previous L

  } # Need to iterate separately for each record


  ### End step 5

  # Take real values to remove any complex values that arise from missing data or NaN.
  C_D=Re(C_D);
  C_E=Re(C_E);
  C_H=Re(C_H);
  z_0=Re(z_0);
  z_E=Re(z_E);
  z_T=Re(z_T);

  # Compute evaporation [mm/day]
  Evap = 86400*1000*E/(rho_w*L_v);

  sensible = H
  return( sensible * (-1))
}
# latent <- function(Tair, Twater, Uw, p2, pa, ea, RH, A, Cd = 0.013){ # evaporation / latent heat
# Goudsmit
#   Twater = Twater + 273.15
#   Tair = Tair + 273.15
#   Pressure = pa / 100
#   fu = 4.4 + 1.82 * Uw + 0.26 *(Twater - Tair)
#   fw = 0.61 * (1 + 10^(-6) * Pressure * (4.5 + 6 * 10^(-5) * Twater**2))
#   ew = fw * 10 * ((0.7859+0.03477* Twater)/(1+0.00412* Twater))
#   latent = fu * p2 * (ew - ea) # * 1.33) #* 1/6
#   return((-1) * latent * 0.8)
# }
# latent <- function(Tair, Twater, Uw, p2, pa, ea, RH, A, Cd = 0.013){ # evaporation / latent heat HOSTETLER
#   # https://ascelibrary.org/doi/epdf/10.1061/JYCEAJ.0004424
#
#   A = max(A)
#   Lv = 2264 # J/g
#   Pressure = pa / 100
#   Nmt = 3.367 * 10^(-9) * A^(-0.05)
#   tk = 1 - (373.15 / (Twater + 273.15))
#   tka =  1 - (373.15 / (Tair + 273.15))
#   eo = 101.325 * exp(13.3185 * tk - 1.976 * tk^2 - 0.6445 * tk^3 - 0.1229 * tk^4)
#   es = 101.325 * exp(13.3185 * tka - 1.976 * tka^2 - 0.6445 * tka^3 - 0.1229 * tka^4) * RH/100
#
#   latent = Lv * Nmt * Uw * (eo - es)
#
#   b0 = 1.84 * 10^(-7)
#
#   Tvwater = Twater / (1 - 0.378 * es /Pressure)
#   Tvair = Tair / (1 - 0.378 * ea/ Pressure)
#
#   Tvirtual = (Tvwater - Tvair)
#   if (Tvirtual < 0){
#     Tvirtual = 0
#   }
#
#   latent = (b0 * (Tvirtual)^(1/3) + Nmt * Uw) * (eo - es)
#
#   return((-1) * latent)
# }
latent <- function(Tair, Twater, Uw, p2, pa, ea, RH, A, Cd = 0.013){ # evaporation / latent heat
  # https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2009JD012839

  # Tair =0
  # Twater = 0
  # Uw = 0.01
  # pa = 98393
  # ea = 6.079572
  # A = 31861
  # Cd = 0.0037
  #


  const_SpecificHeatAir = 1005;           # Units : J kg-1 K-1
  const_vonKarman = 0.41;                 # Units : none
  const_Gravity = 9.81;                   # Units : m s-2
  const_Charnock = Cd;

  U_Z = Uw
  if (Uw <= 0){
    U_Z = 1e-3
  }
  T = Tair
  if (Tair == 0){
    T = runif(1, min = 1e-7, max = 1e-5)
  }
  T0 = Twater
  if (Twater == 0){
    T0 = runif(1, min = 1e-7, max = 1e-5)
  }
  Rh=RH
  p=pa/100
  z=2

  # Step 2c - Compute saturated vapour pressure at air temperature
  e_s = 6.11*exp(17.27*T/(237.3+T)); # Units : mb ##REF##
  # Step 2d - Compute vapour pressure
  e_a = Rh*e_s/100; # Units : mb
  ### End step 2

  ### Step 3 - Compute other values used in flux calculations
  # Step 3a - Compute specific humidity
  q_z = 0.622*e_a/p; # Units: kg kg-1
  # Step 3b - Compute saturated vapour pressure at water temperature
  e_sat = 6.11*exp(17.27*T0/(237.3+T0)); # Units : mb ##REF##
  # Step 3c - Compute humidity at saturation (Henderson-Sellers 1986 eqn 36)
  q_s = 0.622*e_sat/p; # Units: kg kg-1
  # Step 3d - Compute latent heat of vaporisation
  L_v = 2.501e6-2370*T0; # Units : J kg-1 ** EQUATION FROM PIET ##REF##
  # Step 3e - Compute gas constant for moist air
  R_a = 287*(1+0.608*q_z); # Units : J kg-1 K-1
  # Step 3f - Compute air density
  rho_a = 100*p/(R_a*(T+273.16)); # Units : kg m-3
  # Step 3g - Compute kinematic viscosity of air
  v = (1./rho_a)*(4.94e-8*T + 1.7184e-5); # Units : m2 s-1
  # Step 3h - Compute virtual air temperature and virtual air-water temperature difference
  T_v = (T+273.16)*(1+0.61*q_z); # Units - K
  T_ov = (T0+273.16)*(1+0.61*q_s); # Units - K
  del_theta = T_ov - T_v;
  # Step 3h - Compute water density
  rho_w = 1000*(1-1.9549*0.00001*abs(T0-3.84)^1.68);
  ### End step 3

  # step 4
  u_star = U_Z *sqrt(0.00104+0.0015/(1+exp((-U_Z+12.5)/1.56))); # Amorocho and DeVries, initialise ustar using U_Z

  if (u_star == 0){
    u_star = 1e-6
  }

  z_0 = (const_Charnock*u_star^2./const_Gravity) + (0.11*v/u_star);
  z_0_prev=z_0*1.1; # To initiate the iteration
  for (i1 in 1:length(U_Z)){
    while (abs((z_0[i1] - z_0_prev[i1]))/abs(z_0_prev[i1]) > 0.000001){ # Converge when z_0 within 0.0001# of previous value
      u_star[i1]=const_vonKarman*U_Z[i1]/(log(z/z_0[i1]));  # Compute u_star
      dummy = z_0[i1]; # Used to control while loop
      z_0[i1]=(const_Charnock*u_star[i1]^2./const_Gravity) + (0.11*v[i1]/u_star[i1]); # Compute new roughness length
      z_0_prev[i1] = dummy; # Used to control while loop
    }
  }

  # Step 4d - Compute initial neutral drag coefficient
  C_DN = (u_star^2)/(U_Z^2); # Units - none
  # Step 4e - Compute roughness Reynolds number
  Re_star = u_star*z_0/v; # Units - none
  # Step 4f - Compute initial roughness length for temperature
  z_T = z_0*exp(-2.67*(Re_star)^(1/4) + 2.57); # Units - m
  z_T = Re(z_T); # Get real components, and NaN can create imag component despite no data
  # Step 4g - Compute initial roughness length for vapour
  z_E = z_0*exp(-2.67*(Re_star)^(1/4) + 2.57); # Units - m
  z_E = Re(z_E); # Get real components, and NaN can create imag component despite no data
  # Step 4h - Compute initial neutral sensible heat transfer coefficient
  C_HN = const_vonKarman*sqrt(C_DN)/(log(z/z_T));
  # Step 4i - Compute initial neutral latent heat transfer coefficient
  C_EN = const_vonKarman*sqrt(C_DN)/(log(z/z_E));
  ### End step 4

  ### Step 5 - Start iteration to compute corrections for atmospheric stability
  for (i1 in 1:length(U_Z)){
    # Step 5a - Compute initial sensible heat flux based on neutral coefficients
    H_initial = rho_a[i1]*const_SpecificHeatAir*C_HN[i1]*U_Z[i1]*(T0[i1]-T[i1]); # Units : W m-2
    # Step 5b - Compute initial latent heat flux based on neutral coefficients
    E_initial = rho_a[i1]*L_v[i1]*C_EN[i1]*U_Z[i1]*(q_s[i1]-q_z[i1]); # Units : W m-2
    # Step 5c - Compute initial Monin-Obukhov length
    L_initial = (-rho_a[i1]*u_star[i1]^3*T_v[i1])/(const_vonKarman*const_Gravity*(H_initial[i1]/const_SpecificHeatAir + 0.61*E_initial[i1]*(T[i1]+273.16)/L_v[i1])); # Units - m
    # Step 5d - Compute initial stability parameter
    zeta_initial = z/L_initial[i1];
    # Step 5e - Compute initial stability function
    psim=PSIM(zeta_initial[i1]); # Momentum stability function
    psit=PSITE(zeta_initial[i1]); # Sensible heat stability function
    psie=PSITE(zeta_initial[i1]); # Latent heat stability function
    # Step 5f - Compute corrected coefficients
    C_D=const_vonKarman*const_vonKarman/(log(z/z_0[i1])-psim[i1])^2;
    C_H=const_vonKarman*sqrt(C_D[i1])/(log(z/z_T[i1])-psit[i1]);
    C_E=const_vonKarman*sqrt(C_D[i1])/(log(z/z_E[i1])-psie[i1]);
    # Step 5g - Start iteration
    L_prev = L_initial[i1];
    L = L_prev*1.1; # Initialise while loop
    count=0;
    while (abs((L[i1] - L_prev))/abs(L_prev) > 0.000001){
      # Iteration counter
      count[i1]=count[i1]+1;
      if (count[i1] > 20){
        break;
      }
      # Step 5i - Compute new z_O, roughness length for momentum
      z_0= (const_Charnock*u_star[i1]^2./const_Gravity) + (0.11*v[i1]/u_star[i1]);
      # Step 5j - Compute new Re_star
      Re_star = u_star[i1]*z_0[i1]/v[i1];
      # Step 5k - Compute new z_T, roughness length for temperature
      z_T = z_0[i1]*exp(-2.67*(Re_star[i1])^(1/4) + 2.57);
      # Step 5l - Compute new z_E, roughness length for vapour
      z_E = z_0[i1]*exp(-2.67*(Re_star[i1])^(1/4) + 2.57);
      # Step 5p - Compute new stability parameter
      zeta = z/L[i1];
      #fprintf('zeta #g\n',zeta[i1]);
      # Step 5q - Check and enforce bounds on zeta
      if (zeta[i1] > 15){
        zeta[i1] = 15}
      else if (zeta[i1] < -15) {
        zeta[i1] = -15}
      # Step 5r - Compute new stability functions
      psim=PSIM(zeta[i1]); # Momentum stability function
      psit=PSITE(zeta[i1]); # Sensible heat stability function
      psie=PSITE(zeta[i1]); # Latent heat stability function
      # Step 5s - Compute corrected coefficients
      C_D=const_vonKarman*const_vonKarman/(log(z/z_0[i1])-psim[i1])^2;
      C_H=const_vonKarman*sqrt(C_D[i1])/(log(z/z_T[i1])-psit[i1]);
      C_E=const_vonKarman*sqrt(C_D[i1])/(log(z/z_E[i1])-psie[i1]);
      # Step 5m - Compute new H (now using corrected coefficients)
      H = rho_a[i1]*const_SpecificHeatAir*C_H[i1]*U_Z[i1]*(T0[i1]-T[i1]);
      # Step 5n - Compute new E (now using corrected coefficients)
      E = rho_a[i1]*L_v[i1]*C_E[i1]*U_Z[i1]*(q_s[i1]-q_z[i1]);
      # Step 5h - Compute new u_star
      u_star=sqrt(C_D[i1]*U_Z[i1]^2);
      # Step 5o - Compute new Monin-Obukhov length
      dummy = L[i1]; # Used to control while loop
      L = (-rho_a[i1]*u_star[i1]^3*T_v[i1])/(const_vonKarman*const_Gravity*(H[i1]/const_SpecificHeatAir + 0.61*E[i1]*(T[i1]+273.16)/L_v[i1]));
      L_prev = dummy; # Used to control while loop
    } # Converge when L within 0.0001# or previous L

  } # Need to iterate separately for each record


  ### End step 5

  # Take real values to remove any complex values that arise from missing data or NaN.
  C_D=Re(C_D);
  C_E=Re(C_E);
  C_H=Re(C_H);
  z_0=Re(z_0);
  z_E=Re(z_E);
  z_T=Re(z_T);

  # Compute evaporation [mm/day]
  Evap = 86400*1000*E/(rho_w*L_v);

  latent = E
  return( latent * (-1))
}


# https://www.r-bloggers.com/2017/08/the-trapezoidal-rule-of-numerical-integration-in-r/
composite.trapezoid <- function(f, a, b, n) {
  if (is.function(f) == FALSE) {
    stop('f must be a function with one parameter (variable)')
  }

  h <- (b - a) / n

  j <- 1:n - 1
  xj <- a + j * h
  approx <- (h / 2) * (f(a) + 2 * sum(f(xj)) + f(b))

  return(approx)
}


integrate_agg_fun <- function(dt, y, int_method){
  N = length(dt)
  if(int_method == "average"){
    out = sum((y[1:(N-1)] + y[2:N])/2 * (dt[2:N] - dt[1:(N-1)])) / (max(dt) - min(dt))
  }
  if(int_method == "integral"){
    out = sum((y[1:(N-1)] + y[2:N])/2 * (dt[2:N] - dt[1:(N-1)]))
  }
  return(out)
}

get_interp_drivers <- function(meteo_all, total_runtime, hydrodynamic_timestep = 86400, dt, method="interp", int_method="average",
                               secchi = T){
  times = seq(1, 1 + total_runtime*hydrodynamic_timestep, dt)
  meteo = matrix(NA, nrow = 9, ncol = length(times))
  if(method == "interp"){
    meteo[1,] = approx(x = meteo_all[[1]]$dt,
                       y = meteo_all[[1]]$Shortwave_Radiation_Downwelling_wattPerMeterSquared,
                       xout = times,
                       method = "linear", rule = 2)$y # 1 = Jsw
    meteo[2,] = approx(x = meteo_all[[1]]$dt,
                       y = meteo_all[[1]]$Longwave_Radiation_Downwelling_wattPerMeterSquared,
                       xout = times,
                       method = "linear", rule = 2)$y # 2 = Jlw
    meteo[3,] = approx(x = meteo_all[[1]]$dt,
                       y = meteo_all[[1]]$Air_Temperature_celsius,
                       xout = times,
                       method = "linear", rule = 2)$y # 3 = Tair
    meteo[4,] = approx(x = meteo_all[[1]]$dt,
                       y = meteo_all[[1]]$ea,
                       xout = times,
                       method = "linear", rule = 2)$y # 4 = ea
    meteo[5,] = approx(x = meteo_all[[1]]$dt,
                       y = meteo_all[[1]]$Ten_Meter_Elevation_Wind_Speed_meterPerSecond,
                       xout = times,
                       method = "linear", rule = 2)$y # 5 = Uw
    meteo[6,] = approx(x = meteo_all[[1]]$dt,
                       y = meteo_all[[1]]$Cloud_Cover,
                       xout = times,
                       method = "linear", rule = 2)$y # 6 = CC
    meteo[7,] = approx(x = meteo_all[[1]]$dt,
                       y = meteo_all[[1]]$Surface_Level_Barometric_Pressure_pascal,
                       xout = times,
                       method = "linear", rule = 2)$y # 7 = Pa
    meteo[8,] = approx(x = meteo_all[[2]]$dt,
                       y = meteo_all[[2]]$kd,
                       xout = times,
                       method = "linear", rule = 2)$y # 8 = kd
    meteo[9,] = approx(x = meteo_all[[1]]$dt,
                       y = meteo_all[[1]]$Relative_Humidity_percent,
                       xout = times,
                       method = "linear", rule = 2)$y # 9 = RH
  }
  if(method == "integrate"){
    meteo_all[[1]]$dt = as.numeric(meteo_all[[1]]$dt)
    # get times at exact dt intervals
    x_dt = data.frame(dt = times) %>% left_join(meteo_all[[1]])
    # duplicate above and add "add_to_group" so averages/integrals are calculated from endpoint to endpoint
    x_dt_2 = bind_rows(x_dt %>% mutate(add_group = -1), x_dt %>% mutate(add_group = 0))
    # get any measurements that weren't at dt intervals
    measurements = meteo_all[[1]] %>%
      dplyr::filter(!(dt %in% x_dt$dt))
    # join and sort above
    comb = full_join(x_dt_2, measurements %>% mutate(dt = as.numeric(dt))) %>%
      arrange(dt, add_group) %>%
      dplyr::filter(dt <= max(times))
    # linearly interpolate to present dt's so missing values are filled
    cols_interp_met = c("Shortwave_Radiation_Downwelling_wattPerMeterSquared", "Longwave_Radiation_Downwelling_wattPerMeterSquared", "Air_Temperature_celsius", "ea", "Ten_Meter_Elevation_Wind_Speed_meterPerSecond", "Cloud_Cover", "Surface_Level_Barometric_Pressure_pascal", "Relative_Humidity_percent")
    for(i in 1:length(cols_interp_met)){
      comb[, cols_interp_met[i]] = approx(comb$dt, comb[, cols_interp_met[i]], comb$dt, method="linear", rule=2)$y
    }
    if (secchi){
      comb[, "kd"] = approx(meteo_all[[2]]$dt, meteo_all[[2]]$kd, comb$dt, method="linear", rule=2)$y
    } else {
      comb[, "kd"] = -999
    }
    # add group column
    dt_hold = dt
    comb = comb %>%
      mutate(group = dt %/% dt_hold) %>%
      mutate(group = ifelse(!is.na(add_group), group + add_group, group)) %>%
      dplyr::filter(group >=0 )
    # aggregate to group
    integral = comb %>%
      arrange(group, dt) %>%
      group_by(group) %>%
      summarise(across(all_of(c(cols_interp_met, "kd")), ~ integrate_agg_fun(dt, ., int_method)))
    # format into matrix to be returned
    cols_interp_ordered = c("Shortwave_Radiation_Downwelling_wattPerMeterSquared", "Longwave_Radiation_Downwelling_wattPerMeterSquared", "Air_Temperature_celsius", "ea", "Ten_Meter_Elevation_Wind_Speed_meterPerSecond", "Cloud_Cover", "Surface_Level_Barometric_Pressure_pascal", "kd", "Relative_Humidity_percent")
    meteo = integral %>%
      select(cols_interp_ordered) %>%
      as.matrix() %>%
      t()
  }

  rownames(meteo) = c("Jsw", "Jlw", "Tair", "ea", "Uw", "CC", "Pa", "kd", "RH")

  return(meteo)
}

run_thermalmodel <- function(u, # initial temperature profile
                             startTime, # starting time
                             endTime, # ending time
                             ice = FALSE, # is ice exisitng?
                             Hi = 0, # ice thickness
                             iceT = 6, # movign average for ice formation
                             supercooled = 0, # how many layers below 0 deg C?
                             scheme = 'implicit', # numerical diffusion scheme
                             kd_light = NULL, # light extinction
                             densThresh = 1e-3, # threshold for convective overturn
                             reflect = 0.3, # light reflection
                             infra = 0.7, # infrared absorption
                             albedo = 0.1, # water albedo
                             eps = 0.97, # epsilon
                             emissivity = 0.97, # water emissivity
                             sigma = 5.67 * 10^(-8), # constant
                             p2 = 1, # calibraiton value
                             B = 0.61,
                             g = 9.81, # gravity
                             Cd = 0.0013, # momentum coefficient (wind)
                             sw_factor = 1.0, # multiplier for short-wave radiation
                             meltP = 1, # multiplier for melting energy
                             dt_iceon_avg = 0.8, # determines ice formation
                             Hgeo = 0.1, # geothermal heat flux
                             KEice = 1/1000, # wind attenuation for ice cover
                             Ice_min = 0.1, # minimum ice thickness
                             area, # area
                             depth, # depth
                             volume, # volume
                             zmax, # maximum lake depth
                             nx, # number of layers we will have
                             dt, # time step
                             dx, # space step
                             daily_meteo # meteorology
                             ){

  print('What a beautiful day to run a lake model.')

  um <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
  kzm <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
  n2m <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
  mix <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))#(floor(endTime/dt - startTime/dt)))
  therm.z <- rep(NA, length =length(seq(startTime, endTime, dt)/dt))
  mix.z <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))
  Him <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))

  SW <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))
  LW_in <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))
  LW_out <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))
  LAT <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))
  SEN <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))

  pb = txtProgressBar(min = startTime, max = endTime/dt, initial = 0)
  stepi = 1
  start.time <- Sys.time()

  ## modeling code for vertical 1D mixing and heat transport
  for (n in seq(startTime, endTime/dt, 1)){

    if (!is.null(kd_light)){
      kd = kd_light
    }else{
      kd = daily_meteo["kd",n]
    }

    un = u # prior temperature values

    kz = eddy_diffusivity(calc_dens(un), depth, 9.81, 998.2, ice, area) / 86400

    if (ice & daily_meteo["Tair",n] <= 0){
      kzn = kz
      absorp = 1 - 0.7
      infra = 1 - absorp
      albedo = 0.7
    } else if (ice & daily_meteo["Tair",n] >= 0){
      kzn = kz
      absorp = 1 - 0.3
      infra = 1 - absorp
      albedo = 0.7
    } else if (!ice) {
      kzn = kz
      absorp = 1 - reflect# 0.3
      infra = 1 - absorp
      albedo = 0.1
    }
    kzm[, n] <- kzn


    ## (1) Heat addition
    # surface heat flux
    Q <- (
      longwave(cc = daily_meteo["CC",n], sigma = sigma, Tair = daily_meteo["Tair",n], ea = daily_meteo["ea",n], emissivity = emissivity, Jlw = daily_meteo["Jlw",n]) +
        backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps) +
        latent(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd) +
        sensible(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd))
    SW[n] <- absorp * daily_meteo["Jsw",n]
    LW_in[n] <- longwave(cc = daily_meteo["CC",n], sigma = sigma, Tair = daily_meteo["Tair",n], ea = daily_meteo["ea",n], emissivity = emissivity, Jlw = daily_meteo["Jlw",n])
    LW_out[n] <- backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps)
    LAT[n] <- latent(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd)
    SEN[n] <- sensible(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd)


    # heat addition over depth
    H =  (1- albedo) * (daily_meteo["Jsw",n] * sw_factor) * exp(-(kd) *depth)

    Hg <- (area-lead(area))/dx * Hgeo/(4181 * calc_dens(un))
    Hg[nx] <- (area[nx-1]-area[nx])/dx * Hgeo/(4181 * calc_dens(un[nx]))



    # add heat to all layers
    ## (2) DIFFUSION
    if (scheme == 'implicit'){
      ## (2a) Boundary heat addition
      # surface layer
      u[1] = un[1] +
        (Q * area[1]/(dx)*1/(4184 * calc_dens(un[1]) ) +
           abs(H[1+1]-H[1]) * area[1]/(dx) * 1/(4184 * calc_dens(un[1]) ) +
           Hg[1]) * dt/area[1]


      # all other layers in between
      for (i in 2:(nx-1)){
        u[i] = un[i] +
          (abs(H[i+1]-H[i]) * area[i]/(dx) * 1/(4184 * calc_dens(un[i]) ) +
             Hg[i])* dt/area[i]

      }

      # bottom layer
      u[nx] = un[nx] +
        (abs(H[nx]-H[nx-1]) * area[nx]/(area[nx]*dx) * 1/(4181 * calc_dens(un[nx])) +
           Hg[nx]/area[nx]) * dt

      ## (2b) Diffusion by Crank-Nicholson Scheme (CNS)
      j <- length(u)
      y <- array(0, c(j,j))

      # all other layers in between
      # Linearized heat conservation equation matrix (diffusion only)
      alpha = (dt/dx**2) * kzn
      az <- -alpha                               #coefficient for i-1
      bz <- 2 * (1 + alpha)                                                       #coefficient for i
      cz <- -alpha                #coefficient for i+1

      #Boundary conditions, surface
      az[1] <- 0
      bz[1]<- 1
      cz[length(cz)] <- 0
      bz[length(bz)] <- 1


      y[0 + 1:(j - 1) * (j + 1)] <- cz[-length(bz)]	# superdiagonal
      y[1 + 0:(j - 1) * (j + 1)] <- bz	# diagonal
      y[2 + 0:(j - 2) * (j + 1)] <- az[-1] 	# subdiagonal

      y[1,2] <- 0#
      y[nrow(y), (ncol(y)-1)] = 0

      mn <- rep(0, j)
      mn[1] = u[1]
      mn[j] = u[j]
      for (g in 2:(j-1)){
        mn[g] = alpha[g] * u[g-1] + 2 * (1-alpha[g])*u[g] + alpha[g] * u[g+1]
      }

      u  <- solve(y, mn)

    }


    # surface layer
    if (scheme == 'explicit'){ # forward time centered space (FTCS)
      u[1] = un[1] +
        (Q * area[1]/(dx)*1/(4184 * calc_dens(un[1]) ) +
           abs(H[1+1]-H[1]) * area[1]/(dx) * 1/(4184 * calc_dens(un[1]) ) +
           Hg[1]) * dt/area[1]

      # all other layers in between
      for (i in 2:(nx-1)){
        u[i] = un[i] +
          (area[i] * kzn[i] * 1 / dx**2 * (un[i+1] - 2 * un[i] + un[i-1]) +
             abs(H[i+1]-H[i]) * area[i]/(dx) * 1/(4184 * calc_dens(un[i]) ) +
             Hg[i])* dt/area[i]
      }
      # bottom layer
      u[nx] = un[nx] +
        abs(H[nx]-H[nx-1]) * area[nx]/(area[nx]*dx) * 1/(4181 * calc_dens(un[nx]) +
                                                           Hg[nx]/area[nx]) * dt
    }



    ## (3) TURBULENT MIXING OF MIXED LAYER
    # the mixed layer depth is determined for each time step by comparing kinetic
    # energy available from wind and the potential energy required to completely
    # mix the water column to a given depth

    Zcv <- depth %*% area / sum(area) # center of volume
    tau = 1.225 * Cd * daily_meteo["Uw",n]^2 # wind shear is air density times wind velocity
    if (daily_meteo["Uw",n] <= 15) {
      c10 = 0.0005 * sqrt(daily_meteo["Uw",n])
    } else {
      c10 = 0.0026
    }
    shear = sqrt((c10 * calc_dens(un[1]))/1.225) *  daily_meteo["Uw",n] # shear velocity
    # coefficient times wind velocity squared
    KE = shear *  tau * dt # kinetic energy as function of wind

    if (ice){
      KE = KE * KEice
    }
    maxdep = 1
    for (dep in 1:(nx-1)){
      if (dep == 1){

        PE = abs(g *   depth[dep] *( depth[dep+1] - Zcv)  *

                   abs(calc_dens(u[dep+1])- mean(calc_dens(u[1:dep]))))
      } else {
        PEprior = PE

        PE = abs(g *   depth[dep] *( depth[dep+1] - Zcv)  *

                   abs(calc_dens(u[dep+1])- mean(calc_dens(u[1:dep])))) + PEprior

      }
      if (PE > KE){
        maxdep = dep-1
        break
      } else if (dep>1 & PE < KE ){
        u[(dep-1):dep] = (u[(dep-1):dep] %*% volume[(dep-1):dep])/sum(volume[(dep-1):dep])
      }
      maxdep = dep
    }

    mix[n] <- KE/PE
    therm.z[n] <- maxdep


    ## (4) DENSITY INSTABILITIES
    # convective overturn: Convective mixing is induced by an unstable density
    # profile. All groups of water layers where the vertical density profile is
    # unstable are mixed with the first stable layer below the unstable layer(s)
    # (i.e., a layer volume weighed means of temperature and other variables are
    # calculated for the mixed water column). This procedure is continued until
    # the vertical density profile in the whole water column becomes neutral or stable.

    dens_u = calc_dens(u)
    diff_dens_u <- (diff(dens_u))
    diff_dens_u[abs(diff(dens_u)) <= densThresh] = 0
    while (any(diff_dens_u < 0)){
      dens_u = calc_dens(u)
      for (dep in 1:(nx-1)){
        if (dens_u[dep+1] < dens_u[dep] & abs(dens_u[dep+1] - dens_u[dep]) >= densThresh){
          u[dep:(dep+1)] = (u[dep:(dep+1)] %*% volume[dep:(dep+1)])/sum(volume[dep:(dep+1)]) #mean(u[dep:(dep+1)])
          break
        }
      }
      dens_u = calc_dens(u)
      diff_dens_u <- (diff(dens_u))
      diff_dens_u[abs(diff(dens_u)) <= densThresh] = 0
    }

    dens_u_n2 = calc_dens(u)
    n2 <- 9.81/mean(calc_dens(u)) * (lead(dens_u_n2) - dens_u_n2)/dx
    max.n2 <- ifelse(max(n2, na.rm = T) > 1E-4, which.max(n2) * dx, dx * nx)
    mix.z[n] <- max.n2



    ## (5) ICE FORMATION
    # according to Hostetler & Bartlein (1990):
    # (1) ice forms when surface water temp <= 1 deg C and melts when > 1 deg
    # (2) rate of ice formation/melting is exponential function of ice thickness
    # (the thicker the ice, the slower the formation rate, and vice versa)
    # (3) heat of fusion is added/subtracted from surface energy balance
    # (4) diffusion below ice only happens on molecular level
    # (5) with ice, surface absorption of incoming solar radiation increases to 85 %

    icep  = max(dt_iceon_avg,  (dt/86400))
    x = (dt/86400) / icep
    iceT = iceT * (1 - x) + u[1] * x

    Him[ n] <- Hi


    if ((iceT <= 0) == TRUE & Hi < Ice_min & daily_meteo["Tair",n] <= 0){
      supercooled <- which(u < 0)
      initEnergy <- sum((0-u[supercooled])*area[supercooled] * dx * 4.18E6)

      if (ice != TRUE) {
        Hi <- Ice_min+(initEnergy/(910*333500))/max(area)
      } else {
        if (daily_meteo["Tair",n] > 0){
          Tice <- 0 #
          Hi = Hi -max(c(0, meltP * dt*((absorp*daily_meteo["Jsw",n])+(longwave(cc = daily_meteo["CC",n], sigma = sigma, Tair = daily_meteo["Tair",n], ea = daily_meteo["ea",n], emissivity = emissivity, Jlw = daily_meteo["Jlw",n]) +
                                                                         backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps) +
                                                                         latent(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n],  RH = daily_meteo["RH",n], A = area, Cd = Cd) +
                                                                         sensible(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd)) )/(1000*333500)))
        } else {
          Tice <-  ((1/(10 * Hi)) * 0 + daily_meteo["Tair",n]) / (1 + (1/(10 * Hi)))
          Hi <- max(Ice_min, sqrt(Hi**2 + 2 * 2.1/(910 * 333500)* (0 - Tice) * dt))
        }
      }
      ice = TRUE
      if (Hi >= 0){
        u[supercooled] = 0
        u[1] = 0
      }
      Him[ n] <- Hi
    } else if (ice == TRUE & Hi >= Ice_min) {
      if (daily_meteo["Tair",n] > 0){
        Tice <- 0 #
        Hi = Hi -max(c(0, meltP * dt*((absorp*daily_meteo["Jsw",n])+(longwave(cc = daily_meteo["CC",n], sigma = sigma, Tair = daily_meteo["Tair",n], ea = daily_meteo["ea",n], emissivity = emissivity, Jlw = daily_meteo["Jlw",n]) +
                                                                        backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps) +
                                                                       latent(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n],  RH = daily_meteo["RH",n], A = area, Cd = Cd) +
                                                                       sensible(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd)) )/(1000*333500)))
      } else {
        Tice <-  ((1/(10 * Hi)) * 0 +  daily_meteo["Tair",n]) / (1 + (1/(10 * Hi)))
        Hi <- max(Ice_min, sqrt(Hi**2 + 2 * 2.1/(910 * 333500)* (0 - Tice) * dt))
      }
      u[supercooled] = 0
      u[1] = 0
      Him[ n] <- Hi
    } else if (ice == TRUE & Hi < Ice_min){
      ice = FALSE
      Him[ n] <- Hi
    }

    if (ice == FALSE){
      Hi = 0
      Him[n] = Hi
    }

    n2m[, n] <- n2
    um[, n] <- u


    stepi = stepi + 1
    setTxtProgressBar(pb,stepi)
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)

  df.sim <- data.frame(cbind(seq(startTime,endTime,dt), t(um)) )
  colnames(df.sim) <- c("datetime", as.character(paste0('wtemp.',seq(1,nrow(um))*dx)))
  df.sim$datetime <-   seq(startTime, endTime, dt)#/24/3600

  ## averaged responses
  bf.sim <- apply(df.sim[,-1], 1, function(x) rLakeAnalyzer::buoyancy.freq(wtr = x, depths = seq(1,nrow(um))*dx))

  bf.sim <- apply(df.sim[,-1], 1, function(x) rLakeAnalyzer::center.buoyancy(wtr = x, depths = seq(1,nrow(um))*dx))

  df.z.df.sim <- data.frame('time' = df.sim$datetime, 'z' = bf.sim)

  avg.epi.sim <- NULL
  avg.hyp.sim <- NULL
  avg.tot.sim <- NULL
  for (j in 1:nrow(df.z.df.sim)){
    d = df.sim[,-1]
    if (is.na(df.z.df.sim$z[j])){
      df.z.df.sim$z[j] = 1
    }
    avg.epi.sim <- append(avg.epi.sim,((as.numeric(d[j,1:df.z.df.sim$z[j]], na.rm = T) %*%
                                          area[1:df.z.df.sim$z[j]] )/
                                         sum(area[1:df.z.df.sim$z[j]])))
    avg.hyp.sim <- append(avg.hyp.sim,((as.numeric(d[j,df.z.df.sim$z[j]:ncol(d)], na.rm = T)%*%
                                          area[df.z.df.sim$z[j]:ncol(d)] )/
                                         sum(area[df.z.df.sim$z[j]:ncol(d)])))
    avg.tot.sim <- append(avg.tot.sim,((as.numeric(d[j,1:ncol(d)], na.rm = T)%*%
                                          area[1:ncol(d)] )/
                                         sum(area[1:ncol(d)])))
  }

  stratFlag = rep(NA, length = ncol(um))
  for (v in 1:length(stratFlag)){
    stratFlag[v] = ifelse((calc_dens(um[nx,v]) - calc_dens(um[1,v])) >= 0.1 &
                            mean(um[,v]) >= 4, 1, 0)
  }


  df.avg.sim <- data.frame('time' = df.sim$datetime,
                           'epi' = avg.epi.sim,
                           'hyp' = avg.hyp.sim,
                           'tot' = avg.tot.sim,
                           'stratFlag' = stratFlag,
                           'thermoclineDep' = bf.sim)
  print("Have a lovely rest of your day!")
  dat = list('temp'  = um,
             'diff' = kzm,
             'mixing' = mix,
             'buoyancy' = n2m,
             'icethickness' = Him,
             'iceflag' = ice,
             'icemovAvg' = iceT,
             'supercooled' = supercooled,
             'mixingdepth' = mix.z,
             'thermoclinedepth' = therm.z,
             'endtime' = endTime,
             'average' = df.avg.sim,
             'SW' = SW,
             'LW_in' = LW_in,
             'LW_out' = LW_out,
             'LAT' = LAT,
             'SEN' = SEN)

  return(dat)


}

