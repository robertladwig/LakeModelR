#' Runs 1D integral energy model
#'
#' Returns a data frame containing model output in the form of vertical temperature
#' profiles, diffusivity values, mixing enery, buoyancy frequencies, ice thickness,
#' logical value for ice, ice moving average temperature, depth of supercooled
#' layers, mixing depth, thermocline depth, end time of model run, averaged
#' temperatures for epilimnion and hypolimnion volumes, incoming short-wave
#' radiation, incoming long-wave radiation, backscattered long-wave radiation,
#' latent and sensible heat fluxes.
#'
#' The mixing algorithms are based on the MINLAKE (Ford and Stefan 1980, Riley
#' and Stefan 1988, Herb and Stefan 2004) and the MyLake (Saloranta and Andersen 2007)
#' models. Implementations to estimate the incoming and outgoing long-wave heat fluxes
#' were taken from Livingstone and Imboden (1989) and Goudsmit et al. (2002). The latent
#' and sensible heat fluxes were calculated taking into account atmospheric stability using
#' the algorithms by Verburg and Antenucci (2010). The ice algorithms from MyLake (Saloranta
#' and Andersen 2007) were applied to simulate ice formation and melting.
#'
#' References:
#' Ford, D.E., and H.G. Stefan. 1980. Thermal predictions using an integral energy model. J. Hydraul. Div. ASCE 106(1). 39-55
#'
#' Goudsmit, G.H., H. Burchard, F. Peeters, and A. Wüst. 2002. Application of k-e turbulence models to enclosed basins: The role of internal seiches. Journal of Geophysical Research 107. C12. 3230. doi:10.1029/2001JC000954
#'
#' Herb, W.R., and H.G. Stefan. 2004. Temperature Stratification and Mixing Dynamics in a Shallow Lake With Submersed Macrophytes. Lake and Reservoir Management 20. 4. 296-308. doi:10.1080/07438140409354159
#'
#' Livingstone, D., and D. Imboden. 1989. Annual heat balance and equilibrium temperature of Lake Aegeri, Switzerland. Aquat. Sci. 51
#'
#' Riley, M., and H.G. Stefan. 1988. MINLAKE: A dynamic lake water quality simulation model. Ecol. Model. 43. 155-182
#'
#' Saloranta, T.M., and T. Andersen. 2007. MyLake – A multi-year lake simulation model code suitable for uncertainty and sensitivity analysis simulation. Ecol. Model. 207. 1. 45-60. doi:https://doi.org/10.1016/j.ecolmodel.2007.03.018
#'
#' Verburg, P., and J.P. Antenucci. 2010. Persistent unstable atmospheric boundary layer enhances sensible and latent heat loss in a tropical great lake: Lake Tanganyika. Journal of Geophysical Research Atmospheres 115. D11. doi:https://doi.org/10.1029/2009JD012839
#'
#' @param u vector; initial temperature profile
#' @param startTime integer; starting time inseconds of the meteorology driver data
#' @param endTime integer; starting time inseconds of the meteorology driver data
#' @param ice logical; determines if ice is existing or not at model start. Defaults to FALSE
#' @param Hi double; ice thickness (m) at model start. Defaults to 0
#' @param iceT double; moving average of air temperature for ice formation. Defaults to 6 deg C
#' @param supercooled integer; amount of layers below 0 deg C at model start. Defaults to 0
#' @param scheme character; numerical diffusion scheme, either 'explicit' or 'implicit'
#' @param kd_light double; if NULL then the dataframe from daily_meteo is used, otherwise a DOUBLE value for light extinction (m-1)
#' @param densThresh double; density difference threshold between neighboring layers for convective overturn. Defaults to 1E-3
#' @param albedo double; water albedo. Defaults to 0.1
#' @param eps double; epsilon. Defaults to 0.97
#' @param emissivity double; water emissivity. Defaults to 0.97
#' @param sigma double; constant sigma value. Defaults to 5.67 * 10^(-8)
#' @param p2 double; calibration multiplier. Defaults to 1.0
#' @param B double; constant. Defaults to 0.61
#' @param g double; gravity. Defaults to 9.81 m/s2
#' @param Cd double; momentum drag coefficient. Defaults to 0.0013
#' @param sw_factor double; calibration multiplier for short-wave radiation. Defaults to 1.0
#' @param meltP double; calibration multiplier for ice melting. Defaults to 1.0
#' @param dt_iceon_avg double; determines weighting of ice onset. Defaults to 0.8
#' @param Hgeo double; geothermal heat flux. Defaults to 0.1 W/m2
#' @param KEice double; TKE multiplier due to ice cover. Defaults to 0
#' @param Ice_min double; minimum ice thickness in m. Defaults to 0.1
#' @param area vector; area over depth in m2
#' @param depth vector; depths in m references from surface
#' @param volume vector; volumes over depth in m3
#' @param zmax double; maximum lake depth in m
#' @param nx double; maximum number of grid cells
#' @param dt double; model time step in seconds
#' @param dx double; model spatial discretisation in m
#' @param daily_meteo dataframe; meteorological driver data
#' @param Cdplant double; drag coefficient macrophytes. Defaults to 0.0
#' @param ahat double; surface to volume ratio macrophytes. Defaults to 0.4
#' @param km double; ight extinction macrophytes in m-1. Defaults to 0.0
#' @param rho_plant double; density macrophytes. Defaults to 30 kg/m3
#' @param canpy dataframe; canophy height time series, NULL refers to no macrophyte simulations. Defaults to NULL
#' @param biomass dataframe; biomass time series, NULL refers to no macrophyte simulations. Defaults to NULL
#' @param water.quality logical; if TRUE runs water quality (need a defined function 'water_quality_boundary_conditions'). Defaults to FALSE
#' @param wq vector; initial vector of water quality variable profile, needs water.quality = TRUE. Defaults to NULL
#' @param agents vector; initial vector of individual phytoplankton cells. Defaults to NULL
#' @param Hs double; snow thickness (m) at model start. Defaults to 0
#' @param rho_snow double; snow density (kg/m3). Defaults to 250
#' @param Hsi double; snow-ice thickness (m) at model starts. Default to 0
#' @param rho_ice double; ice density (kg/m3). Defaults to 910
#' @param rho_fw double; freshwater density (kg/m3). Defaults to 1000
#' @param rho_new_snow double; new snow density (kg/m3). Defaults to 250
#' @param rho_max_snow double; maximum snow density (kg/m3). Defaults to 450
#' @param K_ice double; thermal conductivity of ice (W/K/m). Defaults to 2.1
#' @param Cw double; volumetric heat capacity of water (J/K/m3). Defaults to 4.18E6
#' @param L_ice double; latent heat of freezing (J/kg). Defaults to 333500
#' @param kd_snow double; snow albedo. Default to 0.9
#' @param kd_ice double; ice albedoa. Defaults to 0.7
#' @param pgdl_model logical; get additional diagnostics. Defaults to FALSE
#' @author Robert Ladwig
#'
#' @examples
#' \dontrun{
#' res <-  run_thermalmodel(u = u_ini,
#' startTime = startTime,
#' endTime =  endTime,
#' ice = FALSE,
#' Hi = 0,
#' iceT = 6,
#' supercooled = 0,
#' kd_light = 0.5,
#' sw_factor = 1.0,
#' zmax = zmax,
#' nx = nx,
#' dt = dt,
#' dx = dx,
#' area = hyps_all[[1]], # area
#' depth = hyps_all[[2]], # depth
#' volume = hyps_all[[3]], # volume
#' daily_meteo = meteo,
#' Cd = 0.0013,
#' scheme = 'implicit')
#' }
#'
#' @export
#'
run_thermalmodel <- function(u, # initial temperature profile
                             startTime, # starting time
                             endTime, # ending time
                             ice = FALSE, # is ice exisitng?
                             Hi = 0, # ice thickness
                             iceT = 6, # moving average for ice formation
                             supercooled = 0, # how many layers below 0 deg C?
                             scheme = 'implicit', # numerical diffusion scheme
                             kd_light = NULL, # light extinction
                             densThresh = 1e-3, # threshold for convective overturn
                             albedo = 0.1, # water albedo
                             eps = 0.97, # epsilon
                             emissivity = 0.97, # water emissivity
                             sigma = 5.67 * 10^(-8), # constant
                             p2 = 1.0, # calibration value
                             B = 0.61,
                             g = 9.81, # gravity
                             Cd = 0.0013, # momentum coefficient (wind)
                             sw_factor = 1.0, # multiplier for short-wave radiation
                             meltP = 1.0, # multiplier for melting energy
                             dt_iceon_avg = 0.8, # determines ice formation
                             Hgeo = 0.1, # geothermal heat flux
                             KEice = 0, # wind attenuation for ice cover
                             Ice_min = 0.1, # minimum ice thickness
                             area, # area
                             depth, # depth
                             volume, # volume
                             zmax, # maximum lake depth
                             nx, # number of layers we will have
                             dt, # time step
                             dx, # space step
                             daily_meteo, # meteorology
                             Cdplant = 0.0, # drag coefficient macrophytes
                             ahat = 0.4, # surface to volume macrophytes
                             km = 0.0, # light extinction macrophytes
                             rho_plant = 30, # density macrophytes
                             canpy = NULL, # canophy height time series
                             biomass = NULL, # biomass time series
                             water.quality = FALSE, # check if water quality simulation should be run
                             wq = NULL, # initial water quality profile
                             agents = NULL, # initial profile of phytoplankton individuals
                             Hs = 0,
                             rho_snow = 250,
                             Hsi = 0,
                             rho_ice = 910,
                             rho_fw = 1000,
                             rho_new_snow = 250,
                             rho_max_snow = 450,
                             K_ice = 2.1,
                             Cw = 4.18E6,
                             L_ice = 333500,
                             kd_snow = 0.9,
                             kd_ice = 0.7,
                             pgdl_mode = FALSE
){

  greet <- data.frame(greet = c('What a beautiful day to run a lake model.',
                                'What is the deepest lake in the world?',
                                'What is the highest elevation lake in the world?',
                                'Which country has the most lakes in the world?'),
                      bye = c('Have a lovely rest of your day!',
                              'Lake Baikal is the deepest lake of the world (1,620 meters [5,315 feet])',
                              'Ojos del Salado is the highest active volcano and fresh waterbody of the world, at 6,390 meters (20,965 feet).',
                              'Canada has the most lakes in the world: an estimated 879,800. Russia comes second with about 201,200 lakes.'))

  whichgreet <- sample(x = 1:nrow(greet), size = 1) # Sample one greeting message

  cat(greet[whichgreet, 1])

  if (is.null(canpy)){
    canpy = data.frame('dt' = c(0, endTime/dt), 'mean_canopy' = c(0,0))
  }
  if (is.null(biomass)){
    biomass = data.frame('dt' = c(0, endTime/dt), 'mean_biomass' = c(0,0))
  }
  macroheight <- approxfun(x = canpy$dt, y = canpy$mean_canopy, method = "linear", rule = 2)
  macrobiomss <- approxfun(x = biomass$dt, y = biomass$mean_biomass, method = "linear", rule = 2)

  um <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
  wqm <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
  kzm <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
  n2m <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
  mix <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))#(floor(endTime/dt - startTime/dt)))
  therm.z <- rep(NA, length =length(seq(startTime, endTime, dt)/dt))
  mix.z <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))
  Him <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))
  Hsm <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))
  Hsim <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))
  Ticem <- rep(NA, length = length(seq(startTime, endTime, dt)/dt))
  magents <- matrix(NA, ncol  = length(seq(startTime, endTime, dt)/dt), nrow = length(agents))
  mloc <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)

  if (pgdl_mode){
    um_initial <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
    um_heat <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
    um_diff <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
    um_mix <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
    um_conv <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
    um_ice <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
    n2_pgdl <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = nx)
    meteo_pgdl <- matrix(NA, ncol = length(seq(startTime, endTime, dt)/dt), nrow = 9)
  }

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

    kz = eddy_diffusivity(calc_dens(un), depth, 9.81,
                          mean(calc_dens(un)), ice, area) / 86400
    if (pgdl_mode){
      rho <- calc_dens(un)
      buoy = rep(1, (nx)) * 7e-5
      buoy[1:(nx-1)] = abs(rho[2:nx] - rho[1:(nx-1)]) / (depth[2:nx] - depth[1:(nx-1)]) * g/mean(rho)
      buoy[nx] = ( abs(rho[nx-1] - rho[nx]) / abs(depth[nx-1] - depth[nx]) *
                     g/mean(rho) )
      n2_pgdl[, n] = buoy
    }

    if (ice & daily_meteo["Tair",n] <= 0){
      kzn = kz
      albedo = 0.3
      IceSnowAttCoeff <- exp(-kd_ice * Hi) * exp(-kd_snow * (rho_fw/rho_snow)* Hs)
    } else if (ice & daily_meteo["Tair",n] >= 0){
      kzn = kz
      albedo = 0.3
      IceSnowAttCoeff <- exp(-kd_ice * Hi) * exp(-kd_snow * (rho_fw/rho_snow)* Hs)
    } else if (!ice) {
      kzn = kz
      albedo = 0.1
      IceSnowAttCoeff <- 1
    }
    kzm[, n] <- kzn


    ## (1) Heat addition
    # surface heat flux
    if (!ice){
    Q <- (
      longwave(cc = daily_meteo["CC",n], sigma = sigma, Tair = daily_meteo["Tair",n], ea = daily_meteo["ea",n], emissivity = emissivity, Jlw = daily_meteo["Jlw",n]) +
        backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps) +
        latent(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd) +
        sensible(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd))
    } else {
      Q <- 0
    }
    SW[n] <- (1-albedo) * daily_meteo["Jsw",n]
    LW_in[n] <- longwave(cc = daily_meteo["CC",n], sigma = sigma, Tair = daily_meteo["Tair",n], ea = daily_meteo["ea",n], emissivity = emissivity, Jlw = daily_meteo["Jlw",n])
    LW_out[n] <- backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps)
    LAT[n] <- latent(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd)
    SEN[n] <- sensible(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd)


    Hmacrophytes <- seq(dx, zmax, length.out = nx)
    Hmacrophytes <- ifelse(Hmacrophytes <= ((max(Hmacrophytes) - macroheight(n) )), 0, 1)

    P = macrobiomss(n)

    # heat addition over depth
    if (ice){
      H =  (1- albedo) * (daily_meteo["Jsw",n] * sw_factor) * exp(-(kd + km * P * Hmacrophytes) *depth)
    } else {
      H =  IceSnowAttCoeff * (daily_meteo["Jsw",n] * sw_factor) * exp(-(kd + km * P * Hmacrophytes) *depth)
    }


    Hg <- (area-lead(area))/dx * Hgeo/(4181 * calc_dens(un))
    Hg[nx] <- (area[nx-1]-area[nx])/dx * Hgeo/(4181 * calc_dens(un[nx]))


    if (pgdl_mode){
      proxy_u <- un

      um_initial[, n] <- proxy_u
      # surface layer
      proxy_u[1] = un[1] +
        (Q * area[1]/(dx)*1/(4184 * calc_dens(un[1]) ) +
           abs(H[1+1]-H[1]) * area[1]/(dx) * 1/(4184 * calc_dens(un[1]) ) +
           Hg[1]) * dt/area[1]


      # all other layers in between
      for (i in 2:(nx-1)){
        proxy_u[i] = un[i] +
          (abs(H[i+1]-H[i]) * area[i]/(dx) * 1/(4184 * calc_dens(un[i]) ) +
             Hg[i])* dt/area[i]

      }

      # bottom layer
      proxy_u[nx] = un[nx] +
        (abs(H[nx]-H[nx-1]) * area[nx]/(area[nx]*dx) * 1/(4181 * calc_dens(un[nx])) +
           Hg[nx]/area[nx]) * dt

      um_heat[, n] <- proxy_u
    }

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

    if (pgdl_mode){
      um_diff[, n] <- u
    }

    if (water.quality){

      wq <- water_quality_boundary_conditions(WQ = wq, TEMP = u, WIND = daily_meteo["Uw", n],
                                              AREA = area, VOLUME = volume, ICE = ice, dt = dt,
                                              dx = dx, nx = nx)

      j <- length(wq)
      y <- array(0, c(j,j))

      alpha = (dt/dx**2) * kzn
      az <- -alpha
      bz <- 2 * (1 + alpha)
      cz <- -alpha

      az[1] <- 0

      bz[1]<- 1

      cz[length(cz)] <- 0
      bz[length(bz)] <- 1
      y[0 + 1:(j - 1) * (j + 1)] <- cz[-length(bz)]
      y[1 + 0:(j - 1) * (j + 1)] <- bz
      y[2 + 0:(j - 2) * (j + 1)] <- az[-1]

      y[1,2] <- 0
      y[nrow(y), (ncol(y)-1)] = 0

      mn <- rep(0, j)
      mn[1] = wq[1]
      mn[j] = wq[j]
      for (g in 2:(j-1)){
        mn[g] = alpha[g] * wq[g-1] + 2 * (1-alpha[g])*wq[g] + alpha[g] * wq[g+1]
      }

      wq  <- solve(y, mn)
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
        DKE = Hmacrophytes[dep]*(rho_plant* ahat * Cdplant) * daily_meteo["Uw",n]^3 * dt  *dx

      } else {
        PEprior = PE

        PE = abs(g *   depth[dep] *( depth[dep+1] - Zcv)  *
                   abs(calc_dens(u[dep+1])- mean(calc_dens(u[1:dep])))) + PEprior

        DKEprior = DKE
        DKE = Hmacrophytes[dep]*(rho_plant * ahat * Cdplant) * daily_meteo["Uw",n]^3 * dt  *dx + DKEprior

        KE = KE - DKE

      }
      if (PE > KE){
        maxdep = dep-1

        agents[which(agents <= maxdep)] = runif(length(agents[which(agents <= maxdep)]),
                                                min(seq(1,nrow(um))*dx), maxdep)

        break
      } else if (dep>1 & PE < KE ){
        u[(dep-1):dep] = (u[(dep-1):dep] %*% volume[(dep-1):dep])/sum(volume[(dep-1):dep])

        if (water.quality){
          wq[(dep-1):dep] = (wq[(dep-1):dep] %*% volume[(dep-1):dep])/sum(volume[(dep-1):dep])
        }
      }
      maxdep = dep

    }

    mix[n] <- KE/PE
    therm.z[n] <- maxdep

    if (pgdl_mode){
      um_mix[, n] <- u
    }


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
          u[dep:(dep+1)] = (u[dep:(dep+1)] %*% volume[dep:(dep+1)])/sum(volume[dep:(dep+1)])

          if (water.quality){
            wq[dep:(dep+1)] = (wq[dep:(dep+1)] %*% volume[dep:(dep+1)])/sum(volume[dep:(dep+1)])
          }
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
    if (pgdl_mode){
      um_conv[, n] <- u
    }


    ## (5) ICE FORMATION
    # according to Saloranta & Andersen (2007) and ice growth due to Stefan's law
    # (Leppäranta 1991)

    icep  = max(dt_iceon_avg,  (dt/86400))
    x = (dt/86400) / icep
    iceT = iceT * (1 - x) + u[1] * x

    # Him[ n] <- Hi
    K_snow <- 2.22362 * (rho_snow/1000)^1.885
    Tice = 0

    if ((iceT <= 0) == TRUE & Hi < Ice_min & daily_meteo["Tair",n] <= 0 & ice == FALSE){
      supercooled <- which(u < 0)
      initEnergy <- sum((0-u[supercooled])*area[supercooled] * dx * Cw)

      # if (ice != TRUE) {
      Hi <- Ice_min+(initEnergy/(910*L_ice))/max(area)
      # } else {
      #   if (daily_meteo["Tair",n] > 0){
      #     Tice <- 0 #
      #     Hi = Hi -max(c(0, meltP * dt*(((1-albedo)*daily_meteo["Jsw",n])+(longwave(cc = daily_meteo["CC",n], sigma = sigma, Tair = daily_meteo["Tair",n], ea = daily_meteo["ea",n], emissivity = emissivity, Jlw = daily_meteo["Jlw",n]) +
      #                                                                    backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps) +
      #                                                                    latent(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n],  RH = daily_meteo["RH",n], A = area, Cd = Cd) +
      #                                                                    sensible(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd)) )/(1000*L_ice)))
      #   } else {
      #     Tice <-  ((1/(10 * Hi)) * 0 + daily_meteo["Tair",n]) / (1 + (1/(10 * Hi)))
      #     Hi <- max(Ice_min, sqrt(Hi**2 + 2 * 2.1/(910 * L_ice)* (0 - Tice) * dt))
      #   }
      # }
      ice = TRUE
      if (Hi >= 0){
        u[supercooled] = 0
        u[1] = 0
      }
      # Him[ n] <- Hi
    } else if (ice == TRUE & Hi >= Ice_min) {

      Q_surf <- (u[1] - 0) * Cw * dx
      # u[supercooled] = 0
      u[1] = 0



      if (daily_meteo["Tair",n] > 0){
        Tice <- 0
        dHsnew <- 0

        if (Hs  > 0){
          dHs <- -max(c(0, meltP * dt * (((1 - IceSnowAttCoeff) * daily_meteo["Jsw",n])+(longwave(cc = daily_meteo["CC",n], sigma = sigma, Tair = daily_meteo["Tair",n], ea = daily_meteo["ea",n], emissivity = emissivity, Jlw = daily_meteo["Jlw",n]) +
                                                                                     backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps) +
                                                                                     latent(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n],  RH = daily_meteo["RH",n], A = area, Cd = Cd) +
                                                                                     sensible(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd)) )/(rho_fw*L_ice)))
          if ((Hs + dHs) < 0) {
            Hi_new <- Hi + (Hs + dHs) * (rho_fw/rho_ice)
          } else {
            Hi_new <- Hi
          }
        } else {
          dHs <- 0

          Hi_new = Hi -max(c(0, meltP * dt*(((1-IceSnowAttCoeff)*daily_meteo["Jsw",n])+(longwave(cc = daily_meteo["CC",n], sigma = sigma, Tair = daily_meteo["Tair",n], ea = daily_meteo["ea",n], emissivity = emissivity, Jlw = daily_meteo["Jlw",n]) +
                                                                                      backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps) +
                                                                                      latent(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n],  RH = daily_meteo["RH",n], A = area, Cd = Cd) +
                                                                                      sensible(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd)) )/(rho_ice*L_ice)))
          Hsi <- Hsi - max(c(0, meltP * dt *(((1-IceSnowAttCoeff)*daily_meteo["Jsw",n])+(longwave(cc = daily_meteo["CC",n], sigma = sigma, Tair = daily_meteo["Tair",n], ea = daily_meteo["ea",n], emissivity = emissivity, Jlw = daily_meteo["Jlw",n]) +
                                                                                          backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps) +
                                                                                          latent(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n],  RH = daily_meteo["RH",n], A = area, Cd = Cd) +
                                                                                          sensible(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n], RH = daily_meteo["RH",n], A = area, Cd = Cd)) )/(rho_ice*L_ice)))
          if (Hsi <= 0){Hsi <- 0}
        }

      } else {

        if (Hs > 0){
          K_snow <- 2.22362 * (rho_snow/1000)^1.885
          p <- (K_ice/K_snow) * (((rho_fw/rho_snow) * Hs ) / Hi)
          dHsi <- max(c(0, Hi * (rho_ice/rho_fw -1) + Hs))
          Hsi <- Hsi + dHsi
        } else {
          p = 1/(10 * Hi)
          dHsi <- 0
        }

        Tice <-  (p * 0 +  daily_meteo["Tair",n]) / (1 + p)
        Hi_new <- sqrt((Hi + dHsi)**2 + 2 * K_ice/(rho_ice * L_ice)* (0 - Tice) * dt)

        dHsnew <- daily_meteo["PP",n] * 1/(1000 * 86400) * dt
        dHs <- dHsnew  - dHsi * (rho_ice/rho_fw)
        dHsi <- 0


      }

      Hi <- Hi_new - (Q_surf/(rho_ice * L_ice))


      Q_surf <- 0
      Hs <- Hs + dHs

      if (Hi < Hsi){
        Hsi <- max(0, Hi)
      }

      if (Hs <= 0){
        Hs <- 0
        rho_snow <- rho_new_snow
      } else {
        rho_snow <- rho_snow * (Hs - dHsnew)/Hs + rho_new_snow * dHsnew/Hs
      }

      # Him[ n] <- Hi
    } else if (ice == TRUE & Hi < Ice_min){
      ice = FALSE
      # Him[ n] <- Hi
    }

    if (ice == FALSE){
      Hi = 0
      Hs = 0
      Hsi = 0
      # Him[n] = Hi
    }

    Him[ n] <- Hi
    Hsm[n] <- Hs
    Hsim[n] <- Hsi
    Ticem[n] <- Tice

    n2m[, n] <- n2
    um[, n] <- u

    if (water.quality){
      wqm[, n] <- wq
    }

    if (pgdl_mode){
      um_ice[, n] <- u

      meteo_pgdl[1, n] <-  daily_meteo["Tair",n]
      meteo_pgdl[2, n] <-   longwave(cc = daily_meteo["CC",n], sigma = sigma, Tair = daily_meteo["Tair",n], ea = daily_meteo["ea",n], emissivity = emissivity, Jlw = daily_meteo["Jlw",n]) -
        backscattering(emissivity = emissivity, sigma = sigma, Twater = un[1], eps = eps)
      meteo_pgdl[3, n] <-   latent(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n],  RH = daily_meteo["RH",n], A = area, Cd = Cd)
      meteo_pgdl[4, n] <-   sensible(Tair = daily_meteo["Tair",n], Twater = un[1], Uw = daily_meteo["Uw",n], p2 = p2, pa = daily_meteo["Pa",n], ea=daily_meteo["ea",n],  RH = daily_meteo["RH",n], A = area, Cd = Cd)
      meteo_pgdl[5, n] <-   daily_meteo["Jsw",n]
      meteo_pgdl[6, n] <-   kd
      meteo_pgdl[7, n] <-   shear
      meteo_pgdl[8, n] <-   tau
      meteo_pgdl[9, n] <-   max(area, na.rm = T)
    }


    ## (6) Individual-based Modeling
    # Modeled using passive transport approach from Hellweger & Bucci (2009)
    dKdz = rep(1, (nx))
    for (i in seq(2, nx-1)){
      dKdz[i] = ( abs(kzn[i+1] - kzn[i-1]) / (2 * dx) )
    }
    dKdz[nx] = ( abs(kzn[nx] - kzn[nx-1]) / abs(depth[nx] - depth[nx-1]) )
    dKdz[1] = ( abs(kzn[1+1] - kzn[1]) / abs(depth[1+1] - depth[1]) )

    intdKdz = approx(seq(1,nx)*dx, dKdz, agents)$y
    intkzn = approx(seq(1,nx)*dx, kzn,
                    agents + 0.5 * intdKdz * dt)$y

    agents <- agents + intdKdz * dt +
      rnorm(length(agents),0,1)  * sqrt(2 * intkzn * dt)

    agents[which(agents >= max(seq(1,nx)*dx))] = rep(max(seq(1,nx)*dx), length(agents[agents >= max(seq(1,nx)*dx)] )) + (rep(max(seq(1,nx)*dx),   length( agents[agents >= max(seq(1,nx)*dx)] )) - agents[which(agents >= max(seq(1,nx)*dx))])
    agents[which(agents <= min(seq(1,nx)*dx))] = rep(min(seq(1,nx)*dx), length(agents[agents <= min(seq(1,nx)*dx)] )) + (rep(min(seq(1,nx)*dx),  length(agents[agents <= min(seq(1,nx)*dx)] )) - agents[which(agents <= min(seq(1,nx)*dx))])

    magents[, n] <- agents

    loc <- as.data.frame(table(factor(round(agents), levels = seq(1, nx) * dx)), responseName = 'n')
    mloc[, n] <- loc$n

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

  df.avg.sim.do = NULL
  if (water.quality){
    df.sim.do <- data.frame(cbind(seq(startTime,endTime,dt), t(wqm)) )
    colnames(df.sim.do) <- c("datetime", as.character(paste0('do',seq(1,nrow(wqm))*dx)))
    df.sim.do$datetime <-   seq(startTime, endTime, dt)#/24/3600


    avg.epi.sim.do <- NULL
    avg.hyp.sim.do <- NULL
    avg.tot.sim.do <- NULL
    for (j in 1:nrow(df.z.df.sim)){
      d = df.sim.do[,-1]
      if (is.na(df.z.df.sim$z[j])){
        df.z.df.sim$z[j] = 1
      }
      avg.epi.sim.do <- append(avg.epi.sim.do,((as.numeric(d[j,1:df.z.df.sim$z[j]], na.rm = T) %*%
                                            area[1:df.z.df.sim$z[j]] )/
                                           sum(area[1:df.z.df.sim$z[j]])))
      avg.hyp.sim.do <- append(avg.hyp.sim.do,((as.numeric(d[j,df.z.df.sim$z[j]:ncol(d)], na.rm = T)%*%
                                            area[df.z.df.sim$z[j]:ncol(d)] )/
                                           sum(area[df.z.df.sim$z[j]:ncol(d)])))
      avg.tot.sim.do <- append(avg.tot.sim.do,((as.numeric(d[j,1:ncol(d)], na.rm = T)%*%
                                            area[1:ncol(d)] )/
                                           sum(area[1:ncol(d)])))
    }

    stratFlag = rep(NA, length = ncol(um))
    for (v in 1:length(stratFlag)){
      stratFlag[v] = ifelse((calc_dens(um[nx,v]) - calc_dens(um[1,v])) >= 0.1 &
                              mean(um[,v]) >= 4, 1, 0)
    }


    df.avg.sim.do <- data.frame('time' = df.sim$datetime,
                             'epi' = avg.epi.sim.do,
                             'hyp' = avg.hyp.sim.do,
                             'tot' = avg.tot.sim.do,
                             'stratFlag' = stratFlag,
                             'thermoclineDep' = bf.sim)
  }
  cat(greet[whichgreet, 2])
  dat = list('temp'  = um,
             'diff' = kzm,
             'mixing' = mix,
             'buoyancy' = n2m,
             'icethickness' = Him,
             'snowthickness' = Hsm,
             'snowicethickness' = Hsim,
             'ice_temp' = Ticem,
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
             'SEN' = SEN,
             'water.quality' = wqm,
             'average_do' = df.avg.sim.do,
             'agents' = magents,
             'location' = mloc)

  if (pgdl_mode){
    dat = list('temp'  = um,
               'diff' = kzm,
               'mixing' = mix,
               'buoyancy' = n2m,
               'icethickness' = Him,
               'snowthickness' = Hsm,
               'snowicethickness' = Hsim,
               'ice_temp' = Ticem,
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
               'SEN' = SEN,
               'water.quality' = wqm,
               'average_do' = df.avg.sim.do,
               'agents' = magents,
               'location' = mloc,
               'temp_initial' = um_initial,
               'temp_heat' = um_heat,
               'temp_diff' = um_diff,
               'temp_mix' = um_mix,
               'temp_conv' = um_conv,
               'temp_ice' = um_ice,
               'meteo_input' = meteo_pgdl,
               'buoyancy_pgdl' = n2_pgdl)
  }

  return(dat)


}

