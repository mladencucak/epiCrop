#' BlightR model.
#'
#' This function calculates potato late blight risk using BlightR model.
#' The risk needs to be calculated for a single location, so if the
#' calculation is necessary for multiple locations, data for each location needs to
#' run separately.
#'
#' @param data The weather data formated as \code{data frame}.
#' @param max_na Maximum proportion of missing values. Set to 0.01 by default.
#' @param temporal_res By default, the teporal resolution of the output is daily. By changing the argument \code{temporal_res = "hourly"} user will get daily values attached at 12.
#' @param model_parameters resoulution of the final data to be returned, daily or hourly. If hourly is selected, outputst are returned at noon.
#' @import stringr dplyr lubridate zoo
#' @keywords potato late blight, model, decision support
#' @return This function returns a \code{data.frame} including columns:
#' \itemize{
#'  \item date Date formated as "yyyy-mm-dd"
#'  \item spor Sporulation risk.
#'  \item spor_cond Conditions for sporulation.
#'  \item inf Infection risk.
#'  \item surv_prob Spore mortality.
#'  \item risk_si Risk as sum of infection and sporulation reduced by moratlity.
#'  \item risk_mi
#'  \item risk Risk as product of sporulation, spore survival and infection.
#' }
#' @examples
#' \donttest{
#' library(epiCrop)
#' weather <- epiCrop::weather
#' out <- system.time(BlightR(weather))
#' head(out)
#' }


BlightR <- function(data,
                    max_na = NULL,
                    temporal_res = "daily",
                    model_parameters = "default"
) {

  if(is.data.frame(data)==FALSE){
    stop("Weather data is not a data frame!")
  }

  #If the data has more than max_na proportion of missing values stop the function
  if(is.null(max_na)) max_na <- 0.01
  if (mean(is.na(data[, c("temp", "rhum")])) > max_na){
    stop("Percentage of missing values for relative humidity and temperature is higher than ", max_na*100, "%!")
  }

  #Check the format of the data frame
  if(any(apply(data[ , c("temp", "rhum", "sol_rad")],2, is.numeric))==FALSE ){
    stop("Temperature, relative humidty, solar radiation need to be numeric variables!")
  }

  # Sort columns
  if(!"doy"%in% colnames(data)) data$doy <- lubridate::yday(data$short_date)

  if(is.Date(data[ , "short_date"]) == FALSE){
    data[ , "short_date"]<-
      base::as.Date( as.character(data[ , "short_date"]),  "%Y-%m-%d")
  }

  if(nrow(data)/unique(data[ , "short_date"]) %>% length()!=24){
    stop("The weather data needs to have 24 rows per day!")
  }
  colnames <- c("short_date","doy","temp", "rhum")
  if(all(colnames %in% colnames(data))==FALSE) stop("Rename column names."); rm(colnames)


  # Parameters
  if(model_parameters == "default"){
  parameters <-
    data.frame(  TminInf=6,
                 ToptInf=14,
                 TmaxInf=26,
                 RfactInf=1,
                 ShapeInf=15,
                 TminInfDir=6,
                 ToptInfDir=21,
                 TmaxInfDir=26,
                 RfactInfDir=0.6,
                 ShapeInfDir=15,
                 RhminInf=86,
                 RhoptInf=95,
                 B0=2.37,
                 B1=0.45,
                 TminSpor=7,
                 ToptSpor=21,
                 TmaxSpor=25,
                 RfactSpor=1,
                 ShapeSpor=4,
                 KSpor=97049.81,
                 n0Spor=6.05E-04,
                 rSpor=1.734924,
                 spor_dur=10,
                 hr_before_spor=5,
                 hr_after_spor=5,
                 hr_after_inf=5)
  }else {
    if(!is.data.frame(model_parameters)){
      stop("Custom parameters need to be provided as a data frame.
           Each column name is one parameter and the fist row contains the parameter value.")
    }
  }



  # General functions
  ExtractCol <- function(data, column){
    if(!is.data.frame(data)){stop("The data is not of class data.frame!")}

    if(length(names(data)[str_detect(names(data), fixed(column, ignore_case=TRUE))])>0){
      col_name <- names(data)[str_detect(names(data), fixed(column, ignore_case=TRUE))]
    } else{
      stop("Variable doesnot exist in the data set.")
    }
    if (length(col_name)>1){
      col_name <-
        col_name [ column   == col_name]
    }
    if (column == "rh"){
      col_name <- "rhum"
    }
    return(col_name)
  }

  # Calculate the intersection of two curves
  # Function taken from
  # https://rdrr.io/github/andrewheiss/reconPlots/#vignettes
  # Calculate where two lines or curves intersect. Curves are defined as data
  # frames with x and y columns providing cartesian coordinates for the lines.
  # This function works on both linear and nonlinear curves.
  #

  curve_intersect <- function(curve1, curve2, empirical=TRUE, domain=NULL) {
    if (!empirical & missing(domain)) {
      stop("'domain' must be provided with non-empirical curves")
    }

    if (!empirical & (length(domain) != 2 | !is.numeric(domain))) {
      stop("'domain' must be a two-value numeric vector, like c(0, 10)")
    }

    if (empirical) {
      # Approximate the functional form of both curves
      curve1_f <- stats::approxfun(curve1$x, curve1$y, rule = 2)
      curve2_f <- stats::approxfun(curve2$x, curve2$y, rule = 2)

      # Calculate the intersection of curve 1 and curve 2 along the x-axis
      point_x <- stats::uniroot(function(x) curve1_f(x) - curve2_f(x),
                         c(min(curve1$x), max(curve1$x)))$root

      # Find where point_x is in curve 2
      point_y <- curve2_f(point_x)
    } else {
      # Calculate the intersection of curve 1 and curve 2 along the x-axis
      # within the given domain
      point_x <- stats::uniroot(function(x) curve1(x) - curve2(x), domain)$root

      # Find where point_x is in curve 2
      point_y <- curve2(point_x)
    }

    return(list(x = point_x, y = point_y))
  }
  ########################################################
  #Sporulation
  ########################################################
  Sporulation <-
    function(temp,
             rh,
             parameters
    ) {

      #Import parameters
      #Temp factor
      TminSpor <- parameters[, "TminSpor"] %>% as.numeric()
      ToptSpor <- parameters[, "ToptSpor"]%>% as.numeric()
      TmaxSpor <- parameters[, "TmaxSpor"]%>% as.numeric()
      RfactSpor <- parameters[, "RfactSpor"]%>% as.numeric()
      ShapeSpor <- parameters[, "ShapeSpor"]%>% as.numeric()

      #RH factor
      KSpor <- parameters[, "KSpor"]%>% as.numeric()
      n0Spor <- parameters[, "n0Spor"]%>% as.numeric()
      rSpor <- parameters[, "rSpor"]%>% as.numeric()


      #Calculate temp factor
      #Function to cacluate the hourly temperature sporulation risk
      SPORtemp <-
        sapply(temp, function(x) {
          sporulation_temperature <-
            RfactSpor * ((TmaxSpor - x) / (TmaxSpor - ToptSpor) * ((x - TminSpor) / (ToptSpor - TminSpor)) ^ ((ToptSpor - TminSpor) / (TmaxSpor - ToptSpor))) ^  ShapeSpor
          sporulation_temperature = ifelse(x < TminSpor | x > TmaxSpor, 0, sporulation_temperature)
          return(sporulation_temperature)
        }) %>% unlist()

      #Function to cacluate the hourly temperature sporulation risk
      rh <- ifelse(rh >=80, rh - 80, 0) #resize rh scale to 1-20
      SPORrh <-
        sapply(rh, function(x) {
          spor_rh_hour <-
            KSpor * n0Spor * exp(rSpor * x) / (KSpor + n0Spor * (exp(rSpor * x) - 1))
          # Change to proportion; sporulation is divided by sporulation at the equilibrium
          spor_rh_hour <-as.numeric(spor_rh_hour)/ as.numeric(KSpor)
          return(spor_rh_hour)
        })

      #Calculate the hourly sporulation risk
      SPOR <- c(SPORtemp * SPORrh)

      return(SPOR)
    }

  ########################################################
  #Survival
  ########################################################
  SolSurv <-
    function(sol,
             params_solsurv
    ) {

      B0 <- params_solsurv[, "B0"] %>% as.numeric()
      B1 <- params_solsurv[, "B1"] %>% as.numeric()

      Survival <- function(x) {
        Pr <- 1 / (1 + exp(-(B0 - B1 * x)))
      }

      daySURV <-  Survival(sol)

      return(daySURV)
    }



  ########################################################
  #Infection
  ########################################################

  Infection <-
    function(temp,
             rh,
             params_inf
    ) {

      #Import parameters
      #Temp factor zoospore
      TminInf <- params_inf[, "TminInf"] %>% as.numeric()
      ToptInf <- params_inf[, "ToptInf"]%>% as.numeric()
      TmaxInf <- params_inf[, "TmaxInf"]%>% as.numeric()
      RfactInf <- params_inf[, "RfactInf"]%>% as.numeric()
      ShapeInf <- params_inf[, "ShapeInf"]%>% as.numeric()
      #Temp factor direct
      TminInfDir <- params_inf[, "TminInf"]%>% as.numeric()
      ToptInfDir <- params_inf[, "ToptInfDir"]%>% as.numeric()
      TmaxInfDir <- params_inf[, "TmaxInf"]%>% as.numeric()


      RfactInfDir <- params_inf[, "RfactInfDir"]%>% as.numeric()
      ShapeInfDir <- params_inf[, "ShapeInfDir"]%>% as.numeric()


      # Temperature intersect between functions for mechanisms of infection
      CalcIntersect <- function(TminInf,ToptInf,TmaxInf,RfactInf,ShapeInf,
                                TminInfDir,ToptInfDir,TmaxInfDir,RfactInfDir,ShapeInfDir){
        temp <- c(0:34)

        zooINFtemp <- sapply(temp, function(x) {
          Infection_temperature <-
            RfactInf * ((TmaxInf - x) / (TmaxInf - ToptInf) * ((x - TminInf) / (ToptInf - TminInf)) ^ ((ToptInf - TminInf) / (TmaxInf - ToptInf))) ^  ShapeInf
          Infection_temperature = ifelse(x < TminInf | x > TmaxInf, 0, Infection_temperature)
        }) %>% unlist() %>% as.numeric()

        #Function to calculate the Infection temperature factor
        dirINFtemp <- sapply(temp, function(x) {
          Infection_temperature <-
            RfactInfDir * ((TmaxInfDir - x) / (TmaxInfDir - ToptInfDir) * ((x - TminInfDir) / (ToptInfDir - TminInfDir)) ^ ((ToptInfDir - TminInfDir) / (TmaxInfDir - ToptInfDir))) ^  ShapeInfDir
          Infection_temperature = ifelse(x < TminInfDir | x > TmaxInfDir, 0, Infection_temperature)
        }) %>% unlist()

        # Find the curve intersect
        x <- which(zooINFtemp==RfactInf) #peek of zoospore germiantion
        y <- which(dirINFtemp==RfactInfDir) #peek of direct germination

        # Straight lines (empirical)
        line1 <- data.frame(x = temp[x:y], y = dirINFtemp[x:y])
        line2 <- data.frame(x = temp[x:y], y = zooINFtemp[x:y])

        intersect <-  curve_intersect(line1, line2) %>% as.data.frame()

        InterTemp <- intersect[, "x"]
        return(InterTemp)
      }

      Tint <-
        CalcIntersect(TminInf,ToptInf,TmaxInf,RfactInf,ShapeInf,
                      TminInfDir,ToptInfDir,TmaxInfDir,RfactInfDir,ShapeInfDir)


      #Calculate temp factor
      #Function to calculate the Infection temperature factor
      zooINFtemp <- sapply(temp, function(x) {
        Infection_temperature <-
          RfactInf * ((TmaxInf - x) / (TmaxInf - ToptInf) * ((x - TminInf) / (ToptInf - TminInf)) ^ ((ToptInf - TminInf) / (TmaxInf - ToptInf))) ^  ShapeInf
        Infection_temperature = ifelse(x < TminInf |
                                         x > TmaxInf, 0, Infection_temperature)
      }) %>% unlist()

      #Function to cacluate the Infection temperature factor
      dirINFtemp <- sapply(temp, function(x) {
        Infection_temperature <-
          RfactInfDir * ((TmaxInfDir - x) / (TmaxInfDir - ToptInfDir) * ((x - TminInfDir) / (ToptInfDir - TminInfDir)) ^ ((ToptInfDir - TminInfDir) / (TmaxInfDir - ToptInfDir))) ^  ShapeInfDir

        Infection_temperature <-
          ifelse(x < TminInfDir |x > TmaxInfDir, 0, Infection_temperature)
      })%>% unlist()


      INFtemp <-
        ifelse(temp <= Tint, zooINFtemp, dirINFtemp )

      #RH factor
      RhminInf <- params_inf[, "RhminInf"]
      RhoptInf <- params_inf[, "RhoptInf"]


      # Calculate the hourly RH infection factor
      #fit linear model
      x <- c(RhminInf, RhoptInf) %>% as.numeric()
      coefs <- stats::lm(c(0, 1) ~ x)[["coefficients"]]

      INFlw <- sapply(rh, function (rhum_val) {
        if (rhum_val > RhminInf & rhum_val <= RhoptInf) {
          SporRh <-  coefs[["(Intercept)"]] + coefs[["x"]] * rhum_val
        } else if (rhum_val > RhoptInf & rhum_val <= 100) {
          SporRh <- 1
        } else {
          SporRh <- 0
        }
      })

      #Calculate the daily infection
      dayINF = INFtemp * INFlw
      return(dayINF)
    }



  ########################################################
  #Sun info
  ########################################################

  #https://github.com/cran/HelpersMG/blob/master/R/sun.info.R


  SunInfo <- function(date, latitude, longitude,
                      #Ireland uses Irish Standard Time (IST, UTC+01:00 in the summer months
                      #and Greenwich Mean Time (UTC+0) in the winter period.
                      UTC_zone
  ){

    d <- as.numeric(as.POSIXlt(date)$yday)+1
    Lat <- latitude
    Long <- longitude

    ## d is the day of year
    ## Lat is latitude in decimal degrees
    ## Long is longitude in decimal degrees (negative == West)

    ##This method is copied from:
    ##Teets, D.A. 2003. Predicting sunrise and sunset times.
    ##  The College Mathematics Journal 34(4):317-321.

    ## At the default location the estimates of sunrise and sunset are within
    ## seven minutes of the correct times (http://aa.usno.navy.mil/data/docs/RS_OneYear.php)
    ## with a mean of 2.4 minutes error.

    ## Function to convert degrees to radians
    rad <- function(x) pi*x/180

    ##Radius of the earth (km)
    R=6378

    ##Radians between the xy-plane and the ecliptic plane
    epsilon=rad(23.45)

    ##Convert observer's latitude to radians
    L=rad(Lat)

    ## Calculate offset of sunrise based on longitude (min)
    ## If Long is negative, then the mod represents degrees West of
    ## a standard time meridian, so timing of sunrise and sunset should
    ## be made later.
    timezone = -4*(abs(Long)%%15)*sign(Long)

    ## The earth's mean distance from the sun (km)
    r = 149598000

    theta = 2*pi/365.25*(d-80)

    z.s = r*sin(theta)*sin(epsilon)
    r.p = sqrt(r^2-z.s^2)

    t0 = 1440/(2*pi)*acos((R-z.s*sin(L))/(r.p*cos(L)))

    ##a kludge adjustment for the radius of the sun
    that = t0+5

    ## Adjust "noon" for the fact that the earth's orbit is not circular:
    n = 720-10*sin(4*pi*(d-80)/365.25)+8*sin(2*pi*d/365.25)

    ## now sunrise and sunset are:
    sunrise = (n-that+timezone)/60
    sunset = (n+that+timezone)/60

    UTC <- (((7.5+Long)%%360)) %/% 15
    if (UTC>12) {UTC <- 12-UTC;tz <- "Etc/GMT"} else {tz <- "Etc/GMT+"}
    tz <- paste0(tz, UTC)

    df <-
      data.frame(
        sunrise = sunrise,
        sunset = sunset,
        day.length = sunset - sunrise,
        date.time.sunrise = as.POSIXlt(format(date, "%Y-%m-%d"), tz =
                                         tz) + sunrise * 60 * 60,
        date.time.sunset = as.POSIXlt(format(date, "%Y-%m-%d"), tz =
                                        tz) + sunset * 60 * 60
      )

    sunrise.UTC <-
      as.POSIXlt(
        format(df$date.time.sunrise, format = "%Y-%m-%d %H:%M:%S"),
        tz = "UTC",
        use.tz = TRUE
      )
    sunrise.UTC.dec <-
      sunrise.UTC$hour + sunrise.UTC$min / 60 + sunrise.UTC$sec / 3600
    sunset.UTC <-
      as.POSIXlt(
        format(df$date.time.sunset, format = "%Y-%m-%d %H:%M:%S"),
        tz = "UTC",
        use.tz = TRUE
      )
    sunset.UTC.dec <-
      sunset.UTC$hour + sunset.UTC$min / 60 + sunset.UTC$sec / 3600

    #Return
    df <-
      data.frame(
        # df,
        sunrise = sunrise.UTC + c(UTC_zone*60*60),
        sunset = sunset.UTC  + c(UTC_zone*60*60)
        # time.sunrise.UTC = sunrise.UTC.dec,
        # time.sunset.UTC = sunset.UTC.dec
      )

    return(df)
  }


  ########################################################
  #Calculation of cut-off times
  ########################################################

  GetTimes <- function(data) {

    # Set the sunset and sunrise manually at 20/6
    # data[data$hour == 20, "daytime"] <- "sunset"
    # data[data$hour == 6, "daytime"] <- "sunrise"
    if(all(!str_detect(colnames(data), fixed("lon", ignore_case=TRUE)))){  stop("No Longitude reference or it is not named: 'lon'!")}
    if(all(!str_detect(colnames(data), fixed("lat", ignore_case=TRUE)))){  stop("No Latitude reference or it is not named: 'lat'!")}

    lat <- data[[ExtractCol(data, "lat")]]
    lon <- data[[ ExtractCol(data, "lon")]]

    #Extract dates and geo-coordinates
    tempdf <-
      data.frame(
        date = unique(data$short_date),
        lat = unique(lat),
        long = unique(lon)
      )

    fun_ls <-  list()
    for (i in seq_along(1:nrow(tempdf))) {
      fun_ls[[i]] <- SunInfo(tempdf[i,"date"],
                             tempdf[i, "lat"],
                             tempdf[i,"long"],
                             UTC_zone=1)
    }

    tempdf <-
      fun_ls %>%
      dplyr::bind_rows() %>%
      dplyr::bind_cols(tempdf, .) %>%
      dplyr::mutate(sunset_hr = lubridate::hour(sunset),
             sunrise_hr = lubridate::hour(sunrise),
             doy = lubridate::yday(date)) %>%
      dplyr::select(c("doy", "sunrise_hr", "sunset_hr"))

    #Retun sunset and sunrise times as a data frame
    return(tempdf)

  }


  ################################
  #Model Run
  #################################

  #Calculate initiation/termiantion points for sporulation and infection
  # Calculate an approximate sunrise and sunset times
  #https://www.timeanddate.com/sun/@2963597
  timedf <- GetTimes(data = data)

  #Define number of days for risk estimation
  # Risk can be calculated for from day 2, because sporulaiton is calculated for
  # the previous night. Same goes for the infection - it extends to the following morning
  begin <- unique(data$doy)[2]
  end <- unique(data$doy)[length(unique(data$doy))-1]

  result_ls <- list()


  #Define the start of sporulation
  # i = 136
  for(i in c(begin:end)){
    sporstart <-
      timedf[ timedf$doy == c(i-1), "sunset_hr"] -
      parameters[ ,"hr_before_spor"] %>% as.numeric()
    infstop <- timedf[ timedf$doy == c(i+1), "sunrise_hr"] + parameters[ ,"hr_after_inf"] %>% as.numeric()

    daydf <-
      do.call("rbind",
              list(
                data[with(data, hour >= sporstart & doy == c(i-1)),],# from the afternoon the day before
                data[with(data,  doy == i),],
                data[with(data, hour <= infstop & doy == c(i+1)),] #ending in the morning on the day after
              )
      )

    #If there are missing values just return rows with NAs
    if (any(is.na(daydf[, c("temp")])|is.na(daydf[, c("rhum")]))) {
      final <- data.frame(
        doy = i,
        spor = NA,
        spor_cond = NA,
        inf = NA,
        surv_prob =NA,
        inf_sol = NA,
        risk_si = NA,
        risk_mi = NA,
        risk = NA
      )

    } else{ #start of if is not NA statement

      ############################
      #Sporulation
      ############################
      #Calculate sporulation
      daydf$spor <-
        Sporulation(daydf$temp, daydf$rhum, parameters)

      # Stop the sporulation n(hr_after_spor) hours after sunrise
      sporstop <-  timedf[timedf$doy == i, "sunrise_hr"] + parameters[ ,"hr_after_spor"] %>% as.numeric()
      daydf[daydf$doy == i & daydf$hour >= sporstop |  daydf$doy == c(i+1), c("spor")] <- 0

      #cumulative sporulation per event
      daydf$spor_sum <-
        stats::ave(dplyr::coalesce(daydf$spor, 0), data.table::rleid(zoo::na.locf(daydf$spor != 0,maxgap = 3)), FUN = cumsum)

      # daydf[daydf$doy == i & daydf$hour < sporstop |  daydf$doy == c(i-1), c("spor", "spor_sum")]

      # check if there is 10 hours for sporulation
      criteria<- as.numeric(daydf$spor>0)

      #cumulative sum of hours that met the criteria with restart at zero
      criteria_sum <-  stats::ave(criteria, cumsum(criteria == 0), FUN = cumsum)

      risk <- rep(0, nrow(daydf))

      hours <- parameters[ ,"spor_dur"] %>% as.numeric()
      criteria_met  <- as.numeric( criteria_sum >= hours )
      idx           <- which(criteria_sum == hours)



      # If the sporulation criteria was not met calculate risk based on mortality and infection periods
      # The infection period starts in the morning and lasts until next morning
      if(sum(criteria_met)==0){

        inf_start <-
          timedf[timedf$doy == i, "sunrise_hr"] + parameters[ ,"hr_after_spor"] %>% as.numeric()
        idx <- which(daydf$doy == i & daydf$hour == inf_start)

        ############################
        #Infection
        ############################
        # Calculate the infection period
        # Starts after the conditions for sporulation have been met
        # Sum of Infection and sporulation for each hour reduced by survival

        daydf$inf <- 0
        daydf[idx:nrow(daydf), "inf"] <-
          Infection(temp = daydf[idx:nrow(daydf), "temp"],
                    rh = daydf[idx:nrow(daydf), "rhum"],
                    params_inf = parameters)



        #Cumulative sum of the infection after sporulation requirement has been met
        #The initial value is total sporulation of the day
        # daydf[idx-1,"inf"] <- max(daydf$spor_sum)
        # daydf[c(idx-1):nrow(daydf),"cumul_inf"] <-
        #   cumsum( daydf[c(idx-1):nrow(daydf),"inf"])

        #Cumulative sum of the infection after sporulation requirement has been met
        #The accumulation breaks if the conditions arent met for more than infstop hours
        infstop <- 3

        infr <-
          daydf[c(idx-1):nrow(daydf),"inf"]
        infwin <- rep(1, infstop)
        infr <- c(infr, infwin) %>% unlist()
        infrr <- infr
        for (k in c(1:c(length(infr)-infstop))){
          # i = 7
          if(sum(infr[k :k+infstop])>0){
            infrr[k] <- infr[k]
          }else{
            infrr[k:length(infr)] <-0
            break
          }
        }

        infrr <- infrr[1:c(length(infr)-infstop)]

        # Finally, Infection is a sum of
        daydf[c(idx-1):nrow(daydf),"cumul_inf"] <- cumsum(infrr)

        ############################
        #Survival
        ############################
        # Calculate the mortality of spores due to solar raditaion if there is solar radiation data

        #Airborne sporangia survival
        #The estimated probablity of spore survival is calculated using
        # The spore load was calculated as a product of total daily sporulation risk and the
        # probability of spornagia survival as a function of solar radiation
        solar_rad <-
          sum(daydf[daydf$doy == i, "sol_rad"])
        surv_prob <- SolSurv(solar_rad, parameters)

        ############################
        #Risk calculation
        ############################
        # Risk estimation based on mortality and infection risk estimation
        risk_mi <-
          max(daydf$cumul_inf, na.rm = TRUE)*
          surv_prob


        final <- data.frame(
          doy = i,
          spor = max(daydf$spor_sum, na.rm = TRUE),
          spor_cond = "no",
          inf = 0,
          surv_prob = 0,
          risk_si = 0,
          risk_mi = risk_mi,
          risk = 0,
          stringsAsFactors =FALSE
        )
        result_ls [[i]]<-final



      }else {

        ############################
        #Infection
        ############################
        # Calculate the infection period
        # Starts after the conditions for sporulation have been met
        # Sum of Infection and sporulation for each hour reduced by survival
        daydf$inf <- 0
        daydf[idx:nrow(daydf), "inf"] <-
          Infection(temp = daydf[idx:nrow(daydf), "temp"],
                    rh = daydf[idx:nrow(daydf), "rhum"],
                    params_inf = parameters)



        #Cumulative sum of the infection after sporulation requirement has been met
        #The accumulation breaks if the conditions arent met for more than infstop hours
        infstop <- 3


        infr <-
          daydf[c(idx-1):nrow(daydf),"inf"]
        infwin <- rep(1, infstop)
        infr <- c(infr, infwin) %>% unlist()
        infrr <- infr
        for (k in c(1:c(length(infr)-infstop))){
          # i = 7
          if(sum(infr[k :k+infstop])>0){
            infrr[k] <- infr[k]
          }else{
            infrr[k:length(infr)] <-0
            break
          }
        }

        infrr <- infrr[1:c(length(infr)-infstop)]

        # Finally, Infection is a sum of hourly infection
        daydf[c(idx-1):nrow(daydf),"cumul_inf"] <- cumsum(infrr)

        ############################
        #Survival
        ############################
        # Calculate the mortality of spores due to solar raditaion if there is solar radiation data

        #Airborne sporangia survival
        #The estimated probablity of spore survival is calculated using
        # The spore load was calculated as a product of total daily sporulation risk and the
        # probability of spornagia survival as a function of solar radiation

        solar_rad <-
          sum(daydf[daydf$doy == i, "sol_rad"])
        surv_prob <- SolSurv(solar_rad, parameters)



        ############################
        #Risk calculation
        ############################
        # Risk is calculated as a product of sporulation, spore mortality
        # (due to the solar radiation) and infection risk
        risk <-
          max(daydf$spor_sum, na.rm = TRUE) *
          max(daydf$cumul_inf, na.rm = TRUE)*
          surv_prob

        # Risk estimation based on sporulation and infection risk estimation
        risk_si <-
          max(daydf$spor_sum, na.rm = TRUE) *
          max(daydf$cumul_inf, na.rm = TRUE)

        # Risk estimation based on mortality and infection risk estimation
        risk_mi <-
          max(daydf$cumul_inf, na.rm = TRUE)*
          surv_prob



        final <-
          data.frame(
            doy = i,
            spor = max(daydf$spor_sum, na.rm = TRUE),
            spor_cond = "yes",
            inf = max(daydf$cumul_inf, na.rm = TRUE),
            surv_prob =surv_prob,
            risk_si = risk_si,
            risk_mi = risk_mi,
            risk = risk,
            stringsAsFactors =FALSE
          )
      }
      result_ls [[i]]<-final
    }
  }

  fin <-
    dplyr::left_join( timedf, bind_rows(result_ls),by = "doy") %>%
    select(-c("sunrise_hr", "sunset_hr"))

  ###############################
  # Final results
  ###############################


  if(temporal_res == "daily") {
    final <-
      data %>%
      dplyr::group_by(doy) %>%
      dplyr::summarise(date  = unique(short_date)) %>%
      dplyr::left_join(., fin, by = "doy")
  }


  DailyAt <-  function(x, time) {
    xx <- c(sapply(x, function(x) c(rep(NA,23),x)))
    x <- c(xx[time:length(xx)], xx[2:time])
    return(x)
  }

  if(temporal_res == "hourly"){
    fin_hourly <-
      data.frame(
        spor =   DailyAt(fin$spor, 12),
        spor_inf =   DailyAt(fin$spor_cond, 12),
        inf = DailyAt(fin$inf, 12),
        surv_prob = DailyAt(fin$surv_prob, 12),
        risk_mi = DailyAt(fin$risk_mi, 12),
        risk_si = DailyAt(fin$risk_si, 12),
        risk =  DailyAt(fin$risk, 12),
      )
    fin_hourly <-
      dplyr::bind_cols(data[, c("short_date", "hour")],fin_hourly)
    final <- fin_hourly
  }

  final <- final %>% select(-doy)
  return(final)
}






