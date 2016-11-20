###############################################################
#
# WindPlume - an R script to conduct Monte Carlo analyses
# with a simple atmospheric dispersion model
#
###############################################################

### definitions ###

options(stringsAsFactors = FALSE)

library(reshape)

require(compiler)
enableJIT(3)

#setwd("C:/Users/wmcnab/Desktop/fall 2016 projects/air dispersion R")
setwd("D:/fall 2016 projects/air dispersion R")


WindVector <- function(intense_bin, orient_bin, orient, compass) {

  # posit wind velocity magnitude
  if (intense_bin == "CALM") {
    v_min <- 0.0
    v_max <- 0.5
    }
  else if (intense_bin == ">46") {
    v_min <- 46.
    v_max <- 46.
    }
  else {
    bounds <- strsplit(intense_bin, "-")[[1]]
    v_min <- as.numeric(bounds[1]) - 0.5
    v_max <- as.numeric(bounds[2]) + 0.5
    }
  v <- runif(1, v_min, v_max) * 0.447  # convert MPH to m/sec
    
  # posit wind direction (origin)
  if (orient_bin != "CALM"){    
    delin <- match(orient_bin, orient)
    angle_center <- compass[delin]
    angle_min <- angle_center - 22.5/2.
    angle_max <- angle_center + 22.5/2.
    }
  else {
    angle_min <- 0.
    angle_max <- 360.
    }
  angle <- runif(1, angle_min, angle_max)
  if (angle < 0.) {angle <- 360. + angle}

  return (c(v, angle*pi/180.))                 # convert angle to radians  
}


Stability <- function(v, day_f, stability){

  # assign atmospheric stability class

  if (v<2) {row <- 1}                 # wind velocity magnitude component
  else if (v>=2 & v<3) {row <- 2} 
  else if (v>=3 & v<4) {row <- 3}   
  else if (v>=4 & v<6) {row <- 4}     
  else {row <- 5}

  if (day_f <= 0.5) {col <- 1}        # solar irradiation component
  else {col <- 2}

  return (stability[row, col])
}


Dispersion <- function(x, stable_class){
  
  # assign transverse and vertical dispersion coefficients (standard terrain)
  
  switch(stable_class,
    A = {
      sigma_y <- 0.22*x / sqrt(1.0 + 0.0001*x)
      sigma_z <- 0.2*x      
      },
    B = {
      sigma_y <- 0.16*x / sqrt(1.0 + 0.0001*x)
      sigma_z <- 0.12*x     
      },   
    C = {
      sigma_y <- 0.11*x / sqrt(1.0 + 0.0001*x)
      sigma_z <- 0.08*x / sqrt(1.0 + 0.0002*x)    
      },    
    D = {
      sigma_y <- 0.08*x / sqrt(1.0 + 0.0001*x)
      sigma_z <- 0.06*x / sqrt(1.0 + 0.0015*x)      
      },    
    E = {
      sigma_y <- 0.06*x / (1.0 + 0.0001*x)
      sigma_z <- 0.03*x / (1.0 + 0.0003*x)       
      },    
    F_ = {
      sigma_y <- 0.04*x / (1.0 + 0.0001*x)
      sigma_z <- 0.016*x / (1.0 + 0.0003*x)      
      }        
    
    )
  
  return (c(sigma_y, sigma_z))
  
}


Rotate <- function(x, y, theta) {
  
  # rotate coordinate system about theta
  x_prime <- x*cos(theta) - y*sin(theta)
  y_prime <- x*sin(theta) + y*cos(theta)

  return (c(x_prime, y_prime))
  
  }


C <- function(x, y, z, u, H, flux, stable_class) {
  
  # calculate time-integrated concentration at x, y, z
  
  if (x > 0) {
    sigma <- Dispersion(x, stable_class)
    sigma_y = sigma[1]
    sigma_z = sigma[2]    
    f <- exp(-(y**2./(2.*sigma_y**2)))
    g1 <- exp(-((z - H)**2./(2.*sigma_z**2)))
    g2 <- exp(-((z + H)**2./(2.*sigma_z**2)))
    conc <- flux/u * f/(sigma_y*sqrt(2.*pi)) * (g1 + g2)/(sigma_z*sqrt(2.*pi))
    }
  else
    {conc <- 0.}
    
  return (conc)
  
  }


WindPlume <- function(){

  ### this is the main function of this script; specifies parameter sets, reads data, calls other functions, writes output

  # basic parameters
  num_trials <- 10000
  d_min <- 5.
  d_max <- 60.
  H <- 5.            
  z <- 2.
  Q <- 2.
  log_avg_C0 <- 1.0
  log_stdev_C0 <- 0.3
  
  # define working parameter sets
  compass <- c(90., 67.5, 45., 22.5, 0., 337.5, 315., 292.5, 270., 247.5, 225., 202.5, 180., 157.5, 135., 112.5, 0.)
  stability_table <- matrix(c("A", "B", "F_", "A", "C", "E", "B", "C", "D", "C", "D", "D", "C", "D", "D"),
                      nrow=5, ncol=3, byrow = TRUE)

  # read in meteorology data
  wind_df <- read.csv(file="wind.txt", sep="\t", header=TRUE)	
      
  # record bin labels for velocity ranges and direction
  orient <- colnames(wind_df)[-1]
  intense <- wind_df[,1]    

  # set up wind bins lookup table
  wind_bins_df <- melt(wind_df)
  colnames(wind_bins_df) <- c("intensity", "direction", "freq")
  wind_bins_df$direction <- as.character(wind_bins_df$direction)
  wind_bins_df <- subset(wind_bins_df, freq > 0.)
  wind_bins_df$cumul <- cumsum(wind_bins_df$freq)
  rownames(wind_bins_df) <- NULL

  # posit source term vector
  log_C0 <- rnorm(num_trials, log_avg_C0, log_stdev_C0)
  C0 <- 10.**log_C0                          
  flux <- Q * C0     

  # posit receptor location vectors
  d <- runif(num_trials, d_min, d_max)
  psi <- runif(num_trials, 0., 2.*pi)
  x <- d * cos(psi)
  y <- d * sin(psi)

  # create (empty) data frame to hold model results, one row per trial
  results_df <- data.frame("C0" = double(num_trials),
                  "x" = double(num_trials),
                  "y" = double(num_trials),
                  "w_origin" = character(num_trials),
                  "theta" = double(num_trials),
                  "v" = double(num_trials),
                  "s_class" = character(num_trials),
                  "conc" = double(num_trials),
                  stringsAsFactors=FALSE)
  
  # generate trials  

  for (i in 1:num_trials){
    
    # select random numbers to choose wind vector and solar irradiation conditions
    r <- runif(1, 0., 1.)
    day_f <- runif(1, 0., 1.)  
    floor <- ifelse(r > wind_bins_df$cumul, 1, 0)  
    bins_index <- sum(floor) + 1
    bins_index <- min(nrow(wind_bins_df), bins_index)

    # posit meteorology
    wvector <- WindVector(wind_bins_df$intensity[bins_index], wind_bins_df$direction[bins_index], orient, compass)
    v <- wvector[1]
    theta <- wvector[2]
    if (theta>pi) {theta_point <- theta - pi}    # point theta 180 degrees away for orientation of local positive x-axis
    else {theta_point <- theta + pi}
    rot_pt <- Rotate(x[i], y[i], -theta_point)
    x_prime <- rot_pt[1]
    y_prime <- rot_pt[2]    
    s <- Stability(v, day_f, stability_table)    

    # run Gaussian plume model
    conc <- C(x_prime, y_prime, z, v, H, flux[i], s)    
    
    # append results to data frame
    wind_dir <- wind_bins_df$direction[bins_index]
    results_df[i, ] <- c(C0[i], x[i], y[i], wind_dir, theta, v, s, conc)
    
    }

  write.csv(results_df, file = "results.csv") 

}

# run the air dispersion script

ptm <- proc.time()
WindPlume()
print (proc.time() - ptm)

# process results ...
num_trials <- 10000
d_min <- 5.
d_max <- 60.
H <- 5.            
z <- 2.
Q <- 2.
log_avg_C0 <- 1.0
log_stdev_C0 <- 0.3
model <- read.csv(file="results.csv", sep=",", header=TRUE)
C0 = 10.**log_avg_C0
count_steps <- c(0.01, 0.001, 0.0001, 1e-5, 1e-6)
num_bins = length(count_steps) + 1

# pie graph to summarize exposure probability
count_bin <- double(num_bins)
count_bin[1] <- sum(ifelse(model$conc>=count_steps[1]*C0, 1, 0))
for (i in 2:(num_bins-1)) {
  count_bin[i] <- sum(ifelse(model$conc>=count_steps[i]*C0, 1, 0)) - count_bin[i-1]
  }
count_bin[num_bins] <- num_trials - count_bin[num_bins-1]
lbls <- c(">0.01", "0.001 - 0.01", "0.0001 - 0.001", "0.00001 - 0.0001", "0.000001 - 0.00001", "<0.000001")
pct <- round(count_bin/num_trials*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(count_bin,labels = lbls, col=rainbow(length(lbls)),
    main="Normalized Conc Exposures")

# augment model results data frame with concentration bin designation (for plotting)
model$bin <- ifelse(model$conc>=count_steps[1]*C0, 1, 0)
for (i in 2:(num_bins-1)) {
  model$bin <- model$bin + ifelse((model$conc>=count_steps[i]*C0) & (model$conc<count_steps[i-1]*C0), i, 0)
}
model$bin <- model$bin + ifelse(model$conc<count_steps[num_bins-1]*C0, num_bins, 0)
write.csv(model, file = "model.csv")

print ("Done.")
