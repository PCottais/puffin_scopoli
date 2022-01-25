# LOADING PACKAGES
library(dplyr)
library(lubridate)
library(amt)

# LOADING DATA
# generalisation dataset
load("CALDIO_MPE_2021/generalisation datasets/grid_oceano_2011.Rdata")
load("CALDIO_MPE_2021/generalisation datasets/grid_oceano_Date_2011.Rdata")

# Scopoli dataset
load("CALDIO_MPE_2021/training datasets/CALDIO_2011_oceano.Rdata")
load("CALDIO_MPE_2021/training datasets/CALDIO_2012_oceano.Rdata")


# PART 1 : iSSA for each bird of 1 colony ----

# select all birds from Giraglia colony in 2011 (Scopoli's babies)
df_G <- ALL2011 %>% 
  filter(Site == "Giraglia") %>% 
  select(x = Longitude, y = Latitude, t = Time, id = ID,
         Bathy, SST1, logCHLA28)  # covariates

# nest data into a list column (Henry & Wickham, 2017)
df_G <- df_G %>% nest(-id)

# build tracks from orignial data frame
df_G <- df_G %>% 
  mutate(trks = lapply(data, function(d){
    amt::make_track(d, x, y, t, crs = "epsg:2154", all_cols = TRUE)
  }))
# inspect the sampling rate of the 30 individuals
df_G %>% mutate(smr = lapply(trks, summarize_sampling_rate)) %>% 
  select(id, smr) %>% unnest(cols = c(smr))

# build steps from tracks
df_G <- df_G %>% 
  mutate(ssf = lapply(trks, function(x){ x %>% 
      amt::track_resample(rate = minutes(10), tolerance = seconds(60)) %>% 
      amt::filter_min_n_burst(min_n = 3) %>% 
      amt::steps_by_burst(keep_col = "end") %>% 
      amt::random_steps(n_control = 9) %>% 
      mutate(log_sl = log(sl_)) %>% 
      amt::fit_issf(case_ ~ sl_ + ta_ + log_sl + strata(step_id_))
  }))

# select models
mod_G <- df_G %>% select(ssf)



# PART 2 : iSSA for all birds of 1 colony ----

# select all birds from Giraglia colony in 2011 (Scopoli's babies)
df_G <- ALL2011 %>% 
  filter(Site == "Giraglia") %>% 
  select(x = Longitude, y = Latitude, t = Time, id = ID,
         Bathy, SST28, logCHLA28)  # covariates

# tracks
trks <- df_G %>% 
  amt::make_track(.x = x, .y = y, .t = t, crs = "epsg:2154", all_cols = TRUE)

summarize_sampling_rate(trks)

# observed steps
stps <- trks %>% 
  amt::track_resample(rate = minutes(5), tolerance = seconds(60)) %>% 
  amt::filter_min_n_burst(min_n = 3) %>% 
  amt::steps_by_burst(keep_col = 'end')

# # deal with duplicate burst_id
# stps %>% as_tibble() %>% group_by(burst_) %>% 
#   summarize(nb_ind = n_distinct(id))
# 
# # 1st idea: modify `burst_` values so it is different for each bird
# stps <- stps %>% 
#   mutate(burst_ = paste(id, burst_, sep = ''))


# observed + random steps
stps <- stps %>% 
  amt::random_steps(n_control = 9) %>% 
  mutate(log_sl = log(sl_),
         burst_ = as.factor(burst_),
         step_id_ = as.factor(step_id_))


get_habitat <- function(steps, var = c("Bathy",'SST28','logCHLA28')){
  
  days <- unique(date(steps$t2_))
  # get habitat covariates for random ENDING points
  for (day in days){
    grid_day <- grid_oceano[[which(grid_Oceano_Date == day)]]
    matrx <- grid_day[, c("Longitude", "Latitude")] %>% as.matrix()
    filtered_stps <- steps[date(steps$t2_)==day,]
    nearest_id <- Rfast::dista(filtered_stps[,c("x2_", "y2_")],
                               matrx, k = 1, index = TRUE)
    steps[date(steps$t2_)==day, var] <- grid_day[nearest_id, var]
  }
  return(steps)
}

# get covariates of available steps from oceano grid
start <- Sys.time()
stps[stps$case_==FALSE,] <- stps %>% 
  filter(case_ == FALSE) %>% 
  get_habitat()
end <- Sys.time()
print(end-start)


# iSSA conditional regression model
mod_G <- stps %>% amt::fit_issf(case_ ~ sl_ + ta_ + log_sl + strata(step_id_) +
                                  Bathy + SST28 + logCHLA28)
mod_G$model



# PART 3 : iSSA for each colony ----

# function for the whole iSSA
iSSA <- function(colony = "Giraglia", year = 2011){
  if(year == 2011){
    data <- ALL2011
  }else{
    dataa <- ALL2012
  }
  
  # select all birds from a colony in 2011 (Scopoli's babies)
  df <- data %>% 
    filter(Site == colony) %>% 
    select(x = Longitude, y = Latitude, t = Time, id = ID,
           Bathy, SST28, logCHLA28)  # covariates
  
  # tracks
  trks <- df %>% 
    amt::make_track(.x = x, .y = y, .t = t, crs = "epsg:2154", all_cols = TRUE)
  
  # observed steps
  stps <- trks %>% 
    amt::track_resample(rate = minutes(5), tolerance = seconds(60)) %>% 
    amt::filter_min_n_burst(min_n = 3) %>% 
    amt::steps_by_burst(keep_col = 'end')
  
  # deal with duplicate burst_id
  stps %>% as_tibble() %>% group_by(burst_) %>%
    summarize(nb_ind = n_distinct(id))
  # 1st idea: modify `burst_` values so it is different for each bird
  stps <- stps %>%
    mutate(burst_ = paste(id, burst_, sep = ''))
  
  # observed + random steps
  stps <- stps %>% 
    amt::random_steps(n_control = 9) %>% 
    mutate(log_sl = log(sl_),
           burst_ = as.factor(burst_),
           step_id_ = as.factor(step_id_))
  
  # get covariates of available steps from oceano grid
  stps[stps$case_==FALSE,] <- stps %>% 
    filter(case_ == FALSE) %>% 
    get_habitat()
  
  # iSSA conditional regression model
  model <- stps %>% amt::fit_issf(case_ ~ sl_ + ta_ + log_sl + strata(step_id_) +
                                       Bathy + SST28 + logCHLA28)
  return(model)
}

# initialisation
models <- list()
start <- Sys.time()
for (site in unique(ALL2011$Site)){
  models[[site]] <- iSSA(colony = site, year = 2011)
}
end <- Sys.time()
print(end-start)
