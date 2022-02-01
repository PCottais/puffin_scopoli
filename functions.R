get_habitat <- function(steps, var = c("Bathy",'SST28','logCHLA28')){
  days <- unique(date(steps$t2_))
  # get habitat covariates for random ENDING points
  for (day in as.character(days)){
    if (substr(day, 1, 4) == "2011"){
      grid_day <- grid_oceano_2011[[which(grid_Oceano_Date_2011 == day)]]
    }else{
      grid_day <- grid_oceano_2012[[which(grid_Oceano_Date_2012 == day)]]
    }
    
    # # get oceanographic values of nearest "sea" cells
    # grid_day <- grid_day %>% filter(!is.na(SST1), !is.na(Bathy))
    
    matrx <- grid_day[, c("Longitude", "Latitude")] %>% as.matrix()
    filtered_stps <- steps[date(steps$t2_)==day,]
    nearest_id <- Rfast::dista(filtered_stps[,c("x2_", "y2_")],
                               matrx, k = 1, index = TRUE)
    steps[date(steps$t2_)==day, var] <- grid_day[nearest_id, var]
  }
  return(steps)
}


iSSA_steps <- function(colony, rmv_near_coast = FALSE,
                 covariates = c("Bathy",'SST28','logCHLA28')){
  
  # select all birds from a colony in 2011 (Scopoli's babies)
  df <- ALL %>% 
    filter(Site == colony) %>% 
    select(x = Longitude, y = Latitude, t = Time, id = ID)
  df[, covariates] <- NA  # covariates
  
  # tracks
  trks <- df %>% 
    amt::make_track(.x = x, .y = y, .t = t, crs = "epsg:2154", all_cols = TRUE)
  
  # observed steps
  stps <- trks %>% 
    amt::track_resample(rate = minutes(3), tolerance = seconds(60)) %>% 
    amt::filter_min_n_burst(min_n = 3) %>% 
    amt::steps_by_burst(keep_col = 'end')
  
  # deal with duplicate burst_id
  # stps %>% as_tibble() %>% group_by(burst_) %>%
  #   summarize(nb_ind = n_distinct(id))
  # 1st idea: modify `burst_` values so it is different for each bird
  stps <- stps %>%
    mutate(burst_ = paste(id, burst_, sep = ''))
  
  # observed + random steps
  stps <- stps %>% 
    amt::random_steps(n_control = 3) %>% 
    mutate(log_sl = log(sl_),
           burst_ = as.factor(burst_),
           step_id_ = as.factor(step_id_))
  
  # get covariates of available steps from oceano grid
  stps <- stps %>% 
    # filter(case_ == FALSE) %>% 
    get_habitat(var = covariates)
  
  return(stps)
}
