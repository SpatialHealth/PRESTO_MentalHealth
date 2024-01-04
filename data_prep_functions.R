library(dplyr)
library(lubridate)


### Generates column names of NDVI measures 
generate_ndvi_colnames <- function(time_measure = "seasonal", stat = "max", years = 2012:2019, buffers = c(50,100,250,500)){
  output <- NA
  if(time_measure == "annual" & stat == "max"){
    #loop through buffers, apply across years 
    for(i in 1:length(buffers)){
      output <- c(output,paste0(stat, "_", buffers[i], "m_", years))
    }
  }
  
  if(time_measure == "seasonal"){
    seasons <- c("winter", "spring", "summer", "fall")
    #loop through seasons
    for(j in 1:length(seasons)){
      #loops through buffers, apply both across years 
      for(i in 1:length(buffers)){
        output <- c(output,paste0(stat, "_", seasons[j], "_", buffers[i], "m_", years))
      }
    }
  }
  return(output[-1])
}


### Shortcut for generating ALL NDVI colnames
generate_all_ndvi_colnames <- function(years = 2012:2019, buffers = c(50,100,250,500)){
  output = generate_ndvi_colnames("annual", stat = "max", years, buffers)
  output = append(output, generate_ndvi_colnames("seasonal", stat = "max", years, buffers))
  output = append(output, generate_ndvi_colnames("seasonal", stat = "mean", years, buffers))
  return(output)
}


### Annual NDVI assignment function 
annual_ndvi <- function(vector, id_colname = "StudyID", date_colname = "b_finisheddate", buffers = c(50,100,250,500)){
  # initialize output
  output <- vector[id_colname]
  
  # index using a standard format "max_[buffer]m_[year]"
  ndvi_string <- paste0("max_", as.character(buffers), "m_", year(as.Date(vector[date_colname])))
  output <- append(output, vector[ndvi_string])
  
  return(as.numeric(output))
}

### Seasonal NDVI assignment function 
seasonal_ndvi <- function(vector, stat="max", id_colname = "StudyID", date_colname = "b_finisheddate", buffers = c(50,100,250,500)){
  # initialize output
  output <- vector[id_colname]
  mon <- month(vector[date_colname])
  seasons <- c("winter", "spring", "summer", "fall")
  
  # sort by seasons, index using a standard format "max_[season]_[buffer]m_[year]"
  if( mon == 12 | mon == 1 | mon == 2 ){
    ndvi_string <- paste0(stat, "_", seasons[1], "_", as.character(buffers), "m_", year(as.Date(vector[date_colname])))
    output <- append(output, vector[ndvi_string])
  }else if(mon >= 3 & mon <= 5){
    ndvi_string <- paste0(stat, "_", seasons[2], "_", as.character(buffers), "m_", year(as.Date(vector[date_colname])))
    output <- append(output, vector[ndvi_string])
  }else if(mon >= 6 & mon <= 8){
    ndvi_string <- paste0(stat,"_", seasons[3], "_", as.character(buffers), "m_", year(as.Date(vector[date_colname])))
    output <- append(output, vector[ndvi_string])
  }else if(mon >= 9 & mon <= 11){
    ndvi_string <- paste0(stat, "_", seasons[4], "_", as.character(buffers), "m_", year(as.Date(vector[date_colname])))
    output <- append(output, vector[ndvi_string])
  }
  return(as.numeric(output))
}


### Apply NDVI assignment functions and format output 
add_ndvi <- function(data, id_colname = "StudyID", date_colname = "b_finisheddate", buffers = c(50,100,250,500)){
  # annual max 
  annual_max <- t(apply(data,1,annual_ndvi, id_colname=id_colname, date_colname=date_colname, buffers=buffers))
  annual_max <- as.data.frame(annual_max)
  colnames(annual_max) <- c(id_colname, paste0("ndvi_", buffers, "_amax") )
  data <- data[, !names(data) %in% c(generate_ndvi_colnames(time_measure="annual"))]
  
  # seasonal max 
  seasonal_max <- t(apply(data,1, seasonal_ndvi, stat="max", id_colname=id_colname, date_colname=date_colname, buffers=buffers))
  seasonal_max <- as.data.frame(seasonal_max)
  colnames(seasonal_max) <-c(id_colname, paste0("ndvi_", buffers, "_smax"))
  data <- data[, !names(data) %in% c(generate_ndvi_colnames(time_measure="seasonal", stat="max"))]
  
  # seasonal max 
  seasonal_mean <- t(apply(data,1, seasonal_ndvi, stat="mean", id_colname=id_colname, date_colname=date_colname, buffers=buffers))
  seasonal_mean <- as.data.frame(seasonal_mean)
  colnames(seasonal_mean) <-c(id_colname, paste0("ndvi_", buffers, "_smean"))
  data <- data[, !names(data) %in% c(generate_ndvi_colnames(time_measure="seasonal", stat="mean"))]
  
  output <- inner_join(data, annual_max, by = id_colname)
  output <- inner_join(output, seasonal_max, by = id_colname)
  output <- inner_join(output, seasonal_mean, by = id_colname)
  
  return(output)
}


### Separate Urban and Non-Urban - return a categorical 1=Urban, 0=Non-urban
add_urban_cat <- function(data){
  data$ct_urban_cat = NA
  for(i in 1:dim(data)[1]){
    if(is.na(data$ct_urban[i])){ 
      data$ct_urban[i] = NA
    }else if(data$Country[i] == "US"){
      data$ct_urban_cat[i] = as.numeric(data$ct_urban[i] >= 90)
    }else if(data$Country[i] == "Canada" & data$ct_urban[i] > 1){
      data$ct_urban_cat[i] = 1
    }else if(data$Country[i] == "Canada"){
      data$ct_urban_cat[i] = data$ct_urban[i]
    }
  }
  return(data)
}

### Compute Quantiles and IQR for all NDVI measures 
process_ndvi <- function(data, outpath){
  if(!is.data.frame(data) || !all(colnames(data) != "")){
    stop("Input must be a data frame with named columns")
  }
  
  ndvi_colnames <- names(data)[grep("ndvi", names(data))]
  
  breaks <- data.frame()
  
  for(i in 1:length(ndvi_colnames)){
    # NDVI quantiles
    q <- quantile(data[ndvi_colnames[i]], probs = c(0, 0.25, 0.5, 0.75, 1), names = TRUE, na.rm=TRUE)
    data[paste0(ndvi_colnames[i], "_quantile")] <- case_when(as.logical(data[,ndvi_colnames[i]] <= q[2]) ~ "1",
                                                             as.logical(data[,ndvi_colnames[i]] <= q[3]) ~ "2",
                                                             as.logical(data[,ndvi_colnames[i]] <= q[4]) ~ "3",
                                                             TRUE ~ "4")
    # NDVI IQR
    data[paste0(ndvi_colnames[i],"_iqr")] <- data[ndvi_colnames[i]] / IQR(unlist(data[ndvi_colnames[i]]), na.rm = TRUE)
    
    # IQR quantiles
    q_iqr <- quantile(data[paste0(ndvi_colnames[i],"_iqr")], probs = c(0, 0.25, 0.5, 0.75, 1), names = TRUE, na.rm=TRUE)
    data[paste0(ndvi_colnames[i], "_iqr_quantile")] <- case_when(as.logical(data[paste0(ndvi_colnames[i],"_iqr")] <= q_iqr[2]) ~ "1",
                                                             as.logical(data[paste0(ndvi_colnames[i],"_iqr")] <= q_iqr[3]) ~ "2",
                                                             as.logical(data[paste0(ndvi_colnames[i],"_iqr")] <= q_iqr[4]) ~ "3",
                                                             TRUE ~ "4")
    # save quantile ranges 
    breaks <- rbind(breaks, as.data.frame(t(c(ndvi_colnames[i], q))))
    breaks <- rbind(breaks, as.data.frame(t(c(paste0(ndvi_colnames[i],"_iqr"), q_iqr))))
  }
  write.csv(breaks, outpath)
  return(data)
}
