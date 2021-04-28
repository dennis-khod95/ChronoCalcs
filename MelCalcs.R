# input packages
library(tidyverse)

#########################################################################
# filter_na
## Removal of rows with NA values in the Melatonin column
## Necessary for the creation of smoothed melatonin values using a spline
## Default na_threshold is 8. Keeps only datasets with 8 or more observations. Based on the 24 hourly samples. 
### na_threshold can be adjusted as necessary. 

filter_na <- function(dataset, id, melatonin, session, na_threshold=8) {

  session <- enquo(session)
  id <- enquo(id)
  melatonin <- enquo(melatonin)
  dataset <- dataset %>% 
    rename(id=!!id,
           session=!!session, 
           melatonin=!!melatonin)

  dat_sub <- dataset[!is.na(dataset$melatonin), ]
  dat_sub <- dat_sub %>%
    dplyr::group_by(session, id) %>%
    dplyr::filter(n() >= na_threshold) %>% 
    dplyr::ungroup()
  return(dat_sub)
}



#########################################################################
# ob_time
## Creates observation time variable, which starts at the designated start time, and sequences by 1
## Allows for the ordering of data based on when samples were taken, rather than by clock time
## Product of this function (ob.hour) is used in most other functions


ob_time <- function(dataset, clock.hour, start.time) {
  clock.hour <- enquo(clock.hour)
  dataset <- dataset %>% 
    rename(clock.hour=!!clock.hour)  
  
  dataset <- dataset %>% 
    mutate(ob.hour = ifelse((clock.hour - start.time) >= 0, 
                            (clock.hour - start.time), 
                            (24 - start.time) + clock.hour))
  return(dataset)
}




#########################################################################
# spline_fit
## Creation of smoothed melatonin values using a spline
## Subets by ID (participant number) and by session (can be season or sampling session)
## spar is set to 0.3 as a default to balance over/underfitting. Can adjust with spar_val.

spline_fit <- function(dataset, id, melatonin, session, ob.hour, spar_val=0.3) {
  session <- enquo(session)
  id <- enquo(id)
  melatonin <- enquo(melatonin)
  ob.hour <- enquo(ob.hour)
  dataset <- dataset %>% 
    rename(id=!!id,
           session=!!session, 
           melatonin=!!melatonin, 
           ob.hour=!!ob.hour)  
  
  counter<-0
    seas.vec<- unique(as.character(dataset$session))
    
    for(l in 1:length(seas.vec)){
      data<- subset(dataset, session==seas.vec[l])
      indiv.vec<- unique(as.character(data$id))
      for(m in 1:length(indiv.vec)){
        sub<- subset(data,id==indiv.vec[m])
        sub<- sub[order(sub$ob.hour), ]
        
        sp <- smooth.spline(sub$ob.hour, sub$melatonin, spar = spar_val)
        
        sub$spline_mel <- sp$y
        
        counter<- counter + 1
        if(counter==1){newdat <- sub}
        else{newdat <- rbind(newdat, sub)}
      }}
    return(newdat)
}

# example
# spline_dat <- spline_fit(dataset = sumtest, id = "id", 
#          melatonin = "melatonin", session = "session", ob.hour = "ob.hour")





#########################################################################
# threshold
## Uses the single maximum value, and the average of the 3 minimum values occuring earlier than the maximum. 
## Recommended to use smoothed melatonin values generated with spline_fit()
## Overall formula: ((maximum - minimum) * threshold) + minimum
## onset_threshold set to default of 0.25, offset threshold set to default of 0.5
### Can be adjusted manually as necessary

threshold <- function(dataset, id, melatonin_val, session, ob.hour, 
                       onset_threshold=0.25, offset_threshold=0.5) {
  session <- enquo(session)
  id <- enquo(id)
  melatonin_val <- enquo(melatonin_val)
  ob.hour <- enquo(ob.hour)
  dataset <- dataset %>% 
    rename(id=!!id,
           session=!!session, 
           melatonin_val=!!melatonin_val, 
           ob.hour=!!ob.hour)  
  
  counter<-0
  seas.vec<- unique(as.character(dataset$session))

for(l in 1:length(seas.vec)){
  data<- subset(dataset,session==seas.vec[l])
  indiv.vec<- unique(as.character(data$id))
  for(m in 1:length(indiv.vec)){
    sub<- subset(data,id==indiv.vec[m])
    
    ##### find single max value
    maxi <- sub[which.max(sub$melatonin_val),]
    max_mel <- maxi$melatonin_val
    max_time <- maxi$ob.hour
    
    ##### find 3 minimum values and average
    # filter to only consider values occuring before the maximum occurs
    submin <- sub %>% 
      dplyr::filter(ob.hour <= max_time & ob.hour >=1)
    # order by melatonin concentration
    submin <- submin %>%
      dplyr::arrange(melatonin_val) 
    # take lowest 3 values (nonconsecutive)
    submini <- submin %>% 
      head(n=3)
    # average of the 3 min values
    mean_min_mel <- mean(submini$melatonin_val) 
    
    df <- data.frame(matrix(ncol = 5, nrow = 1))
    x <- c("participant_id", "session_id", 
           "max_mel", "max_time", 
           "mean_min_mel")
    colnames(df) <- x
    df$participant_id <- paste(indiv.vec[m])
    df$session_id <- paste(seas.vec[l])
    df$max_mel <- max_mel
    df$max_time <- max_time
    df$mean_min_mel <- mean_min_mel
    
    counter<- counter + 1
    if(counter==1){quarterdat <- df}
    else{quarterdat <- rbind(quarterdat, df)}
  }}
# 25% for onset and 50% for offset
quarterdat$onset_thresh <- ((quarterdat$max_mel - quarterdat$mean_min_mel) * onset_threshold) + 
  quarterdat$mean_min_mel
quarterdat$offset_thresh <- ((quarterdat$max_mel - quarterdat$mean_min_mel) * offset_threshold) + 
  quarterdat$mean_min_mel
return(quarterdat)
}

# example
# thresholds <- threshold(dataset = spline_dat, id = "id", 
#          melatonin = "spline_mel", session = "session", ob.hour = "ob.hour")



#########################################################################
# dlmo_onset
## Finds the time of the DLMO (onset)
## output provides DLMO hour and DLMO minute
## dlmo_onset is based on predetermined expectations of melatonin trends
### individuals with abnormal trends must be removed prior to running this function and calculated seperately 
### use the dlmo_prefilter function to determine if any individuals exhibit abnormal trends



#
summm <- spline_dat %>% 
  filter(id != 11)
summm <- summm %>% 
  group_by(id, session) %>% 
  filter(id != 1 | session != "Autumn") %>% 
  ungroup()
#


dlmo_onset <- function(dataset, id, melatonin_val, session, ob.hour, 
                       threshold_data, threshold_id, threshold_session, onset_threshold, max_time) {

  session <- enquo(session)
  id <- enquo(id)
  melatonin_val <- enquo(melatonin_val)
  ob.hour <- enquo(ob.hour)
  
  threshold_session <- enquo(threshold_session)
  threshold_id <- enquo(threshold_id)
  onset_threshold <- enquo(onset_threshold)
  max_time <- enquo(max_time)
  
  dataset <- dataset %>% 
    rename(id=!!id,
           session=!!session, 
           melatonin_val=!!melatonin_val, 
           ob.hour=!!ob.hour)  
  threshold_data <- threshold_data %>% 
    rename(threshold_session=!!threshold_session, 
           threshold_id=!!threshold_id, 
           onset_threshold=!!onset_threshold, 
           max_time=!!max_time)

  counter<-0
  seas.vec<- unique(as.character(dataset$session))

for(l in 1:length(seas.vec)){
  data<- subset(dataset,session==seas.vec[l])
  indiv.vec<- unique(as.character(data$id))
  for(m in 1:length(indiv.vec)){
    sub <- subset(data,id==indiv.vec[m])
    
    thresh <- threshold_data %>% 
      dplyr::filter(threshold_session == seas.vec[l] & threshold_id == indiv.vec[m])
    
    # find 25% threshold value
    threshold <- thresh$onset_threshold
    maxtime <- thresh$max_time
    
    # filter values greater than threshold
    upper <- sub %>% 
      dplyr::filter(melatonin_val >= threshold)
    # find earliest occurance of a value above the threshold
    upper <- upper %>% 
      dplyr::arrange(ob.hour) %>% 
      head(n=1)
    # save time and smoothed melatonin value
    uppertime <- upper$ob.hour
    uppermel <- upper$melatonin_val
    
    # filter values lower than threshold
    # also filtering to values occuring earlier in time than the maximum
    lower <- sub %>% 
      dplyr::filter(melatonin_val <= threshold & ob.hour < maxtime)
    # find latest occurance of a value below the threshold
    lower <- lower %>% 
      dplyr::arrange(ob.hour) %>% 
      tail(n=1)
    # save time and smoothed melatonin value
    lowertime <- lower$ob.hour
    lowermel <- lower$melatonin_val
    
    # calculations 
    # find change in smoothed melatonin
    change_mel <- uppermel - lowermel
    change_time <- (uppertime - lowertime) * 60 
    changepermin <- change_mel / change_time
    
    # create data frame
    df <- data.frame(matrix(ncol = 3, nrow = change_time+1))
    x <- c("mel", "minute", "ob.hour")
    colnames(df) <- x
    # sequencing melatonin value at each minute of the hour based on above calculations
    df$mel <- seq(from=lowermel, to=uppermel, by=changepermin)
    df$minute <- seq(from=0, to=change_time, by=1)
    df$ob.hour <- lowertime
    
    # filter to values below the threshold value and order by time
    df_filt <- df %>% 
      filter(mel <= threshold)
    # find latest occurance of a value below the threshold
    df_filt <- df_filt %>% 
      dplyr::arrange(minute) %>% 
      tail(n=1)
    
    # save time and smoothed melatonin value
    DLMO_hour <- df_filt$ob.hour
    DLMO_min <- df_filt$minute
    
    # save to dataframe
    dat <- data.frame(matrix(ncol = 6, nrow = 1))
    y <- c("participant_id", "session_id", 
           "DLMO_hour", "DLMO_min", "changepermin_on", "max_time")
    colnames(dat) <- y
    dat$participant_id <- paste(indiv.vec[m]) 
    dat$session_id <- paste(seas.vec[l])
    dat$DLMO_hour <- DLMO_hour
    dat$DLMO_min <- DLMO_min
    dat$changepermin_on <- changepermin
    dat$max_time <- maxtime
    
    counter<- counter + 1
    if(counter==1){DLMO <- dat}
    else{DLMO <- rbind(DLMO, dat)}
  }}
  DLMO$thresh_time_on <- DLMO$DLMO_hour + (DLMO$DLMO_min / 60)
return(DLMO)
}



# example
# onset <- dlmo_onset(dataset = summm, id = "id", melatonin = "spline_mel", session = "session", 
#                 ob.hour = "ob.hour", threshold_data = thresholds, threshold_id = "participant_id", 
#                 threshold_session = "session_id", onset_threshold = "onset_thresh", max_time = "max_time")




#########################################################################
# dlmo_offset
## Finds the time of the DLMO (offset)
## output provides DLMOff hour and DLMOff minute
## dlmo_offset is based on predetermined expectations of melatonin trends
### individuals with abnormal trends must be removed prior to running this function and calculated seperately 
### use the dlmo_prefilter function to determine if any individuals exhibit abnormal trends

dlmo_offset <- function(dataset, id, melatonin_val, session, ob.hour, 
                       threshold_data, threshold_id, threshold_session, offset_threshold, max_time) {
  
  session <- enquo(session)
  id <- enquo(id)
  melatonin_val <- enquo(melatonin_val)
  ob.hour <- enquo(ob.hour)
  
  threshold_session <- enquo(threshold_session)
  threshold_id <- enquo(threshold_id)
  offset_threshold <- enquo(offset_threshold)
  max_time <- enquo(max_time)
  
  dataset <- dataset %>% 
    rename(id=!!id,
           session=!!session, 
           melatonin_val=!!melatonin_val, 
           ob.hour=!!ob.hour)  
  threshold_data <- threshold_data %>% 
    rename(threshold_session=!!threshold_session, 
           threshold_id=!!threshold_id, 
           offset_threshold=!!offset_threshold, 
           max_time=!!max_time)

counter<-0
seas.vec<- unique(as.character(dataset$session))

for(l in 1:length(seas.vec)){
  data<- subset(dataset,session==seas.vec[l])
  indiv.vec<- unique(as.character(data$id))
  for(m in 1:length(indiv.vec)){
    sub <- subset(data,id==indiv.vec[m])
    
    thresh <- threshold_data %>% 
      filter(threshold_session == seas.vec[l] & threshold_id == indiv.vec[m])
    
    # find 50% threshold value
    threshold <- thresh$offset_threshold
    maxtime <- thresh$max_time
    
    # filter values lower than threshold
    # also filtering to values occuring later in time than the maximum
    upper <- sub %>% 
      dplyr::filter(melatonin_val <= threshold & ob.hour > maxtime)
    # find first occurance of a value below the threshold
    upper <- upper %>% 
      arrange(ob.hour) %>% 
      head(n=1)
    # save time and smoothed melatonin value
    uppertime <- upper$ob.hour
    uppermel <- upper$melatonin_val
    
    # filter values greater than threshold
    lower <- sub %>% 
      dplyr::filter(melatonin_val >= threshold & ob.hour >= maxtime)
    # find last occurance of a value above the threshold
    lower <- lower %>% 
      arrange(ob.hour) %>% 
      tail(n=1)
    # save time and smoothed melatonin value
    lowertime <- lower$ob.hour
    lowermel <- lower$melatonin_val
    
    # calculations 
    # find change in smoothed melatonin
    change_mel <- lowermel - uppermel
    change_time <- (uppertime - lowertime) * 60 
    changepermin <- abs(change_mel / change_time)
    
    # create data frame
    df <- data.frame(matrix(ncol = 3, nrow = change_time+1))
    x <- c("mel", "minute", "ob.hour")
    colnames(df) <- x
    # sequencing melatonin value at each minute of the hour based on above calculations
    #df$mel <- seq(from=uppermel, to=lowermel, by=-changepermin)
    #df$minute <- seq(from=0, to=change_time, by=1)
    #df$ob.hour <- lowertime
    
    ####################### seq by doesnt work with negative numbers
    # need alternate way
    df$minute <- seq(from=0, to=change_time, by=1)
    df$mel <- (lowermel - (df$minute * changepermin))
    df$ob.hour <- lowertime
    #######################
    
    # filter to values below the threshold value and order by time
    df_filt <- df %>% 
      dplyr::filter(mel <= threshold)
    # find earliest occurance of a value below the threshold
    df_filt <- df_filt %>% 
      arrange(minute) %>% 
      head(n=1)
    
    # save time and smoothed melatonin value
    DLMO_hour <- df_filt$ob.hour
    DLMO_min <- df_filt$minute
    
    # save to dataframe
    dat <- data.frame(matrix(ncol = 5, nrow = 1))
    y <- c("participant_id", "session_id", 
           "DLMOff_hour", "DLMOff_min", "changepermin_off")
    colnames(dat) <- y
    dat$participant_id <- paste(indiv.vec[m]) # for the loop
    dat$session_id <- paste(seas.vec[l])
    dat$DLMOff_hour <- DLMO_hour
    dat$DLMOff_min <- DLMO_min
    dat$changepermin_off <- changepermin
    
    counter<- counter + 1
    if(counter==1){DLMOff <- dat}
    else{DLMOff <- rbind(DLMOff, dat)}
  }}
DLMOff$thresh_time_off <- DLMOff$DLMOff_hour + (DLMOff$DLMOff_min / 60)
return(DLMOff)
}

# example 
# offset <- dlmo_offset(dataset = summm, id = "id", melatonin_val = "spline_mel", session = "session", 
#                    ob.hour = "ob.hour", threshold_data = thresholds, threshold_id = "participant_id", 
#                    threshold_session = "session_id", offset_threshold = "offset_thresh", max_time = "max_time")




#########################################################################
# dlmo_prefilter
## Used to determine if any individuals exhibit abnormal melatonin trends incompatible with the dlmo calculations
### Individuals displaying a negative value for change_time_offset or change_time_onset 
### need to be removed prior to the corresponding DLMO calculations

dlmo_prefilter <- function(dataset, id, melatonin_val, session, ob.hour, 
                           threshold_data, threshold_id, threshold_session, 
                           offset_threshold, onset_threshold, max_time) {
  
  session <- enquo(session)
  id <- enquo(id)
  melatonin_val <- enquo(melatonin_val)
  ob.hour <- enquo(ob.hour)
  
  threshold_session <- enquo(threshold_session)
  threshold_id <- enquo(threshold_id)
  offset_threshold <- enquo(offset_threshold)
  onset_threshold <- enquo(onset_threshold)
  max_time <- enquo(max_time)
  
  dataset <- dataset %>% 
    rename(id=!!id,
           session=!!session, 
           melatonin_val=!!melatonin_val, 
           ob.hour=!!ob.hour)  
  threshold_data <- threshold_data %>% 
    rename(threshold_session=!!threshold_session, 
           threshold_id=!!threshold_id, 
           offset_threshold=!!offset_threshold, 
           onset_threshold=!!onset_threshold, 
           max_time=!!max_time)

  counter<-0
  seas.vec<- unique(as.character(dataset$session))

for(l in 1:length(seas.vec)){
  data<- subset(dataset,session==seas.vec[l])
  indiv.vec<- unique(as.character(data$id))
  for(m in 1:length(indiv.vec)){
    sub <- subset(data,id==indiv.vec[m])
    
    thresh <- threshold_data %>% 
      dplyr::filter(threshold_session == seas.vec[l] & threshold_id == indiv.vec[m])
    
    # find 25% threshold value
    threshold <- thresh$onset_threshold
    maxtime <- thresh$max_time
    
    # filter values greater than threshold
    upper <- sub %>% 
      dplyr::filter(melatonin_val >= threshold)
    # find earliest occurance of a value above the threshold
    upper <- upper %>% 
      dplyr::arrange(ob.hour) %>% 
      head(n=1)
    # save time and smoothed melatonin value
    uppertime <- upper$ob.hour
    uppermel <- upper$melatonin_val
    
    # filter values lower than threshold
    # also filtering to values occuring earlier in time than the maximum
    lower <- sub %>% 
      dplyr::filter(melatonin_val <= threshold & ob.hour < maxtime)
    # find latest occurance of a value below the threshold
    lower <- lower %>% 
      dplyr::arrange(ob.hour) %>% 
      tail(n=1)
    # save time and smoothed melatonin value
    lowertime <- lower$ob.hour
    lowermel <- lower$melatonin_val
    
    # calculations 
    # find change in smoothed melatonin
    change_mel_onset <- uppermel - lowermel
    change_time_onset <- (uppertime - lowertime) 
    
    
    # offset
    threshold <- thresh$offset_threshold
    maxtime <- thresh$max_time
    
    # filter values lower than threshold
    # also filtering to values occuring later in time than the maximum
    upper <- sub %>% 
      dplyr::filter(melatonin_val <= threshold & ob.hour > maxtime)
    # find first occurance of a value below the threshold
    upper <- upper %>% 
      arrange(ob.hour) %>% 
      head(n=1)
    # save time and smoothed melatonin value
    uppertime <- upper$ob.hour
    uppermel <- upper$melatonin_val
    
    # filter values greater than threshold
    lower <- sub %>% 
      dplyr::filter(melatonin_val >= threshold & ob.hour >= maxtime)
    # find last occurance of a value above the threshold
    lower <- lower %>% 
      arrange(ob.hour) %>% 
      tail(n=1)
    # save time and smoothed melatonin value
    lowertime <- lower$ob.hour
    lowermel <- lower$melatonin_val
    
    # calculations 
    # find change in smoothed melatonin
    change_mel_offset <- lowermel - uppermel
    change_time_offset <- (uppertime - lowertime)
    
    # save to dataframe
    dat <- data.frame(matrix(ncol = 6, nrow = 1))
    y <- c("participant_id", "session_id", 
           "change_mel_offset", "change_time_offset", 
           "change_mel_onset", "change_time_onset")
    colnames(dat) <- y
    dat$participant_id <- paste(indiv.vec[m]) 
    dat$session_id <- paste(seas.vec[l])
    dat$change_time_offset <- change_time_offset
    dat$change_time_onset <- change_time_onset
    dat$change_mel_onset <- change_mel_onset
    dat$change_mel_offset <- change_mel_offset
    
    counter<- counter + 1
    if(counter==1){DLMO <- dat}
    else{DLMO <- rbind(DLMO, dat)}
  }}
return(DLMO)
}

# example
# dlmo_check <- dlmo_prefilter(dataset = spline_dat, id = "id", melatonin_val = "spline_mel", session = "session", 
#                      ob.hour = "ob.hour", threshold_data = thresholds, threshold_id = "participant_id", 
#                      threshold_session = "session_id", offset_threshold = "offset_thresh", 
#                      onset_threshold = "onset_thresh", max_time = "max_time")




#########################################################################
# veri_plot
## Visualization meant to verify the accuracy of melatonin threshold calculuations and dlmo calculations
## Based on ggplot
## Input is complicated, requiring fitted melatonin data, dlmo onset/offset, and melatonin threshold datasets



veri_plot <- function(melatonin_data, melatonin_raw, melatonin_fit, 
                           melatonin_id, melatonin_session, current_session, observation.hour, clock.hour, 
                           onset_data, onset_threshold_time, onset_id, onset_session, 
                           offset_data, offset_threshold_time, offset_id, offset_session, 
                           threshold_data, offset_threshold, onset_threshold, threshold_id, threshold_session) {
  
melatonin_session <- enquo(melatonin_session)
#current_session <- enquo(current_session)
melatonin_id <- enquo(melatonin_id)
melatonin_raw <- enquo(melatonin_raw)
melatonin_fit <- enquo(melatonin_fit)
observation.hour <- enquo(observation.hour)
clock.hour <- enquo(clock.hour)

onset_id <- enquo(onset_id)
onset_session <- enquo(onset_session)
offset_id <- enquo(offset_id)
offset_session <- enquo(offset_session)
threshold_id <- enquo(threshold_id)
threshold_session <- enquo(threshold_session)

offset_threshold_time <- enquo(offset_threshold_time)
onset_threshold_time <- enquo(onset_threshold_time)

offset_threshold <- enquo(offset_threshold)
onset_threshold <- enquo(onset_threshold)

melatonin_data <- melatonin_data %>% 
  dplyr::rename(melatonin_id=!!melatonin_id,
                melatonin_session=!!melatonin_session, 
                melatonin_raw=!!melatonin_raw, 
                melatonin_fit=!!melatonin_fit, 
                clock.hour=!!clock.hour, 
                observation.hour=!!observation.hour)    
onset_data <- onset_data %>% 
  dplyr::rename(onset_id=!!onset_id,
                onset_session=!!onset_session, 
                onset_threshold_time=!!onset_threshold_time)    
offset_data <- offset_data %>% 
  dplyr::rename(offset_id=!!offset_id,
                offset_session=!!offset_session, 
                offset_threshold_time=!!offset_threshold_time)   
threshold_data <- threshold_data %>% 
  dplyr::rename(threshold_id=!!threshold_id,
                threshold_session=!!threshold_session, 
                offset_threshold=!!offset_threshold, 
                onset_threshold=!!onset_threshold) 

mel_seas <- melatonin_data %>%
  dplyr::filter(melatonin_session == current_session)
on_seas <- onset_data %>%
  dplyr::filter(onset_session == current_session)
off_seas <- offset_data %>%
  dplyr::filter(offset_session == current_session)
thresh_seas <- threshold_data %>%
  dplyr::filter(threshold_session == current_session)

indi.vec<-unique(as.character(on_seas$onset_id))

plot_lista = list()
for(m in 1:length(indi.vec)){
  sub <- subset(mel_seas, melatonin_id==indi.vec[m])
  f <- on_seas %>% 
    filter(onset_id == indi.vec[m])
  g <- off_seas %>% 
    filter(offset_id == indi.vec[m])
  h <- thresh_seas %>% 
    filter(threshold_id == indi.vec[m])
  
  sub$on_time <- f$onset_threshold_time
  sub$off_time <- g$offset_threshold_time
  
  sub$new25 <- h$onset_threshold
  sub$new50 <- h$offset_threshold

  q <- ggplot(sub, aes(x=observation.hour, y=melatonin_raw)) +
    geom_point(size=0.5) +    
#    scale_x_continuous(labels = clock.hour) + 
    
#    scale_x_continuous(labels = c(start.time, start.time+4, start.time+8, 
#                       start.time+12, start.time+16, start.time+20, start.time+24)) + 
  
    xlab("Time (observation hour)") + 
    ylab("Plasma Melatonin") + 
    geom_line(aes(x=observation.hour, y=melatonin_fit), color="dodgerblue4") +
    
    geom_line(aes(x=observation.hour, y=new25, color="Onset")) + 
    geom_line(aes(x=observation.hour, y=new50, color="Offset")) +
    
    geom_vline(xintercept = sub$on_time, linetype="dashed") + 
    geom_vline(xintercept = sub$off_time, linetype="dashed") + 
    
    ggtitle(paste(indi.vec[m], current_session, sep = "_"))
  plot_lista[[m]] <- q
}
return(plot_lista)
}


# sumplots <- veri_plot(melatonin_data = summm, melatonin_raw = "melatonin", melatonin_fit = "spline_mel", 
#          observation.hour = "ob.hour", clock.hour = "clock_hour", melatonin_id = "id", 
#          melatonin_session = "session", current_session = "Summer", onset_data = onset, 
#          onset_threshold_time = "thresh_time_on", onset_id = "participant_id", onset_session = "session_id", 
#          offset_data = offset, offset_threshold_time = "thresh_time_off", offset_id = "participant_id", 
#          offset_session = "session_id", threshold_data = thresholds, offset_threshold = "offset_thresh", 
#          onset_threshold = "onset_thresh", threshold_id = "participant_id", threshold_session = "session_id")




#########################################################################
# veri_plot_print
## Saves the output of veri_plot function as a pdf
## plot_list is the output from veri_plot_print
## pdf_name is your desired name for the pdf file (format: "filename.pdf")
### Make sure to end with .pdf and put the entire name in quotation marks


veri_plot_print <- function(plot_list, pdf_name) {
  pdf(pdf_name)
  for (i in 1:length(plot_list)) {
    print(plot_list[[i]])}
  dev.off()
}

# veri_plot_print(plot_list = sumplots, pdf_name = "plot_print_test.pdf")





#########################################################################
# wright_plot
## Useful visualization to show distribution of DLMO, melatonin maximums, and DLMOff across the dataset
## saved output of the wright_plot function can be edited with standard ggplot commands

wright_plot <- function(onset_data, offset_data, id, session, 
                           offset_threshold_time, onset_threshold_time, max_time) {

  session <- enquo(session)
  id <- enquo(id)
  offset_threshold_time <- enquo(offset_threshold_time)
  onset_threshold_time <- enquo(onset_threshold_time)
  max_time <- enquo(max_time)
  
  onset_data <- onset_data %>% 
    dplyr::rename(id=!!id,
           session=!!session, 
           max_time=!!max_time, 
           onset_threshold_time=!!onset_threshold_time)    
  
  offset_data <- offset_data %>% 
    dplyr::rename(id=!!id,
           session=!!session, 
           offset_threshold_time=!!offset_threshold_time)   
  
dat <- merge(onset_data, offset_data, by=c("id", "session"))

max <- dat %>% 
  select(id, session, max_time)
max$timing <- max$max_time
max$max_time <- NULL
max$measure <- paste("Max")

on <- dat %>% 
  select(id, session, onset_threshold_time)
on$timing <- on$onset_threshold_time
on$onset_threshold_time <- NULL
on$measure <- paste("Onset")

off <- dat %>% 
  select(id, session, offset_threshold_time)
off$timing <- off$offset_threshold_time
off$offset_threshold_time <- NULL
off$measure <- paste("Offset")

d <- rbind(max, on)
plotdat <- rbind(d, off)

plotdat <- plotdat %>%
  mutate(time_order = ifelse(measure == "Onset", 1, 
                             ifelse(measure == "Max", 2, 3)))

plotdat$measure  <- with(plotdat, reorder(measure, time_order))


wplot <- ggplot(plotdat, aes(x=session, y=timing, fill=measure)) + 
  geom_boxplot(position="identity", width=0.3) + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position="top", 
        legend.title = element_blank()) + 
  ylab("") + 
  xlab("") + 
  scale_y_continuous(breaks = seq(4,20,4),
                     labels = paste0(c("17:00", "21:00", "01:00", "05:00", "09:00"))) + 
  scale_fill_manual(values=c("#7BCBB4", # onset 
                             "#4050A2", # max
                             "#E69F00")) 

return(wplot)
}


wrightplot <- wright_plot(onset_data = onset, offset_data = offset, id = "participant_id", 
                          session = "session_id", offset_threshold_time = "thresh_time_off", 
                          onset_threshold_time = "thresh_time_on", max_time = "max_time")




