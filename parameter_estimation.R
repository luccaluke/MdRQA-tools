rm(list=ls())
source("estimate_emb_del.R")
source("mdrqa.R")
Data <- read.table("matb_fnirs/matb_fnirs_ROIs.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
channels <- c('DLPFC', 'PMC', 'STG', 'FPA', 'OFC')
trials <- as.matrix(unique(Data[, c("participant", "block", "task")]))

params <- data.frame(delay=numeric(), embedding=numeric(), radius=numeric(), RR=numeric())

## Delay & Embedding
for (i in seq(1, nrow(trials))){
  participant <- trials$participant[i]
  block <- trials$block[i]
  task <- trials$task[i]
  
  fnirsData <- na.omit(data.frame(Data[Data$participant==participant & Data$block==block & Data$task==task, channels]))
  
  DE <- estimate_del_emb(Data=fnirsData, max_delay=30, delay_inc=1, max_dim=20, dim_inc=1, fnn_thresh=.01, theiler_window=1, max_lost=90)
  
  params[i, "delay"] <- DE$del
  params[i, "embedding"] <- DE$emb
}


delay <- unique(params$delay)[which.max(tabulate(match(params$delay, unique(params$delay))))]
embedding <- unique(params$embedding)[which.max(tabulate(match(params$embedding, unique(params$embedding))))]


## Radius
found <- FALSE #flag tracking if estimate has been found
inc <- TRUE #keeps track of if we want to increase or decrease radius during search iterations
rec_list <- NULL #list of recurrence rates for radius being tested

#radius search parameters
rad <- .25 #radius starting estimate
rad_inc <- .1 #initial value by which to increment/decrement radius each iteration of search
rec_target <- 5 #target RR
rec_target_low <- 4.9 #lower bound RR
rec_target_high <- 5.1 #upper bound RR

#while not found, keep searching
while(!found) {
  rec_list <- NULL #reinitialize
  print(paste("rad:", rad))
  
  #loop through all teams to get median RR for current radius
  for (i in seq(1, nrow(trials))){
    participant <- trials$participant[i]
    block <- trials$block[i]
    task <- trials$task[i]
    fnirsData <- na.omit(data.frame(Data[Data$participant==participant & Data$block==block & Data$task==task, channels]))
    
    #get RR
    RQ <- mdrqa(data=fnirsData, emb=embedding, del=delay, rad=rad, norm = 'euc')
    rr <- RQ$REC
    rec_list <- c(rec_list, rr)
  }
  avg_rec <- median(rec_list, na.rm=TRUE) #median RR
  
  if (avg_rec < rec_target_low){ #rr too low
    if (!inc) {
      rad_inc <- rad_inc / 2 #decrease increment value for fine-grained search
      inc <- TRUE 
    }
    if (rad + rad_inc > 1) {
      found <- TRUE
    } else {
      rad <- rad + rad_inc #increment radius
    }
  }
  if (avg_rec > rec_target_high){ #rr too high
    if (inc){
      rad_inc <- rad_inc / 2 #decrease increment value for fine-grained search
      inc <- FALSE
    }
    if (rad < rad_inc) {
      found <- TRUE
    } else {
      rad <- rad - rad_inc #decrement radius
    }
  }
  if (avg_rec > rec_target_low & avg_rec < rec_target_high){ #target RR found
    found <- TRUE
  }
  print(paste("Avg RR:", avg_rec)) #print progress
}

params$radius <- rad #final radius estimate
params$radius <- rad
params$RR <- rec_list

write.csv(params, "matb_fnirs_mdrqa_parameters.csv")