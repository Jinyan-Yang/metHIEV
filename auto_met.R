make.met.func(startDate = format(Sys.Date() - 1,"%Y-%m-01"),
              endDate = Sys.Date() - 1)


# # get packages
# packages <- c("zoo", "HIEv", "dplyr", "doBy", "tidyr",'devtools')
# if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
#   install.packages(setdiff(packages, rownames(installed.packages())))  
# }
# 
# if("HIEv" %in% rownames(installed.packages()) == FALSE){
#   library(devtools)
#   install_bitbucket("remkoduursma/HIEv")
# }
# 
# library(zoo)
# library(HIEv)
# library(dplyr)
# library(doBy)
# library(tidyr)
# 
# setToken()
# if(!dir.exists('download/')){dir.create('download/')}
# download.path <- file.path("download/")
# setToPath(download.path)
# 
# # function to calculate sturate vapor pressure based on air temperature
# saturate.vp.func <- function(Tc,a=6.12,m=7.59,Tn=240.73){
#   # T = Temperature (°C)
#   VPS <- a*10^(m*Tc/(Tc+Tn))
#   return(VPS)
# }
# 
# # function to get met file
# get.met.func <- function(startDate,
#                          endDate){
#   
#   # get met data
#   # ppt from rain out shelter
#   ros15 <- downloadTOA5("ROS_WS_Table15",
#                         startDate = startDate,
#                         endDate = endDate,
#                         keepFiles = TRUE)
#   # other met from each ring
#   ros15_30 <- as.data.frame(summarize(group_by(ros15,DateTime=nearestTimeStep(DateTime,30)),
#                                       Rain=sum(Rain_mm_Tot, na.rm=TRUE)))
#   # get co2 data
#   CaData.ls <- list()
#   
#   for (Ring in 1:6){
#     
#     fn <- sprintf("FACE_AUTO_R%s_FCPLOGG_R",Ring)
#     Rawdata <- downloadTOA5(fn, 
#                             maxnfiles = 600, 
#                             rowbind = FALSE,
#                             keepFiles = TRUE,
#                             startDate = startDate,
#                             endDate = endDate)      
#     
#     #make all wrong format to numeric so that we can rbine
#     Rawdata <- lapply(Rawdata, function(x){
#       x$WindSpeed <- as.numeric(as.character(x$WindSpeed))
#       return(x)
#       
#     })
#     CO2Data <- do.call(rbind, Rawdata)
#     
#     #set limits (350-1000) for CA *baddata reports can be found in separet files
#     CO2Data$Concentration.1Min[CO2Data$Concentration.1Min > 1000] <- 1000
#     CO2Data$Concentration.1Min[CO2Data$Concentration.1Min < 350] <- 350
#     
#     CO2Data$DateTime <- nearestTimeStep(CO2Data$DateTime, 30, "ceiling")
#     
#     CO2Data$WindSpeed[CO2Data$WindSpeed < 0] <- NA
#     
#     CO2Data$IRGA.Pressure[CO2Data$IRGA.Pressure < 900] <- NA
#     CO2Data$IRGA.Pressure[CO2Data$IRGA.Pressure > 1100] <- NA
#     
#     sumCO2 <- data.table(CO2Data)[,list(CO2=mean(Concentration.1Min, na.rm=TRUE),
#                                         PRESS = mean(IRGA.Pressure, na.rm=TRUE),
#                                         WindSpeed=mean(WindSpeed, na.rm=TRUE),
#                                         PPFD=mean(PPFD, na.rm=TRUE),
#                                         Tair=mean(Air.Temp, na.rm=TRUE),
#                                         RH=mean(IRGA.Vapor.Pressure/saturate.vp.func(Air.Temp),na.rm=TRUE)
#                                         
#     ),
#     by = DateTime]
#     # give co2 data datetime
#     date.df <- data.frame(DateTime = 
#                             seq.POSIXt(as.POSIXct(startDate,tz = "UTC"),
#                                        as.POSIXct(paste0(startDate,' 23:30:00'),tz = "UTC"),
#                                        '30 min'))
#     
#     CaData.ls[[Ring]] <- merge(date.df, sumCO2, all=TRUE)
#     
#     CaData.ls[[Ring]]$Ring <- Ring
#   }
#   ca.df <- do.call(rbind,CaData.ls)
#   
#   # set ca for ambient and elevated rings
#   ca.df$treat <- NA
#   ca.df$treat[ca.df$Ring %in% c(2,3,6)] <- 'A'
#   ca.df$treat[ca.df$Ring %in% c(1,4,5)] <- 'E'
#   
#   ca.df <- ca.df[order(ca.df$DateTime),]
#   
#   ca.ppt.df <- merge(ca.df,ros15_30)
#   
#   ca.by.treat.df <- summaryBy(.~DateTime + treat,
#                               data = ca.ppt.df, FUN = median,keep.names = TRUE,
#                               na.rm=TRUE)
#   
#   tem.df <- spread(ca.by.treat.df,treat,CO2)
#   
#   met.df <- summaryBy(.~ DateTime,
#                       data = tem.df, 
#                       FUN = mean,
#                       keep.names = TRUE,
#                       na.rm=TRUE)
#   
#   # RH is 0-1
#   met.df$RH <- met.df$RH / 100
#   met.df$Tair[met.df$Tair > 55] <- 55
#   met.df$RH[met.df$RH > 1] <- 1
#   
#   # fill up nas with nearest data
#   # there should be na only when sensors in all rings are broken
#   na.col <- which(sapply(met.df, function(x) sum(is.na(x))) > 0)
#   
#   if (length(na.col) > 0){
#     for (i in seq_along(na.col)){
#       try(met.df[,na.col[i]] <- na.locf(met.df[,na.col[i]]))
#       try(met.df[,na.col[i]] <- na.locf(met.df[,na.col[i]],fromLast = T))
#     }
#   }
#   
#   # give met datatime
#   met.full.date <- data.frame(DateTime=seq.POSIXt(as.POSIXct(startDate,tz = "UTC"),
#                                                   as.POSIXct(paste(startDate,' 23:30:00'),tz = "UTC"), by="30 min"))
#   
#   met <- merge(met.full.date, met.df, all=TRUE)
#   
#   met <- subset(met,select = -c(Ring))
#   
#   names(met) <- c("DATE","PRESS","WIND","PAR",
#                   "TAIR","RH","PPT","Ca.A",'Ca.E')
#   
#   met$DATE <- as.Date(met$DATE)
#   
#   return(met)
# }
# 
# # function to upload file to hiev from Remko Duursma
# uploadToHIEv <- function(filename, experiment_id, type=c('RAW', 'UNKNOWN', 'CLEANSED', 'PROCESSED'),
#                          description=NULL, tag_names=NULL, parent_filenames=NULL, 
#                          start_time = '2018-03-01 00:00:00',
#                          end_time ='2018-03-01 00:00:00',
#                          label_names=NULL, 
#                          quiet=FALSE){
#   
#   
#   if(!file.exists(filename))stop("Provided filename does not exist.")
#   if(is.na(as.numeric(experiment_id)))stop("Experiment ID must be a number, find it on HIEv.")
#   type <- match.arg(type)
#   
#   if(!is.null(description) && length(description) > 1){
#     description <- paste(description, collapse="\r\n")
#   }
#   
#   if(is.null(getToken())){
#     # Attempt to set token by default location
#     setit <- setToken(interactive=FALSE)
#     if(!setit)stop("Set your HIEv API token with setToken() first.")
#   }
#   hievtoken <- get("hievtoken", envir=HIEv:::.hievenv)
#   
#   if(is.null(getRepository())){
#     setRepository()  # default value
#   }
#   repo <- getRepository()
#   url <- paste0(repo, "/data_files/api_create?auth_token=", hievtoken)
#   
#   postresult <- postForm(url,
#                          experiment_id=experiment_id,
#                          type=type,
#                          description=description,
#                          parent_filenames=parent_filenames,
#                          label_names=label_names,
#                          tag_names=tag_names,
#                          start_time = start_time,
#                          end_time = end_time,
#                          creator_email = 'jinyan.yang@uws.edu.au',
#                          file=fileUpload(filename), .checkParams=FALSE)
#   
#   postresult <- fromJSON(postresult)
#   
#   if(!quiet){
#     message("HIEv says:")
#     message(postresult$messages,"\n")
#     message("Filename on HIEv:")
#     message(postresult$file_name)
#   }
#   
#   return(invisible(postresult))
# }
# 
# # function to process the files
# make.met.func <- function(startDate,endDate){
#   
#   # 
#   met.df <- get.met.func(startDate,endDate)
#   
#   # 
#   out.fn <- sprintf('FACE_AUTO_RA_MET_L2_%s.dat',
#                     paste0(format(as.Date(startDate),'%Y%m%d'),'-',
#                            format(as.Date(endDate),'%Y%m%d')))
#   write.table(met.df,out.fn, row.names = F)
#   
#   uploadToHIEv(out.fn, experiment=21, type='PROCESSED',
#                description = c('Met data aggregated from ROS (rainfall seonsor) and',
#                                'EUcFACE (PAR (LI-190, Li-cor, Lincoln, NE, U.S.)',
#                                'wind speed (Wincap Ultrasonic WMT700 Vaisala, Vantaa, Finland)',
#                                'humidity, and temperature sensors (HUMICAP HMP 155 Vaisala, Vantaa, Finland))',
#                                " ",
#                                "Variables:",
#                                'DATE',"PRESS - air pressure - millibar",
#                                "WIND - windspeed - m/s",
#                                "PAR -photosynthetically active radiation - mumol m-2 s-1",
#                                "TAIR -air temperature - degree C",
#                                "RH - relative humidity - 0-1","PPT -precipitation - mm/d",
#                                "Ca.A - atmospheric CO2 convertation in ambient rings - mumol/mol",
#                                'Ca.E - atmospheric CO2 convertation in ambient rings - mumol/mol',
#                                " ",
#                                'Data used sources from ROS_WS_Table15 and FACE_AUTO_RingNumber_FCPLOGG.',
#                                'Original data were aggregated to mean of ambient and elevated groups every 30 min.'),
#                start_time = paste0(startDate,' 00:00:00'),
#                end_time = paste0(endDate,' 00:00:00'),
#                label_names= c('EucFACE'))
#   file.remove(out.fn)
# }
