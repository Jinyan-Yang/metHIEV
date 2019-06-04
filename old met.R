for(year.num in c(2013)){
  for(mon.num in sprintf("%02d", 1:6)){
    
    startDate = paste0(year.num,'-',mon.num,'-01')
    startDate <- as.Date(startDate)
    endDate = format(startDate,paste0(year.num,'-',
                                      mon.num,
                                      sprintf('-%s',days_in_month(startDate))))
    endDate <- as.Date(endDate)
    
    make.met.func(startDate,endDate)
  }
}
