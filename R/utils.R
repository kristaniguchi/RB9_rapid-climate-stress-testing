### Donny Kim.


### This is the function for creating climate input for RCST from compiled RF predictor csv files which was compiled by SCCWRP.
RCST_inputs_from_compiled_RFpredictors = function(predictor_input, comid_list_arg){
  
  # Need dplyr
  require(dplyr)
  
  comid_list <- comid_list_arg
  dir.create('../hist') # this returns a simple warning message, which is not an error.
  
  # Looping through comids, to write CSV files (P and T) for each comids
  for (i in 1:length(comid_list)){
    temp = predictor_input %>% filter(comid == comid_list[i])
    
    year_list <- c(); month_list <- sprintf('%02d', rep(1:12, nrow(temp)))
    P_list <- c(); T_list <- c();
    
    for (j in 1:nrow(temp)){
      year_list <-  append(year_list, rep(temp$wayr[j]-1, 12) ) # Mismatch in case. (raw = WaYr) but (compiled = wayr)
      P_list <- append(P_list, c(temp[j,16:24], temp[j, 13:15])) #pwy months and then wy months
      T_list <- append(T_list, c(temp[j,37:45], temp[j, 34:36])) #pwy months and then wy
    }
    year_list <-  append(year_list, rep(temp$wayr[length(temp$wayr)], 9) )
    month_list <- append(month_list, sprintf('%02d', rep(1:9)))
    P_list <- append(P_list, temp[length(temp$wayr),4:12])
    T_list <- append(T_list, temp[length(temp$wayr),25:33])
    
    precip_monthly = data.frame(matrix(ncol=3, nrow=9+nrow(temp)*12));   colnames(precip_monthly) = c("Year", "Month", comid_list[i]);
    T_monthly = data.frame(matrix(ncol=3, nrow=9+nrow(temp)*12));   colnames(T_monthly) = c("Year", "Month", comid_list[i]);
    
    precip_monthly$Year <- year_list; T_monthly$Year <- year_list
    precip_monthly$Month <- month_list; T_monthly$Month <- month_list
    precip_monthly[,3] <- as.numeric(P_list); T_monthly[,3] <- as.numeric(T_list)
    
    write.csv(precip_monthly, paste0('../hist/precip_monthly_', comid_list[i], ".csv"),  row.names = FALSE, quote=FALSE)
    write.csv(T_monthly, paste0('../hist/T_monthly_', comid_list[i], ".csv"),  row.names = FALSE, quote=FALSE)
  }
  comid_list <- paste0("comid_", comid_list)
  write.csv(comid_list, paste0('../hist/comid_list.csv'),  row.names = FALSE, quote=FALSE)
  
}




### This is the function for creating climate input for RCST from raw RF predictor csv files created by Ted Grantham.
RCST_inputs_from_raw_RFpredictors = function(comid_arg){
  
  # Need dplyr
  require(dplyr)
  
  # Read CSV file first using comid_arg
  predictor_input <- read.csv(paste0('../CA_All_COMID_descriptors_Ted/', comid_arg, '.csv'))[,1:45]
  
  temp = predictor_input %>% filter(COMID == comid_arg) # Not necessary but just in case. # Mismatch in case. (raw = COMID) but (compiled = comid)
  rm(predictor_input)
    
  year_list <- c(); month_list <- sprintf('%02d', rep(1:12, nrow(temp)))
  P_list <- c(); T_list <- c();
    
  for (j in 1:nrow(temp)){
    year_list <-  append(year_list, rep(temp$WaYr[j]-1, 12) ) # Mismatch in case. (raw = WaYr) but (compiled = wayr)
    P_list <- append(P_list, c(temp[j,16:24], temp[j, 13:15])) #pwy months and then wy
    T_list <- append(T_list, c(temp[j,37:45], temp[j, 34:36])) #pwy months and then wy
  }
  # appending the "last WY's data"
  year_list <-  append(year_list, rep(temp$WaYr[length(temp$WaYr)], 9) ) # Mismatch in case. (raw = WaYr) but (compiled = wayr)
  month_list <- append(month_list, sprintf('%02d', rep(1:9)))
  P_list <- append(P_list, temp[length(temp$WaYr),4:12]) #pwy months and then wy
  T_list <- append(T_list, temp[length(temp$WaYr),25:33]) #pwy months and then wy
  
  precip_monthly = data.frame(matrix(ncol=3, nrow=9+nrow(temp)*12));   colnames(precip_monthly) = c("Year", "Month", comid_arg);
  T_monthly = data.frame(matrix(ncol=3, nrow=9+nrow(temp)*12));   colnames(T_monthly) = c("Year", "Month", comid_arg);
  
  precip_monthly$Year <- year_list; T_monthly$Year <- year_list
  precip_monthly$Month <- month_list; T_monthly$Month <- month_list
  precip_monthly[,3] <- as.numeric(P_list); T_monthly[,3] <- as.numeric(T_list)
    
  
  write.csv(precip_monthly, paste0('../hist_raw/precip_monthly_', comid_arg, ".csv"),  row.names = FALSE, quote=FALSE)
  write.csv(T_monthly, paste0('../hist_raw/T_monthly_', comid_arg, ".csv"),  row.names = FALSE, quote=FALSE)
}