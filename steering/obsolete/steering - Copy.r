# R steering file for "rapid-climate-stress-testing" package, written by Donny Kim, San Diego State University, dkim8398@sdsu.edu (donnykim27@gmail.com)
# It calls Matlab/Octave scripts written by Keirnan Fowler, University of Melbourne, fowler.k@unimelb.edu.au
# Overall workflow or sequence of calling functions follows the example_run.m file written by Fowler.

# Integrated framework for rapid climate stress testing on a monthly timestep
# by Keirnan Fowler, Natasha Ballis, Avril Horne, Andrew John, Rory Nathan and Murray Peel 
# Licence: CC BY 3.0 - see https://creativecommons.org/licenses/by/3.0/au/


# ========================================================================================
### This R script is designed to govern the work flow of using modified version of Fowler's scripts (RCST).
# This R script will use system command to call the Octave scripts to run.

### Initial environment set up
# 1. Install GNU Octave (https://octave.org/download.html)
#    Once Octave is installed in your device, you need to add it to system PATH. I assume the OS to be Windows 10.
# 2. Identify the path to the folder containing the binary of Octave (e.g. "C:\Program Files\GNU Octave\Octave-7.1.0\mingw64\bin")
# 3. Multiple options exists to add new PATH. FYI: https://stackoverflow.com/questions/9546324/adding-a-directory-to-the-path-environment-variable-in-windows
#    HOWEVER, the most safest way would be using Windows GUI.
# 4. Click Windows "Start" button -> type "Edit the System Environment Variables" to search and open it -> check you are under the "Advanced" tab
#    -> click "Environment Variables" -> Either in "User Variables" or "System Variables" Click variable called “Path” and click “Edit…”
#    -> Click “New” -> Copy and paste the path to Octave binary folder (e.g. "C:\Program Files\GNU Octave\Octave-7.1.0\mingw64\bin").
# 5. Now we are ready.
# ========================================================================================


### ============================================================================
### Step 1: reformatting RF predictor for RCST inputs
### No need to run after processed in the first place
library(tidyverse)

# Set working directory
setwd('C:/Workspace/rapid-climate-stress-testing-main/steering/'); getwd();

# Source the functions
source('../R/utils.R')
# Note that directories for saving csv files are hard-coded in each function. To change the directory, edit "../R/utils.R".

# Creating inputs from RFpredictors compiled by SCCWRP (RB9)
RFpredictor_df <- read.csv('../RF_predictor/2022-03-28_RF_predictors_RB9_COMIDS.csv')[,1:45] # we only need column 1 to 45.
comid_list <- (RFpredictor_df %>% distinct(comid))$comid[1:5] # Trying on 5 comids for test purpose
RCST_inputs_from_compiled_RFpredictors(RFpredictor_df, comid_list) # Run it. Ignore the dir.create warning.# In testing stage, I limited it to only 5 COMIDS.

# Creating inputs from raw RFpredictors created by Ted Grantham
comid_list_raw <- list.files(path = '../CA_All_COMID_descriptors_Ted/', full.names = F, recursive = T, ignore.case = TRUE) %>% 
  sub(pattern = "(.*)\\..*$", replacement = "\\1")
dir.create('../hist_raw') # this returns a simple warning message, which is not an error.
write.csv(comid_list_raw, paste0('../hist_raw/comid_list_raw.csv'),  row.names = FALSE, quote=FALSE, col.names = FALSE)

# Loop it.
for (i in 1:length(comid_list_raw[1:5])){
  RCST_inputs_from_raw_RFpredictors(comid_list_raw[i])
}
### Step 1 end==================================================================




### ============================================================================
### Step 2: Calling Octave wrapper that governs the RCST code run.

### For the SCCWRP project, RCST code is modified from original version.
# Followings are the summary of modifications:
# 1. Wrapper function (run_wrapper.m): Overall process is reorganized to fit the purpose of SCCWRP RF project and built as a function. 
#    R can call this wrapper function and pass down the simple argument(s).
# 2. Any steps or components related to rainfall-runoff model are removed or inactive.
# 3. RCST saves final outputs (stochastically generated & perturbed climate variables) as CSV files (TS_P & TS_T) for each comid.
#    Modified scripts can process PET, if monthly PET input is available.
# 4. Additional Octave functions are implemented to handle the "table" datatype from MATLAB

### RCST code is complicated;
# Numerous small functions invoke the other functions in sequence.
# Therefore, it is not easy to keep track of the actual workflow by reading through the scripts.
# When making any edits, be extra cautious how one function can invoke/affect other functions (and their inputs, arguments, and outputs). 


# Set working directory
setwd('C:/Workspace/rapid-climate-stress-testing-main/steering/'); getwd();
library(processx)

#system(paste0("Octave example_run.m"),intern = T) # Most basic way to call the function.
#system(paste0("Octave run_wrapper_test.m comid_20332660"),intern = T) # Most basic way to call the function.
proc <- processx::process$new("Octave", c("./simple_test.m", "comid_20325361", "hist_raw", 111), stdout = "|")
proc <- processx::process$new("Octave", c("./simple_test.m", "comid_20325361"), stdout = "|")
proc$read_all_output()


# A little more advance way.
library(processx)

output <- character(0)
{
  proc <- processx::process$new("Octave", c("./run_wrapper.m", "comid_10000042", "hist_raw", 3), stdout = "|")
  while (proc$is_alive()) {
    Sys.sleep(1)
    now <- Sys.time()
    tmstmp <- sprintf("# [%s]", format(now, format = "%T"))
    thisout <- proc$read_output_lines()
    if (length(thisout)) {
      output <- c(output, thisout)
      message(tmstmp, " Msg from Ocatve:\n", paste0("#> ", thisout), "  ")
    }
  }
}


call_octave_scripts(comid_list[1])

# Make it into function
call_octave_scripts = function(comid_arg){
  require(processx)
  
  proc <- processx::process$new("Octave", c("./run_wrapper_test.m", paste0("comid_", comid_arg)), stdout = "|")
  
  #If you want real time messages...
  while (proc$is_alive()) {
    output <- character(0)
    Sys.sleep(1)
    now <- Sys.time()
    tmstmp <- sprintf("# [%s]", format(now, format = "%T"))
    thisout <- proc$read_output_lines()
    if (length(thisout)) {
      output <- c(output, thisout)
      message(tmstmp, " Msg from Ocatve:\n", paste0("#> ", thisout), "  ")
    }
  }
}

#proc1 <- processx::process$new("Octave", c("./run_wrapper_test.m", paste0("comid_", comid_list[1])), stdout = "|")
#proc2 <- processx::process$new("Octave", c("./run_wrapper_test.m", paste0("comid_", comid_list[2])), stdout = "|")
#proc3 <- processx::process$new("Octave", c("./run_wrapper_test.m", paste0("comid_", comid_list[3])), stdout = "|")
#proc4 <- processx::process$new("Octave", c("./run_wrapper_test.m", paste0("comid_", comid_list[4])), stdout = "|")
#proc5 <- processx::process$new("Octave", c("./run_wrapper_test.m", paste0("comid_", comid_list[5])), stdout = "|")


###### testing parallel run
library(parallel)
library(doParallel)
library(tictoc)
library(foreach)
cl <- parallel::makeCluster(5)
registerDoParallel(cl)

foreach(i = seq_along(comid_list), .combine = 'c') %dopar% {
  call_octave_scripts(comid_list[i])
}

stopCluster(cl) # VERY important to do this!!!!



#====test
cl <- parallel::makeCluster(3)
registerDoParallel(cl)

proc <- processx::process$new("Octave", c("./run_wrapper.m", "comid_10000042", "hist_raw", 1), stdout = "|")
proc2 <- processx::process$new("Octave", c("./run_wrapper.m", "comid_10000042", "hist_raw", 2), stdout = "|")
proc3 <- processx::process$new("Octave", c("./run_wrapper.m", "comid_10000042", "hist_raw", 3), stdout = "|")

stopCluster(cl) # VERY important to do this!!!!





# Alternative
# https://stackoverflow.com/questions/55067932/r-dopar-foreach-on-chunks-instead-of-per-line
if (!require("furrr")) install.packages("furrr"); library("furrr") # furrr package is essential.
temp = data.frame(comid_list_raw, rep('hist_raw', 5), rep(2, 5))
colnames(temp) = c('comid_arg', 'input_dir_arg', 'split_parm_arg')
future::plan(future::multisession(workers = 3))
furrr::future_pmap(temp[c('comid_arg', 'input_dir_arg', 'split_parm_arg')], Call_RCST_Octave_Scripts)


# Alternative 2
doParallel::registerDoParallel(cl = 2, cores = 2)

parallel::mclapply(
  comid_list_raw,
  Call_RCST_Octave_Scripts2)

parallel::stopCluster(cl) # VERY important to do this!!!!