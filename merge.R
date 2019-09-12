# For General Purpose ##
## This file is designed to merge 
## MCPT simulation results from Cluster
## Author: Xueheng Shi ##
## Date: 11/01/2018 ##
## Version: V0
## Updated on 11/09/2018


library(data.table)
parent_path = getwd() 


## Get all subfolders where the results are stored to
subfolders = list.dirs(parent_path, recursive=T)[-1] 
folder.path = file.path(subfolders)


## Initialize empty data frame for storing merged results
#MCPT.bind = data.frame(matrix(0, nrow=5, ncol=22, byrow=T))

MCPT.bind = NULL

for (folder in folder.path){
  # obtain file.path
  file.path = paste(folder, list.files(folder, pattern="*.csv"), sep="/")
  data = read.csv(file.path, header=T,stringsAsFactors = F)
  MCPT.bind = rbind(MCPT.bind, data)
}

#View(MCPT.bind)

write.csv(MCPT.bind, "MCPT_Merged.csv",row.names=T) 
