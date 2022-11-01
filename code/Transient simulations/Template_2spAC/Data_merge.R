########################################################################
### R Script
###   to merge the 40 chuncks of data for each combination of species mass ratios
###   & merge all data for each bmr 
########################################################################
# Create a path directing where will be created all "BMRs" directory folders. Example: "./Template_2spAC";
Path <- "ZZ HERE ZZ";

# Create a directory where all the chunks of data are stored: e.g. "Databank" here.
dir1 <- paste(Path, "Databank",sep="/");
########################################################################
Fw <-'2spAC'             # Food web nomenclature, no invasion in Consumer C-resource A system
Gamma <- c(1, 2, 2.5, 3, 3.5, 3.75, 4, 5, 7.5, 10, 15, 20, 25, 50, 100);
####### 1. Merge the 40 subdata for each body mass ====
for(gamma in Gamma){
  setwd(paste(dir0, paste("BMRs", gamma, sep='_'), 'out', sep='/'));
  Sec_Name <- c("EQ", "ODE", "POP");
  
  data1 <- NULL; data2 <- NULL; data3 <- NULL; 
  for(i in 1:40){
    subdata1 <- read.table( paste(Fw, gamma, paste(i, Sec_Name[1], 'txt', sep='.'), sep='_'), h=T);
    subdata2 <- read.table( paste(Fw, gamma, paste(i, Sec_Name[2], 'txt', sep='.'), sep='_'), h=T);
    subdata3 <- read.table( paste(Fw, gamma, paste(i, Sec_Name[3], 'txt', sep='.'), sep='_'), h=T);
    
    data1 <- rbind(data1, subdata1);
    data2 <- rbind(data2, subdata2);
    data3 <- rbind(data3, subdata3);
  }
  setwd(dir1); # Go in the directory containing all data chunks
  write.table(data1, paste(paste(Fw, gamma, Sec_Name[1], sep='_'), 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
  write.table(data2, paste(paste(Fw, gamma, Sec_Name[2], sep='_'), 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
  write.table(data3, paste(paste(Fw, gamma, Sec_Name[3], sep='_'), 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
}
####### 2. Merge all body mass data together ====
# ODE
setwd(dir1);
Resdata <- data.frame();
for(gamma in Gamma){
  subdata <- read.table( paste(Fw, gamma, paste('ODE', 'txt', sep='.'), sep='_'), h=T)
  Resdata <- rbind(Resdata,subdata)
}
write.table(Resdata, '2spAC.txt', dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)

# POP
setwd(dir1);
Resdata <- data.frame();
for(gamma in Gamma){
  subdata <- read.table( paste(Fw, gamma, paste('POP', 'txt', sep='.'), sep='_'), h=T)
  Resdata <- rbind(Resdata,subdata)
}
write.table(Resdata, '2spAC_Pop.txt', dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)
