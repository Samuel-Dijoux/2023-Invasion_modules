########################################################################
### R Script
###   to merge the 40 chuncks of data for each combination of species mass ratios
###   & merge all data for each bmr 
########################################################################
########################################################################
# Create a path directing where will be created all "BMRs" directory folders. Example: "./Template_2spAB";
Path <- "ZZ HERE ZZ";

# Create a directory where all the chunks of data are stored: e.g. "Databank" here.
dir1 <- paste(Path, "Databank",sep="/");
########################################################################
Fw <-'2spAB'             # Food web nomenclature, no invasion in Consumer B-resource A system
Alpha <- c(1, 2, 5, 10);# BMR B/A
####### 1. Merge the 40 subdata for each body mass ====
for(alpha in Alpha){
  setwd(paste(Path, paste("BMRs", alpha, sep='_'), 'out', sep='/'));
  Sec_Name <- c("EQ", "ODE", "POP");
  
  data1 <- NULL; data2 <- NULL; data3 <- NULL; 
  for(i in 1:40){
    subdata1 <- read.table( paste(Fw, alpha, paste(i, Sec_Name[1], 'txt', sep='.'), sep='_'), h=T);
    subdata2 <- read.table( paste(Fw, alpha, paste(i, Sec_Name[2], 'txt', sep='.'), sep='_'), h=T);
    
    data1 <- rbind(data1, subdata1);
    data2 <- rbind(data2, subdata2);
  }
  setwd(dir1); # Go in the directory containing all data chunks
  write.table(data1, paste(paste(Fw, alpha, Sec_Name[1], sep='_'), 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
  write.table(data2, paste(paste(Fw, alpha, Sec_Name[2], sep='_'), 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
}
####### 2. Merge all body mass data together ====
setwd(dir1);
Resdata <- data.frame();
for(alpha in Alpha){
  subdata <- read.table(paste(Fw, alpha, paste('ODE', 'txt', sep='.'), sep='_'), h=T)
  Resdata <- rbind(Resdata,subdata)
}
write.table(Resdata, '2spAB.txt', dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)