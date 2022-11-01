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
Fw <-'IGPB'             # Food web nomenclature, invasion of consumer B in Consumer C-resource A system
Alpha <- c(1, 2, 5, 10);# BMR B/A
Beta <- c(1, 2, 5, 10);# BMR C/B
####### 1. Merge the 40 subdata for each body mass ====
for(alpha in Alpha){
  for(beta in Beta){
    setwd(paste(Path, paste("BMRs", paste(alpha, beta, sep='.'), sep='_'), 'out', sep='/'));
    Sec_Name <- c("EQ", "ODE", "POP");
    
    data1 <- NULL; data2 <- NULL; data3 <- NULL; 
    for(i in 1:40){
      subdata1 <- read.table( paste(Fw, paste(alpha, beta, sep='.'), paste(i, Sec_Name[1], 'txt', sep='.'), sep='_'), h=T);
      subdata2 <- read.table( paste(Fw, paste(alpha, beta, sep='.'), paste(i, Sec_Name[2], 'txt', sep='.'), sep='_'), h=T);
      subdata3 <- read.table( paste(Fw, paste(alpha, beta, sep='.'), paste(i, Sec_Name[3], 'txt', sep='.'), sep='_'), h=T);
      
      data1 <- rbind(data1, subdata1);
      data2 <- rbind(data2, subdata2);
      data3 <- rbind(data3, subdata3);
    }
    setwd(dir1); # Go in the directory containing all data chunks
    write.table(data1, paste(paste(Fw, paste(alpha, beta, sep='.'), Sec_Name[1], sep='_'), 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
    write.table(data2, paste(paste(Fw, paste(alpha, beta, sep='.'), Sec_Name[2], sep='_'), 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
    write.table(data3, paste(paste(Fw, paste(alpha, beta, sep='.'), Sec_Name[3], sep='_'), 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
  }
}
####### 2. Merge all body mass data together ====
# ODE
setwd(dir1);
Invdata <- data.frame();
for(i in 1:length(Alpha)){
  for(j in 1:length(Beta)){
    gamma <- (Alpha[i] * Beta[j]);
    if(gamma<1){ next
    }else{
      subdata <- read.table(paste(paste(Fw, paste(Alpha[i], Beta[j], sep='.'), "ODE", sep='_'), 'txt', sep='.'), h=T,  sep = "\t");
      Invdata <- rbind(Invdata, subdata)
    }
  }
}
write.table(Invdata, paste(Fw, 'data', 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)

# POP
setwd(dir1);
Invdata <- data.frame();
for(i in 1:length(Alpha)){
  for(j in 1:length(Beta)){
    gamma <- (Alpha[i] * Beta[j]);
    if(gamma<1){ next
    }else{
      subdata <- read.table(paste(paste(Fw, paste(Alpha[i], Beta[j], sep='.'), "POP", sep='_'), 'txt', sep='.'), h=T,  sep = "\t");
      Invdata <- rbind(Invdata, subdata)
    }
  }
}
write.table(Invdata, "IGPB_Pop.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)