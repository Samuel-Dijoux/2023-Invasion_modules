###     Script Expl. Competition						###
###     Transient_EC.R											###
###     version with external functions			###
###     Basal species = A                   ###
###     Invading consumer = B               ###
###     Resident consumer = C               ###
################        Description
## (7) Food web nomenclatures used in our simulations:
# 2spAB and 2spAC for resident consumer-resource systems,
# and Fw ('TC', 'EC', 'AC', 'IGP') for different invasion type in consumer-resource systems
## For 16-25 species mass ratio and
## 80 200 combinations of abiotic conditions (401 T * 200 K) based on Sentis et al. (2017) E.Letters

## Please note that this script is an example to generate a subset (1/40) of the data for an unique species mass ratio.
## Visit the "Template_EC: folder to generate automatically all the scripts (R and Bash), 

################ Required packages and functions
library(deSolve)
library(rootSolve)

source("./Functions.R")
#################################### Script settings ====
Fw <-'EC'             # Food web nomenclature, invasion of a consumer species

Sim = 1;              # Rank order of the sub-simulation
Nsim = 40;            # Number of sub-simulation
Time = 5000;          # Time period of the simulation in year                     
Alpha = 1;           # Body mass ratio between basal and invading consumer
Beta = 1;            # Body mass ratio between invading (B) and resident (C) consumer
Gamma = Alpha*Beta;        # Body mass ratio between basal and resident consumer

Mx = 0.001;           # Mass of basal resource A
My = Mx*Alpha;        # Mass of invading consumer B
Mz = My*Beta;         # Mass of resident consumer C

Tmin <- c(0, seq(1.1, 39.1, by = 40/Nsim));
Tmax <- seq(1, 40, by=40/Nsim);
Tr<-10^(-12)          # Minimum biomass threshold under which species are considered extinct

#################################### Section EQ: Species biomass at equilibrium ====
PPMR_Name <- paste(Alpha, Beta, sep='.');
Sec_Name <- paste(Sim, 'EQ', sep='.');
Fw_Name <- paste(Fw, PPMR_Name, Sec_Name, sep='_');

S1 <- EQ_2spAC(Tmin[Sim], Tmax[Sim], Gamma=Gamma);
data <- rbind(S1);
data$Bx[data$Bx == "NaN"] <- 0;    ##remove NaN
data$Bz[data$Bz == "NaN"] <- 0;
data$Bx[data$Bx < 0]      <- 0;    ##remove negative biomasses
data$Bz[data$Bz < 0]      <- 0;
data<-as.data.frame(data);

write.table(data, paste(Fw_Name, 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#################################### Section ODE: System dynamic ====
Sec_Name <- paste(Sim, 'ODE', sep='.');
Fw_Name <- paste(Fw, PPMR_Name, Sec_Name, sep='_');

#### Creating the outcome data from the 80 200 simulations
outa <- mapply(ExpComp, Time=Time, data$Temp[1:dim(data)[1]], data$Car[1:dim(data)[1]],
		Bx=data$Bx[1:dim(data)[1]], By=1E-6, Bz=data$Bz[1:dim(data)[1]],
		data$Mx[1:dim(data)[1]], rep(My,dim(data)[1]), data$Mz[1:dim(data)[1]], 
		rep(Alpha, dim(data)[1]), rep(Beta,dim(data)[1]), data$Gamma[1:dim(data)[1]]);
outb <- data.frame(matrix(unlist(outa), nrow=157, ncol=dim(data)[1], byrow=F), stringsAsFactors=FALSE);
outc <- t(outb);
bigdata <- cbind.data.frame(data$Temp[1:dim(data)[1]], data$Car[1:dim(data)[1]], rep(Alpha, dim(data)[1]), rep(Beta,dim(data)[1]), data$Gamma[1:dim(data)[1]], outc);

#### Data Part One: System's state ====
EC <- as.data.frame(bigdata[,1:12]);
colnames(EC)<-c("Temperature", "Car", "Alpha", "Beta", "Gamma", "Time", "A", "B", "C", "T_extA", "T_extB", "T_extC");
## Data correction
EC$Afin <- (ifelse(EC$A>Tr,1,0));
EC$Bfin <- (ifelse(EC$B>Tr,1,0));
EC$Bfin[which(EC$Afin == 0)]  <- 0;
EC$Cfin <- (ifelse(EC$C>Tr,1,0));
EC$Cfin[which(EC$Afin == 0)]  <- 0;
EC$Tot  <- (EC$Afin+EC$Bfin+EC$Cfin);  

EC_cor <-cbind.data.frame( EC$Time, EC$Temperature, EC$Car,
                        EC$A, EC$B, EC$C, EC$Afin, EC$Bfin, EC$Cfin, EC$Tot,
                        EC$T_extA, EC$T_extB, EC$T_extC,
                        data$Mx , rep(My, dim(data)[1]), data$Mz,
                        EC$Alpha, EC$Beta, EC$Gamma);

colnames(EC_cor)<-c("Time", "Temperature", "Carr",
                      "A", "B", "C", "Afin", "Bfin", "Cfin", "Tot",
                      "ExtTimeA", "ExtTimeB", "ExtTimeC",
                      "M_A", "M_B", "M_C",
                      "Alpha", "Beta", "Gamma");

write.table(EC_cor, paste(Fw_Name, 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#################################### Section pop: metric ====
## Calculus of populations minima, mean, maxima, sd and cv
## for the last 10 years of the simulation
Sec_Name <- paste(Sim, 'POP', sep='.');
Fw_Name <- paste(Fw, PPMR_Name, Sec_Name, sep='_');

#### Data Part Two: System's metrics ====
sub_data <- data.frame();
for(i in 1:dim(bigdata)[1]){
  arrmat <- as.data.frame(matrix(unlist(bigdata[i,13:162]), ncol=15, byrow=T) );
  corrmat <- cbind( rep(bigdata[i,1],10), rep(bigdata[i,2], 10), seq(1, 10, by=1) , arrmat,
              rep(Mx, 10), rep(My, 10), rep(Mz, 10), rep(Alpha, 10), rep(Beta, 10), rep(Gamma, 10));
  sub_data <- rbind(sub_data, corrmat);
}
EC_pop <- as.data.frame(sub_data);

title<-c("Temperature", "Carr", "Year");
for(j in 1:3){
  title<-append(title,paste(LETTERS[j],c('Min', 'Max', 'Mean', 'Sd', 'CV'), sep='_'));
}
colnames(EC_pop) <- c(title, "M_A", "M_B", "M_C", "Alpha", "Beta", "Gamma");

write.table(EC_pop, paste(Fw_Name, 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)