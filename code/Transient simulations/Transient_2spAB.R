###     Script Consumer(sp. B)-Resource(sp. A)	###
###     Transient_2spAB.R												###
###     version with external functions					###
###     Basal species = A												###
###     Consumer      = B												###
###################################################
################        Description
## (7) Food web nomenclatures used in our simulations:
# 2spAB and 2spAC for resident consumer-resource systems,
# and Fw ('TC', 'EC', 'AC', 'IGP') for different invasion type in consumer-resource systems
## For 16-25 species mass ratio and
## 80 200 combinations of abiotic conditions (401 T * 200 K) based on Sentis et al. (2017) E.Letters

## Please note that this script is an example to generate a subset (1/40) of the data for an unique species mass ratio.
## Visit the "Template_2spAB: folder to generate automatically all the scripts (R and Bash), 

################ Required packages and functions
library(deSolve)
library(rootSolve)

source("./Functions.R")
#################################### Script settings ====
Fw <-'2spAB'          # Food web nomenclature, C-R system without invasion

Sim = 1;              # Rank order of the sub-simulation
Nsim = 40;            # Number of sub-simulation
Time = 5000;          # Time period of the simulation in year                     
Alpha = 1;            # Body mass ratio between basal and consumer species

Mx = 0.001;           # Mass of basal resource A
My = Mx*Alpha;        # Mass of consumer B

Tmin <- c(0, seq(1.1, 39.1, by = 40/Nsim));
Tmax <- seq(1, 40, by=40/Nsim);
Tr<-10^(-12)          # Minimum biomass threshold under which species are considered extinct

#################################### Section EQ: Species biomass at equilibrium ====
PPMR_Name <- Alpha
Sec_Name <- paste(Sim, 'EQ', sep='.')
Fw_Name <- paste(Fw, PPMR_Name, Sec_Name, sep='_')

S1 <- EQ_2spAB(Tmin[Sim],Tmax[Sim], Alpha=Alpha);

data <- rbind(S1);
data$Bx[data$Bx == "NaN"] <- 0;    ##remove NaN
data$By[data$By == "NaN"] <- 0;
data$Bx[data$Bx < 0]      <- 0;    ##remove negative biomasses
data$By[data$By < 0]      <- 0;
data<-as.data.frame(data);

write.table(data, paste(Fw_Name, 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#################################### Section ODE: System dynamic ====
Sec_Name <- paste(Sim, 'ODE', sep='.');
Fw_Name <- paste(Fw, PPMR_Name, Sec_Name, sep='_');

#### Creating the outcome data from the 80 200 simulations
outa <- mapply(TC_2spAB, Time=Time, data$Temp[1:dim(data)[1]], data$Car[1:dim(data)[1]],
		Bx=data$Bx[1:dim(data)[1]]*1.02, By=data$By[1:dim(data)[1]]*1.02,
		data$Mx[1:dim(data)[1]], data$My[1:dim(data)[1]], data$Alpha[1:dim(data)[1]]);
outb <- data.frame(matrix(unlist(outa), nrow=105, ncol=dim(data)[1], byrow=F), stringsAsFactors=FALSE);
outc <- t(outb);
bigdata <- as.data.frame( cbind(data$Temp[1:dim(data)[1]], data$Car[1:dim(data)[1]], data$Alpha[1:dim(data)[1]], outc));

#### Data Part One: System's state ====
bcAB <- as.data.frame(bigdata[,1:8]) 
colnames(bcAB)<-c("Temperature", "Car", "Alpha", "Time", "A", "B", "T_extA", "T_extB")

## Data correction
bcAB$Afin <- (ifelse(bcAB$A>Tr,1,0));
bcAB$Bfin <- (ifelse(bcAB$B>Tr,1,0));
bcAB$Bfin[which(bcAB$Afin == 0)]  <- 0;
bcAB$Tot  <- (bcAB$Afin+bcAB$Bfin);

bcAB_cor <-as.data.frame(cbind( bcAB$Time, bcAB$Temperature, bcAB$Car,
                           	bcAB$A, bcAB$B, bcAB$Afin, bcAB$Bfin, bcAB$Tot, 
				bcAB$T_extA, bcAB$T_extB,
			   	data$Mx , data$My, bcAB$Alpha));

colnames(bcAB_cor)<-c("Time", "Temperature", "Carr",
                "A", "B", "Afin", "Bfin", "Tot", "ExtTimeA", "ExtTimeB", "M_A", "M_B", "Alpha");

write.table(bcAB_cor, paste(Fw_Name, 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#################################### Section pop: metric ====
## Calculus of populations minima, mean, maxima, sd and cv
## for the last 10 years of the simulation
Sec_Name <- paste(Sim, 'POP', sep='.')
Fw_Name <- paste(Fw, PPMR_Name, Sec_Name, sep='_')

sub_data <- data.frame();
for(i in 1:dim(bigdata)[1]){
  arrmat <- as.data.frame(matrix(unlist(bigdata[i,9:108]), ncol=10, byrow=T) );
  corrmat <- cbind( rep(bigdata[i,1],10), rep(bigdata[i,2], 10), seq(1, 10, by=1) , arrmat,
              rep(Mx, 10), rep(My, 10), rep(Alpha, 10));
  sub_data <- rbind(sub_data, corrmat);
}
bcAB_pop <- as.data.frame(sub_data);

title <- c("Temperature", "Carr", "Year");
for(j in 1:2){
  title <- append(title, paste(LETTERS[j],c('Min', 'Max', 'Mean', 'Sd', 'CV'), sep='_'));
}
colnames(bcAB_pop) <- c(title, "M_A", "M_B", "Alpha");

write.table(bcAB_pop, paste(Fw_Name, 'txt', sep='.'), dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)