########################################################################
### R Script                      Tot.data
###
###   Quantifying the influence of species invasion on resident communities:
###   changes in local diversity, species composition & dynamic state
########################################################################
library(deSolve)
library(rootSolve)
library(lattice)
library(gridExtra)
library(RColorBrewer)

source("./Functions.R");
############# 2spAB; C-R system without invasion ====
Res <- read.table("./data/2spAB.txt", h=T,  sep = "\t");       # Data from 2sp C-R (A, B) simulation
Alpha <- c(1, 2, 5, 10);
Tr <- 1e-12;

#### Initiating the dataset ====
Res_BMR <- rep('Alpha', length=dim(Res)[1]);
Rr_i <- rep(0, length=dim(Res)[1]);         # Initial biomass of resource at equilibirum
Cr_i <- rep(0, length=dim(Res)[1]);         # Initial biomass of consumer at equilibirum
Dyn_state_i <- rep(0, length=dim(Res)[1]);  # Dynamic state initial of resident species
Dyn_state_f <- rep(0, length=dim(Res)[1]);  # Dynamic state final, after 5000 years
Zone_i <- c(rep(0, length=dim(Res)[1]));    # Attributed Code given Dyn_state_i
Zone_f <- c(rep(0, length=dim(Res)[1]));    # Attributed Code given Dyn_state_f

Rr_eq <- rep(0, length=dim(Res)[1]);      # Biomass density of resident resource at equilibrium
Cr_eq <- rep(0, length=dim(Res)[1]);      # Biomass density of resident consumer at equilibrium
Tot_eq <- rep(0, length=dim(Res)[1]);     # Total diversity at equilibrium
Dom_Eigenvalue <- rep(NA, length=dim(Res)[1]); # Dominant eigenvalues for the different community structure

EQ_transit_InitFin <- rep(0, length=dim(Res)[1]); # Change in the system stability between the initial and final outcomes
Zone_EQtransit_InitFin <- rep(0, length=dim(Res)[1]); # Corresponding codes of EQ_transit_if

CR_Tot_Alpha <- cbind.data.frame(
  Res$Temperature, Res$Carr, Res_BMR, Res$Alpha, Res$M_A, Res$M_B,
  Rr_i, Cr_i, Res$A, Res$B, Res$Afin, Res$Bfin, Res$Tot, 
  Res$ExtTimeA, Res$ExtTimeB, Res$Time, 
  Rr_eq, Cr_eq, Dom_Eigenvalue, Tot_eq, Dyn_state_i, Zone_i, Dyn_state_f, Zone_f, EQ_transit_InitFin, Zone_EQtransit_InitFin);

colnames(CR_Tot_Alpha)<-c("Temp", "Car", "Res_BMR", "Alpha", "M_Rr", "M_Cr",
                         "Rr_i", "Cr_i", "Rr_f", "Cr_f", "Rr_fin", "Cr_fin", "Tot",
                         "Rr_extime", "Cr_extime", "Sim_duration", 
                         "Rr_eq", "Cr_eq", "Dom_Eigenvalue","Tot_eq",
                         "Dyn_state_i", "Zone_i", "Dyn_state_f", "Zone_f", "EQ_transit_InitFin", "Zone_EQtransit_InitFin");
#### Initial state ====
## Consumer-Resource at equilibrium
dataieq <- as.data.frame(EQUI_2spAB(Temp=CR_Tot_Alpha$Temp, Car=CR_Tot_Alpha$Car, Alpha=CR_Tot_Alpha$Alpha))
# correcting negative biomass
dataieq$By <- ifelse(dataieq$By<Tr,0,dataieq$By)
dataieq$Bx <- ifelse(dataieq$Bx<Tr,0,dataieq$Bx);
dataieq$By <- ifelse(dataieq$Bx<Tr,0,dataieq$By);

CR_Tot_Alpha$Rr_i <- dataieq$Bx;
CR_Tot_Alpha$Cr_i <- dataieq$By;
CR_Tot_Alpha$Rr_i <- ifelse(CR_Tot_Alpha$Rr_i<Tr,1e-6,CR_Tot_Alpha$Rr_i);
CR_Tot_Alpha$Cr_i <- ifelse(CR_Tot_Alpha$Cr_i<Tr,1e-6,CR_Tot_Alpha$Cr_i);

# Calculus of Eigenvalues
data_2spi <- dataieq[which(CR_Tot_Alpha$Cr_i>Tr & CR_Tot_Alpha$Cr_i != 1e-6),];
Dyn <- rep(0, length=dim(data_2spi)[1]); Zone <- rep(0, length=dim(data_2spi)[1]);

lambda_i <- mapply(lambda_2sp_AB, Temp=data_2spi$Temp[1:dim(data_2spi)[1]], Car=data_2spi$Car[1:dim(data_2spi)[1]],
                  Bx=data_2spi$Bx[1:dim(data_2spi)[1]], By=data_2spi$By[1:dim(data_2spi)[1]],
                  Mx=data_2spi$Mx[1:dim(data_2spi)[1]], My=data_2spi$My[1:dim(data_2spi)[1]], Alpha=data_2spi$Alpha[1:dim(data_2spi)[1]]);

data_2spi <- cbind.data.frame(data_2spi, lambda_i, Dyn, Zone);
data_2spi$Dyn[which(data_2spi$lambda_i > 0)] <- paste("Cr", "Rr", "cycles", sep='.');
data_2spi$Dyn[which(data_2spi$lambda_i == 0 | data_2spi$lambda_i < 0)] <- paste("Cr", "Rr", "Eq", sep='.');
data_2spi$Zone[data_2spi$Dyn == paste("Cr", "Rr", "cycles", sep='.')] <- 3;
data_2spi$Zone[data_2spi$Dyn == paste("Cr", "Rr", "Eq", sep='.')] <- 2;

CR_Tot_Alpha$Dyn_state_i[which(CR_Tot_Alpha$Cr_i>Tr & CR_Tot_Alpha$Cr_i != 1e-6)] <- data_2spi$Dyn;
CR_Tot_Alpha$Zone_i[which(CR_Tot_Alpha$Cr_i>Tr & CR_Tot_Alpha$Cr_i != 1e-6)] <- data_2spi$Zone;
## Only R remain at equilibrium
data_1spi <- CR_Tot_Alpha[which(CR_Tot_Alpha$Cr_i == 1e-6),]
data_1spi$Rr_i <- K(Temp=data_1spi$Temp, Mx=data_1spi$M_Rr, Car=data_1spi$Car)

CR_Tot_Alpha$Rr_i[which(CR_Tot_Alpha$Cr_i == 1e-6)] <- data_1spi$Rr_i
CR_Tot_Alpha$Dyn_state_i[which(CR_Tot_Alpha$Cr_i == 1e-6)] <- paste("R","eq",sep='.')
CR_Tot_Alpha$Zone_i[which(CR_Tot_Alpha$Cr_i == 1e-6)] <- 1;
#### Final state ====
data_2spf <- CR_Tot_Alpha[which(CR_Tot_Alpha$Tot==2),]
## Consumer-Resource at equilibrium
datafeq <- as.data.frame(EQUI_2spAB(Temp=data_2spf$Temp, Car=data_2spf$Car, Alpha=data_2spf$Alpha))
# correcting negative biomass
datafeq$By <- ifelse(datafeq$By<Tr,0,datafeq$By)
datafeq$Bx <- ifelse(datafeq$Bx<Tr,0,datafeq$Bx);
datafeq$By <- ifelse(datafeq$Bx<Tr,0,datafeq$By);
Dyn <- rep(0, length=dim(datafeq)[1]); Zone <- rep(0, length=dim(datafeq)[1]);

# Calculus of Eigenvalues
lambda_f <- mapply(lambda_2sp_AB, Temp=datafeq$Temp[1:dim(datafeq)[1]], Car=datafeq$Car[1:dim(datafeq)[1]],
                   Bx=datafeq$Bx[1:dim(datafeq)[1]], By=datafeq$By[1:dim(datafeq)[1]],
                   Mx=datafeq$Mx[1:dim(datafeq)[1]], My=datafeq$My[1:dim(datafeq)[1]], Alpha=datafeq$Alpha[1:dim(datafeq)[1]]);
datafeq <- cbind.data.frame(datafeq, lambda_f, Dyn, Zone);

datafeq$Dyn[which(datafeq$lambda_f > 0)] <- paste("Cr", "Rr", "cycles", sep='.');
datafeq$Dyn[which(datafeq$lambda_f == 0 | datafeq$lambda_f < 0)] <- paste("Cr", "Rr", "eq", sep='.');
datafeq$Zone[which(datafeq$Dyn == paste("Cr", "Rr", "cycles", sep='.'))] <- 3;
datafeq$Zone[which(datafeq$Dyn == paste("Cr", "Rr", "eq", sep='.'))] <- 2;
  
CR_Tot_Alpha$Rr_eq[which(CR_Tot_Alpha$Tot==2)] <- datafeq$Bx;
CR_Tot_Alpha$Cr_eq[which(CR_Tot_Alpha$Tot==2)] <- datafeq$By;
CR_Tot_Alpha$Tot_eq[which(CR_Tot_Alpha$Tot==2 & CR_Tot_Alpha$Cr_eq > Tr)] <- 2;
CR_Tot_Alpha$Tot_eq[which(CR_Tot_Alpha$Tot==2 & CR_Tot_Alpha$Cr_eq < Tr)] <- 1;
CR_Tot_Alpha$Dyn_state_f[which(CR_Tot_Alpha$Tot==2)] <- datafeq$Dyn;
CR_Tot_Alpha$Zone_f[which(CR_Tot_Alpha$Tot==2)] <- datafeq$Zone;

## Only R remain at equilibrium
CR_Tot_Alpha$Tot_eq[which(CR_Tot_Alpha$Tot == 1 & CR_Tot_Alpha$Rr_fin == 1 )] <- 1;
CR_Tot_Alpha$Dyn_state_f[which(CR_Tot_Alpha$Tot == 1 & CR_Tot_Alpha$Rr_fin == 1 )] <- paste("R", "eq", sep=".");
CR_Tot_Alpha$Zone_f[which(CR_Tot_Alpha$Tot == 1 & CR_Tot_Alpha$Rr_fin == 1 )] <- 1;

# No species remaining
CR_Tot_Alpha$Dyn_state_f[which(CR_Tot_Alpha$Tot == 0)] <- paste(0, "species", sep=".");
CR_Tot_Alpha$Zone_f[which(CR_Tot_Alpha $Tot == 0)] <- 0;

#### Changes in system stability in the community (initially vs. after 5000 years) ====
# 0. Initial-Final configurations
conf_name <- c(); initeq <- c("Eq", "Cycles"); fineq <- c("Null", "Eq", "Cycles")
for(i in initeq){conf_name <- append(conf_name, paste(i, fineq, sep='.'))};
conf_code <- c(0:5);
Conf_tab <- cbind.data.frame(conf_name, conf_code);

summary(as.factor(CR_Tot_Alpha$Zone_i))
summary(as.factor(CR_Tot_Alpha$Zone_f))

# A. Prepare the different combinations of stability changes
init_eq_code <- c(1, 2, 3);
init_eq <- c(rep("Eq", 2), "Cycles");
final_eq_code <- c(0:3);
final_eq <- c("Null", rep("Eq", 2), "Cycles");

# Data of observe combinations
I_code <- c(); for(i in init_eq_code){I_code <- append(I_code, rep(i, length(final_eq)))};
F_code <- rep(final_eq_code, length(init_eq_code));

IF_eq <- c(); for(i in init_eq){IF_eq <- append(IF_eq, paste(i, final_eq, sep='.'))}
IF_code <- rep(0, length(IF_eq));
IF <- cbind.data.frame(IF_eq, IF_code, I_code, F_code);

IF$IF_eq <- as.character(IF$IF_eq);
for(i in 1:length(IF_eq)){ IF$IF_code[i] <- Conf_tab$conf_code[which(Conf_tab$conf_name==IF_eq[i])]};

# B. Attribute the AC_Tot_data's outcomes to the corresponding configurations EQ code and name
CR_Tot_Alpha$Zone_i <- as.numeric(CR_Tot_Alpha$Zone_i);
CR_Tot_Alpha$Zone_f <- as.numeric(CR_Tot_Alpha$Zone_f);
CR_Tot_Alpha$EQ_transit_InitFin <- as.character(CR_Tot_Alpha$EQ_transit_InitFin)
CR_Tot_Alpha$Zone_EQtransit_InitFin <- as.numeric(CR_Tot_Alpha$Zone_EQtransit_InitFin);
for(i in 1:dim(IF)[1]){
  CR_Tot_Alpha$EQ_transit_InitFin[which(CR_Tot_Alpha$Zone_i == IF$I_code[i] & CR_Tot_Alpha$Zone_f==IF$F_code[i])] <- IF$IF_eq[i];
  CR_Tot_Alpha$Zone_EQtransit_InitFin[which(CR_Tot_Alpha$Zone_i == IF$I_code[i] & CR_Tot_Alpha$Zone_f==IF$F_code[i])] <- IF$IF_code[i]
}
#### Saving the dataset ====
write.table(CR_Tot_Alpha, "./data/CR_Tot_Alpha.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
############# 2spAC; C-R system without invasion ====
Res <- read.table("2spAC.txt", h=T,  sep = "\t");       # Data from 2sp C-R (A, C) simulation
Gamma <- c(1, 1.5, 2, 2.5, 3, 3.75, 4, 5, 7.5, 10, 15, 20, 25, 50, 100);
Tr <- 1e-12;

#### Initiating the dataset ====
Res_BMR <- rep('Gamma', length=dim(Res)[1]);
Rr_i <- rep(0, length=dim(Res)[1]);         # Initial biomass of resource at equilibirum
Cr_i <- rep(0, length=dim(Res)[1]);         # Initial biomass of consumer at equilibirum
Dyn_state_i <- rep(0, length=dim(Res)[1]);  # Dynamic state initial of resident species
Dyn_state_f <- rep(0, length=dim(Res)[1]);  # Dynamic state final, after 5000 years
Zone_i <- c(rep(0, length=dim(Res)[1]));    # Attributed Code given Dyn_state_i
Zone_f <- c(rep(0, length=dim(Res)[1]));    # Attributed Code given Dyn_state_f

Rr_eq <- rep(0, length=dim(Res)[1]);      # Biomass density of resident resource at equilibrium
Cr_eq <- rep(0, length=dim(Res)[1]);      # Biomass density of resident consumer at equilibrium
Tot_eq <- rep(0, length=dim(Res)[1]);     # Total diversity at equilibrium
Dom_Eigenvalue <- rep(NA, length=dim(Res)[1]); # Dominant eigenvalues for the different community structure

EQ_transit_InitFin <- rep(0, length=dim(Res)[1]); # Change in the system stability between the initial and final outcomes
Zone_EQtransit_InitFin <- rep(0, length=dim(Res)[1]); # Corresponding codes of EQ_transit_if

CR_Tot_Gamma <- cbind.data.frame(
  Res$Temperature, Res$Carr, Res_BMR, Res$Gamma, Res$M_A, Res$M_C,
  Rr_i, Cr_i, Res$A, Res$C, Res$Afin, Res$Cfin, Res$Tot, 
  Res$ExtTimeA, Res$ExtTimeC, Res$Time, 
  Rr_eq, Cr_eq, Dom_Eigenvalue, Tot_eq, Dyn_state_i, Dyn_state_f, Zone_i, Zone_f, EQ_transit_InitFin, Zone_EQtransit_InitFin);

colnames(CR_Tot_Gamma)<-c("Temp", "Car", "Res_BMR", "Gamma", "M_Rr", "M_Cr",
                          "Rr_i", "Cr_i", "Rr_f", "Cr_f", "Rr_fin", "Cr_fin", "Tot",
                          "Rr_extime", "Cr_extime", "Sim_duration", 
                          "Rr_eq", "Cr_eq", "Dom_Eigenvalue","Tot_eq",
                          "Dyn_state_i", "Dyn_state_f", "Zone_i", "Zone_f", "EQ_transit_InitFin", "Zone_EQtransit_InitFin");

#### Initial state ====
## Consumer-Resource at equilibrium
dataieq <- as.data.frame(EQUI_2spAC(Temp=CR_Tot_Gamma$Temp, Car=CR_Tot_Gamma$Car, Gamma=CR_Tot_Gamma$Gamma))
# correcting negative biomass
dataieq$Bz <- ifelse(dataieq$Bz<Tr,0,dataieq$Bz);
dataieq$Bx <- ifelse(dataieq$Bx<Tr,0,dataieq$Bx);
dataieq$Bz <- ifelse(dataieq$Bx<Tr,0,dataieq$Bz);

CR_Tot_Gamma$Rr_i <- dataieq$Bx;
CR_Tot_Gamma$Cr_i <- dataieq$Bz;
CR_Tot_Gamma$Rr_i <- ifelse(CR_Tot_Gamma$Rr_i<Tr,1e-6,CR_Tot_Gamma$Rr_i);
CR_Tot_Gamma$Cr_i <- ifelse(CR_Tot_Gamma$Cr_i<Tr,1e-6,CR_Tot_Gamma$Cr_i);

# Calculus of Eigenvalues
data_2spi <- dataieq[which(dataieq$Bz>Tr),];
Dyn <- rep(0, length=dim(data_2spi)[1]); Zone <- rep(0, length=dim(data_2spi)[1]);

lambda_i <- mapply(lambda_2sp_AC, Temp=data_2spi$Temp[1:dim(data_2spi)[1]], Car=data_2spi$Car[1:dim(data_2spi)[1]],
                   Bx=data_2spi$Bx[1:dim(data_2spi)[1]], Bz=data_2spi$Bz[1:dim(data_2spi)[1]],
                   Mx=data_2spi$Mx[1:dim(data_2spi)[1]], Mz=data_2spi$Mz[1:dim(data_2spi)[1]], Gamma=data_2spi$Gamma[1:dim(data_2spi)[1]]);

data_2spi <- cbind.data.frame(data_2spi, lambda_i, Dyn, Zone);
data_2spi$Dyn <- ifelse(data_2spi$lambda_i > 0, paste("Cr", "Rr", "cycles", sep='.'), paste("Cr", "Rr", "Eq", sep='.'));
data_2spi$Zone <- ifelse(data_2spi$Dyn == paste("Cr", "Rr", "cycles", sep='.'), 3, 2);

CR_Tot_Gamma$Dyn_state_i[which(CR_Tot_Gamma$Cr_i>Tr & CR_Tot_Gamma$Cr_i != 1e-6)] <- data_2spi$Dyn;
CR_Tot_Gamma$Zone_i[which(CR_Tot_Gamma$Cr_i>Tr & CR_Tot_Gamma$Cr_i != 1e-6)] <- data_2spi$Zone;
## Only R remain at equilibrium
data_1spi <- CR_Tot_Gamma[which(CR_Tot_Gamma$Cr_i == 1e-6),]
data_1spi$Rr_i <- K(Temp=data_1spi$Temp, Mx=data_1spi$M_Rr, Car=data_1spi$Car)

CR_Tot_Gamma$Rr_i[which(CR_Tot_Gamma$Cr_i == 1e-6)] <- data_1spi$Rr_i
CR_Tot_Gamma$Dyn_state_i[which(CR_Tot_Gamma$Cr_i == 1e-6)] <- paste("R","eq",sep='.')
CR_Tot_Gamma$Zone_i[which(CR_Tot_Gamma$Cr_i == 1e-6)] <- 1;
#### Final state ====
data_2spf <- CR_Tot_Gamma[which(CR_Tot_Gamma$Tot==2),]
## Consumer-Resource at equilibrium
datafeq <- as.data.frame(EQUI_2spAC(Temp=data_2spf$Temp, Car=data_2spf$Car, Gamma=data_2spf$Gamma))
# correcting negative biomass
datafeq$Bz <- ifelse(datafeq$Bz<Tr,0,datafeq$Bz)
datafeq$Bx <- ifelse(datafeq$Bx<Tr,0,datafeq$Bx);
datafeq$Bz <- ifelse(datafeq$Bx<Tr,0,datafeq$Bz);
Dyn <- rep(0, length=dim(datafeq)[1]); Zone <- rep(0, length=dim(datafeq)[1]);

# Calculus of Eigenvalues
lambda_f <- mapply(lambda_2sp_AC, Temp=datafeq$Temp[1:dim(datafeq)[1]], Car=datafeq$Car[1:dim(datafeq)[1]],
                   Bx=datafeq$Bx[1:dim(datafeq)[1]], Bz=datafeq$Bz[1:dim(datafeq)[1]],
                   Mx=datafeq$Mx[1:dim(datafeq)[1]], Mz=datafeq$Mz[1:dim(datafeq)[1]], Gamma=datafeq$Gamma[1:dim(datafeq)[1]]);
datafeq <- cbind.data.frame(datafeq, lambda_f, Dyn, Zone);

datafeq$Dyn <- ifelse(datafeq$lambda_f > 0, paste("Cr", "Rr", "cycles", sep='.'), paste("Cr", "Rr", "eq", sep='.'));
datafeq$Zone <- ifelse(datafeq$Dyn == paste("Cr", "Rr", "cycles", sep='.'), 3, 2);

CR_Tot_Gamma$Rr_eq[which(CR_Tot_Gamma$Tot==2)] <- datafeq$Bx;
CR_Tot_Gamma$Cr_eq[which(CR_Tot_Gamma$Tot==2)] <- datafeq$Bz;
CR_Tot_Gamma$Tot_eq[which(CR_Tot_Gamma$Tot==2)] <- 2;
CR_Tot_Gamma$Dyn_state_f[which(CR_Tot_Gamma$Tot==2)] <- datafeq$Dyn;
CR_Tot_Gamma$Zone_f[which(CR_Tot_Gamma$Tot==2)] <- datafeq$Zone;

## Only R remain at equilibrium
CR_Tot_Gamma$Dyn_state_f[which(CR_Tot_Gamma$Tot == 1 & CR_Tot_Gamma$Rr_fin == 1 )] <- paste("R", "eq", sep=".");
CR_Tot_Gamma$Zone_f[which(CR_Tot_Gamma$Tot == 1)] <- 1;
CR_Tot_Gamma$Tot_eq[which(CR_Tot_Gamma$Tot==1)] <- 1;

# No species remaining
CR_Tot_Gamma$Dyn_state_f[which(CR_Tot_Gamma$Tot == 0)] <- "Null";
CR_Tot_Gamma$Zone_f[which(CR_Tot_Gamma $Tot == 0)] <- 0;
#### Changes in system stability in the community (initially vs. after 5000 years) ====
# 0. Initial-Final configurations
conf_name <- c(); initeq <- c("Eq", "Cycles"); fineq <- c("Null", "Eq", "Cycles")
for(i in initeq){conf_name <- append(conf_name, paste(i, fineq, sep='.'))};
conf_code <- c(0:5);
Conf_tab <- cbind.data.frame(conf_name, conf_code);

summary(as.factor(CR_Tot_Gamma$Zone_i))
summary(as.factor(CR_Tot_Gamma$Zone_f))

# A. Prepare the different combinations of stability changes
init_eq_code <- c(1, 2, 3);
init_eq <- c(rep("Eq", 2), "Cycles");
final_eq_code <- c(0:3);
final_eq <- c("Null", rep("Eq", 2), "Cycles");

# Data of observe combinations
I_code <- c(); for(i in init_eq_code){I_code <- append(I_code, rep(i, length(final_eq)))};
F_code <- rep(final_eq_code, length(init_eq_code));

IF_eq <- c(); for(i in init_eq){IF_eq <- append(IF_eq, paste(i, final_eq, sep='.'))}
IF_code <- rep(0, length(IF_eq));
IF <- cbind.data.frame(IF_eq, IF_code, I_code, F_code);

IF$IF_eq <- as.character(IF$IF_eq);
for(i in 1:length(IF_eq)){ IF$IF_code[i] <- Conf_tab$conf_code[which(Conf_tab$conf_name==IF_eq[i])]};

# B. Attribute the AC_Tot_data's outcomes to the corresponding configurations EQ code and name
CR_Tot_Gamma$Zone_i <- as.numeric(CR_Tot_Gamma$Zone_i);
CR_Tot_Gamma$Zone_f <- as.numeric(CR_Tot_Gamma$Zone_f);
CR_Tot_Gamma$EQ_transit_InitFin <- as.character(CR_Tot_Gamma$EQ_transit_InitFin)
CR_Tot_Gamma$Zone_EQtransit_InitFin <- as.numeric(CR_Tot_Gamma$Zone_EQtransit_InitFin);
for(i in 1:dim(IF)[1]){
  CR_Tot_Gamma$EQ_transit_InitFin[which(CR_Tot_Gamma$Zone_i == IF$I_code[i] & CR_Tot_Gamma$Zone_f==IF$F_code[i])] <- IF$IF_eq[i];
  CR_Tot_Gamma$Zone_EQtransit_InitFin[which(CR_Tot_Gamma$Zone_i == IF$I_code[i] & CR_Tot_Gamma$Zone_f==IF$F_code[i])] <- IF$IF_code[i]
}
#### Saving the dataset ====
write.table(CR_Tot_Gamma, "./data/CR_Tot_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
############# AC; Apparent competition ====
Inv <- read.table("./data/AC.data.txt", h=T,  sep = "\t");     # Data from AC module simulation
Res <- read.table("./data/CR_Tot_Gamma.txt", h=T,  sep = "\t");# Complete Data from 2sp C-R (A, C) including simulation outcomes and equilibrium analyses

Alpha <- c(0.1, 0.2, 0.5, 0.75, 1, 1.5, 2, 5, 10);# BMR B/A
Beta <- c(1, 2, 5, 10);                           # BMR C/B
Gamma <- c(1, 1.5, 2, 2.5, 3, 3.75, 4, 5, 7.5, 10, 15, 20, 25, 50, 100); 
Tr <- 1e-12;      # Extinction threshold

#### Initiating the dataset ====
InvTP <- rep(1, dim(Inv)[1]); # Trophic position of the invading species (here 1 = basal resource)
Fw <- rep("AC", dim(Inv)[1]); # Food web module in which invasion occurs (here AC = Apparent competition)
Res_BMR <- rep('Gamma', length=dim(Inv)[1]);# species bmr of resident system

Dyn_state_i <- rep('NA', length=dim(Inv)[1]);# Dynamic state initial of resident species
Dyn_state_f <- rep('NA', length=dim(Inv)[1]);# Dynamic state final, after invasion
Zone_i <- rep(0, length=dim(Inv)[1]);     # initial, identical to 2spAC
Zone_f <- rep(0, length=dim(Inv)[1]);     # final
BioD <- rep(0, length=dim(Inv)[1]);       # Biodiversity change after invasion
State <- rep(0, length=dim(Inv)[1]);      # Invasion states
Rr_i <- rep(0, length=dim(Inv)[1]);       # Initial biomass density of resident resource (starting at C-R eq)
Cr_i <- rep(0, length=dim(Inv)[1]);       # Initial biomass density of resident consumer (starting at C-R eq)
Ri_i <- rep(1E-6, length=dim(Inv)[1]);    # Initial biomass density of invading resource

Rr_eq <- rep(0, length=dim(Inv)[1]);      # Biomass density of resident resource at equilibrium
Ri_eq <- rep(0, length=dim(Inv)[1]);      # Biomass density of invading resource at equilibrium
Cr_eq <- rep(0, length=dim(Inv)[1]);      # Biomass density of resident consumer at equilibrium
Tot_eq <- rep(0, length=dim(Inv)[1]);     # Total diversity at equilibrium
Dom_Eigenvalue <- rep(NA, length=dim(Inv)[1]); # Dominant eigenvalues for the different community structure

EQ_transit_InitFin <- rep(0, length=dim(Inv)[1]); # Change in the system stability between the initial and final outcomes
Zone_EQtransit_InitFin <- rep(0, length=dim(Inv)[1]); # Corresponding codes of EQ_transit_if
EQ_transit_ResInv <- rep(0, length=dim(Inv)[1]); # Change in the system stability between the final outcomes of resident and invaded communities
Zone_EQtransit_ResInv <- rep(0, length=dim(Inv)[1]); # Corresponding codes of EQ_transit_invres

Inv <- cbind.data.frame(Inv, Rr_i, Ri_i, Cr_i, Dyn_state_i, Dyn_state_f, Zone_i, Zone_f, BioD, State);

#dynstati<-as.factor(c());
for(i in Alpha){
  for(j in Beta){
    gamma <- i*j;
    if(gamma<1){next}else{
      subinv <- Inv[which(Inv$Alpha==i & Inv$Beta==j),];
      subres <- Res[which(Res$Gamma==gamma),];
      Inv$Rr_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Rr_i;
      Inv$Cr_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Cr_i;
      Inv$Dyn_state_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Dyn_state_i;
      Inv$Zone_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Zone_i;
    }
  }
}

AC_Tot_data <- cbind.data.frame(
  Inv$Temperature, Inv$Carr, Fw, InvTP, Res_BMR, Inv$Alpha, Inv$Beta, Inv$Gamma, Inv$M_A, Inv$M_B, Inv$M_C,
  Inv$Rr_i, Inv$Ri_i, Inv$Cr_i, Inv$A, Inv$B, Inv$C, Inv$Afin, Inv$Bfin, Inv$Cfin, Inv$Tot, 
  Inv$ExtTimeA, Inv$ExtTimeB, Inv$ExtTimeC, Inv$Time, 
  Rr_eq, Ri_eq, Cr_eq, Dom_Eigenvalue, Tot_eq, Inv$Dyn_state_i, Inv$Dyn_state_f, Inv$Zone_i, Inv$Zone_f,
  EQ_transit_InitFin, Zone_EQtransit_InitFin, EQ_transit_ResInv, Zone_EQtransit_ResInv, Inv$BioD, Inv$State);

colnames(AC_Tot_data)<-c("Temp", "Car", "Fw", "Inv.position", "Res_BMR", "Alpha", "Beta", "Gamma", "M_Rr", "M_Ri", "M_Cr",
                         "Rr_i", "Ri_i", "Cr_i", "Rr_f", "Ri_f", "Cr_f", "Rr_fin", "Ri_fin", "Cr_fin", "Tot",
                         "Rr_extime", "Ri_extime", "Cr_extime", "Sim_duration", 
                         "Rr_eq", "Ri_eq", "Cr_eq", "Dom_Eigenvalue","Tot_eq",
                         "Dyn_state_i", "Dyn_state_f", "Zone_i", "Zone_f",
                         "EQ_transit_InitFin", "Zone_EQtransit_InitFin", "EQ_transit_ResInv", "Zone_EQtransit_ResInv", "BioD", "State");
#### Species composition -> Prule, Equilibrium analyses, Calculus of Eigenvalues ====
## 1. Cross-check of species diversity between transient analyses (after 5000 years) and the equilibrium analyses ====
# calculate prule only for 3sp remaining: to detect species exclusion scenarios from coexistence
data3sp <- AC_Tot_data[which(AC_Tot_data$Tot==3),];
PRdata <- as.data.frame(P.Rule.Appcomp(Temp=data3sp$Temp[1:dim(data3sp)[1]], Car=data3sp$Car[1:dim(data3sp)[1]], 
                                   MRr=data3sp$M_Rr[1:dim(data3sp)[1]], MRi=data3sp$M_Ri[1:dim(data3sp)[1]], MCr=data3sp$M_Cr[1:dim(data3sp)[1]]))
# correcting negative biomass
PRdata$C_star_Ri[which(PRdata$C_star_Ri<Tr)] <- 0;
PRdata$C_star_Rr[which(PRdata$C_star_Rr<Tr)] <- 0;
# Calculus of biomass ratio in C eq in C-Ri and C-Rr
C_eqir_ratio  <- PRdata$C_star_Ri/PRdata$C_star_Rr;
C_eqir_diff   <- PRdata$C_star_Ri-PRdata$C_star_Rr;
# species presence/absence
Rr_feq <- ifelse(PRdata$Rr_star > 1e-12, 1, 0);
Ri_feq <- ifelse(PRdata$Ri_star > 1e-12, 1, 0);
Cr_feq <- ifelse(PRdata$C_star_Rr > 1e-12 | PRdata$C_star_Ri > 1e-12, 1, 0);
Tot_eq <- rep(0, length=dim(PRdata)[1]);
Rrstar <- rep(0, length=dim(PRdata)[1]);
Ristar <- rep(0, length=dim(PRdata)[1]);
Crstar <- rep(0, length=dim(PRdata)[1]);
PRdata <- cbind.data.frame(PRdata[,1:9], C_eqir_ratio, C_eqir_diff, Rr_feq, Ri_feq, Cr_feq, Tot_eq, Rrstar, Ristar, Crstar);
# outcomes of P rules
PRdata$Rr_feq[which(PRdata$C_eqir_ratio>1)] <- 0;
PRdata$Ri_feq[which(PRdata$C_eqir_ratio<1)] <- 0;
PRdata$Rr_feq[which(PRdata$C_eqir_ratio==0 & PRdata$C_star_Rr==0)] <- 0;
PRdata$Ri_feq[which(PRdata$C_eqir_ratio==0 & PRdata$C_star_Ri==0)] <- 0;
PRdata$Tot_eq <- apply(PRdata[,c(12:15)], 1, sum);
# Biomass at equilibrium for 3sp system
PRdata$Rrstar[which(PRdata$Tot_eq==3)] <- PRdata$Rr_star[which(PRdata$Tot_eq==3)];
PRdata$Ristar[which(PRdata$Tot_eq==3)] <- PRdata$Ri_star[which(PRdata$Tot_eq==3)];
PRdata$Crstar[which(PRdata$Tot_eq==3)] <- PRdata$C_star_Rr[which(PRdata$Tot_eq==3)];

PRdata$Rrstar[which(PRdata$Tot_eq==2 & PRdata$Rr_feq == 1)] <- PRdata$Rr_star[which(PRdata$Tot_eq==2 & PRdata$Rr_feq == 1)];
PRdata$Crstar[which(PRdata$Tot_eq==2 & PRdata$Rr_feq == 1)] <- PRdata$C_star_Rr[which(PRdata$Tot_eq==2 & PRdata$Rr_feq == 1)];
PRdata$Ristar[which(PRdata$Tot_eq==2 & PRdata$Ri_feq == 1)] <- PRdata$Ri_star[which(PRdata$Tot_eq==2 & PRdata$Ri_feq == 1)];
PRdata$Crstar[which(PRdata$Tot_eq==2 & PRdata$Ri_feq == 1)] <- PRdata$C_star_Ri[which(PRdata$Tot_eq==2 & PRdata$Ri_feq == 1)];

AC_Tot_data$Tot_eq[which(AC_Tot_data$Tot==3)] <- PRdata$Tot_eq;
AC_Tot_data$Rr_eq[which(AC_Tot_data$Tot==3)] <- PRdata$Rrstar;
AC_Tot_data$Ri_eq[which(AC_Tot_data$Tot==3)] <- PRdata$Ristar;
AC_Tot_data$Cr_eq[which(AC_Tot_data$Tot==3)] <- PRdata$Crstar;
## 2. Calculus of Eigenvalues & Defining species composition ====
AC_Tot_data$Dyn_state_f <- as.character(AC_Tot_data$Dyn_state_f);
## Case 0: All species are present
data_3sp <- AC_Tot_data[which(AC_Tot_data$Tot_eq == 3),]
AC_Tot_data$Dyn_state_f[which(AC_Tot_data$Tot_eq == 3)] <- paste(3, "species", sep='.');
AC_Tot_data$Zone_f[which(AC_Tot_data$Dyn_state_f == paste(3, "species", sep='.') )] <- 8;
## Case 1: resident consumer remain on invading resources (Cr & Ri)
# A: From Tot == 3
data_CrRiA <- AC_Tot_data[which(AC_Tot_data$Tot == 3 & AC_Tot_data$Tot_eq == 2 & AC_Tot_data$Rr_eq < Tr),];
lambda_CrRiA <- mapply(lambda_2sp_BC, Temp=data_CrRiA$Temp[1:dim(data_CrRiA)[1]], Car=data_CrRiA$Car[1:dim(data_CrRiA)[1]],
                      By=data_CrRiA$Ri_eq[1:dim(data_CrRiA)[1]], Bz=data_CrRiA$Cr_eq[1:dim(data_CrRiA)[1]],
                      My=data_CrRiA$M_Ri[1:dim(data_CrRiA)[1]], Mz=data_CrRiA$M_Cr[1:dim(data_CrRiA)[1]],
                      Beta=data_CrRiA$Beta[1:dim(data_CrRiA)[1]]);
data_CrRiA <- cbind.data.frame(data_CrRiA, lambda_CrRiA);
AC_Tot_data$Dom_Eigenvalue[which(AC_Tot_data$Tot == 3 & AC_Tot_data$Tot_eq == 2 & AC_Tot_data$Rr_eq < Tr)] <- data_CrRiA$lambda_CrRiA;
# B: From Tot == 2
data_CrRiB <- AC_Tot_data[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Rr_fin == 0),];
data_CrRiBeq <- as.data.frame(EQUI_2spBC(Temp=data_CrRiB$Temp[1:dim(data_CrRiB)[1]], Car=data_CrRiB$Car[1:dim(data_CrRiB)[1]],
                                         Alpha=data_CrRiB$Alpha[1:dim(data_CrRiB)[1]],Beta=data_CrRiB$Beta[1:dim(data_CrRiB)[1]]));
data_CrRiB$Ri_eq <- data_CrRiBeq$By;
data_CrRiB$Cr_eq <- data_CrRiBeq$Bz;

lambda_CrRiB <- mapply(lambda_2sp_BC, Temp=data_CrRiB$Temp[1:dim(data_CrRiB)[1]], Car=data_CrRiB$Car[1:dim(data_CrRiB)[1]],
                      By=data_CrRiB$Ri_eq[1:dim(data_CrRiB)[1]], Bz=data_CrRiB$Cr_eq[1:dim(data_CrRiB)[1]],
                      My=data_CrRiB$M_Ri[1:dim(data_CrRiB)[1]], Mz=data_CrRiB$M_Cr[1:dim(data_CrRiB)[1]],
                      Beta=data_CrRiB$Beta[1:dim(data_CrRiB)[1]]);
data_CrRiB <- cbind.data.frame(data_CrRiB, lambda_CrRiB);
AC_Tot_data$Tot_eq[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Rr_fin == 0)] <- 2;
AC_Tot_data$Ri_eq[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Rr_fin == 0)] <- data_CrRiB$Ri_eq;
AC_Tot_data$Cr_eq[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Rr_fin == 0)] <- data_CrRiB$Cr_eq;
AC_Tot_data$Dom_Eigenvalue[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Rr_fin == 0)] <- data_CrRiB$lambda_CrRiB;

AC_Tot_data$Dyn_state_f[which(AC_Tot_data$Tot_eq ==2 & AC_Tot_data$Rr_eq < Tr & AC_Tot_data$Ri_eq > Tr & AC_Tot_data$Dom_Eigenvalue > 0)] <- paste("Cr", "Ri", "cycles", sep='.');
AC_Tot_data$Dyn_state_f[which(AC_Tot_data$Tot_eq ==2 & AC_Tot_data$Rr_eq < Tr & AC_Tot_data$Ri_eq > Tr & (AC_Tot_data$Dom_Eigenvalue == 0 | AC_Tot_data$Dom_Eigenvalue < 0) )] <- paste("Cr", "Ri", "eq", sep='.');
AC_Tot_data$Zone_f[which(AC_Tot_data$Dyn_state_f == paste("Cr", "Ri", "eq", sep='.') )] <- 4;
AC_Tot_data$Zone_f[which(AC_Tot_data$Dyn_state_f == paste("Cr", "Ri", "cycles", sep='.') )] <- 5;
## Case 2: resident species remain (Cr & Rr)
# A: From Tot == 3
data_CrRrA <- AC_Tot_data[which(AC_Tot_data$Tot == 3 & AC_Tot_data$Tot_eq == 2 & AC_Tot_data$Ri_eq < Tr),];
lambda_CrRrA <- mapply(lambda_2sp_AC, Temp=data_CrRrA$Temp[1:dim(data_CrRrA)[1]], Car=data_CrRrA$Car[1:dim(data_CrRrA)[1]],
                      Bx=data_CrRrA$Rr_eq[1:dim(data_CrRrA)[1]], Bz=data_CrRrA$Cr_eq[1:dim(data_CrRrA)[1]],
                      Mx=data_CrRrA$M_Rr[1:dim(data_CrRrA)[1]], Mz=data_CrRrA$M_Cr[1:dim(data_CrRrA)[1]],
                      Gamma=data_CrRrA$Gamma[1:dim(data_CrRrA)[1]]);
data_CrRrA <- cbind.data.frame(data_CrRrA, lambda_CrRrA);
AC_Tot_data$Dom_Eigenvalue[which(AC_Tot_data$Tot == 3 & AC_Tot_data$Tot_eq == 2 & AC_Tot_data$Ri_eq < Tr)] <- data_CrRrA$lambda_CrRrA;
# B: From Tot == 2
data_CrRrB <- AC_Tot_data[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Ri_fin == 0),];
data_CrRrBeq <- as.data.frame(EQUI_2spAC(Temp=data_CrRrB$Temp[1:dim(data_CrRrB)[1]], Car=data_CrRrB$Car[1:dim(data_CrRrB)[1]],
                                         Gamma=data_CrRrB$Gamma[1:dim(data_CrRrB)[1]]));
data_CrRrB$Rr_eq <- data_CrRrBeq$Bx;
data_CrRrB$Cr_eq <- data_CrRrBeq$Bz;

lambda_CrRrB <- mapply(lambda_2sp_AC, Temp=data_CrRrB$Temp[1:dim(data_CrRrB)[1]], Car=data_CrRrB$Car[1:dim(data_CrRrB)[1]],
                       Bx=data_CrRrB$Rr_eq[1:dim(data_CrRrB)[1]], Bz=data_CrRrB$Cr_eq[1:dim(data_CrRrB)[1]],
                       Mx=data_CrRrB$M_Rr[1:dim(data_CrRrB)[1]], Mz=data_CrRrB$M_Cr[1:dim(data_CrRrB)[1]],
                       Gamma=data_CrRrB$Gamma[1:dim(data_CrRrB)[1]]);
data_CrRrB <- cbind.data.frame(data_CrRrB, lambda_CrRrB);
AC_Tot_data$Tot_eq[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Ri_fin == 0)] <- 2;
AC_Tot_data$Rr_eq[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Ri_fin == 0)] <- data_CrRrB$Rr_eq;
AC_Tot_data$Cr_eq[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Ri_fin == 0)] <- data_CrRrB$Cr_eq;
AC_Tot_data$Dom_Eigenvalue[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Ri_fin == 0)] <- data_CrRrB$lambda_CrRrB;

AC_Tot_data$Dyn_state_f[which(AC_Tot_data$Tot_eq ==2 & AC_Tot_data$Rr_eq > Tr & AC_Tot_data$Ri_eq < Tr & AC_Tot_data$Dom_Eigenvalue > 0)] <- paste("Cr", "Rr", "cycles", sep='.');
AC_Tot_data$Dyn_state_f[which(AC_Tot_data$Tot_eq ==2 & AC_Tot_data$Rr_eq > Tr & AC_Tot_data$Ri_eq < Tr & (AC_Tot_data$Dom_Eigenvalue == 0 | AC_Tot_data$Dom_Eigenvalue < 0) )] <- paste("Cr", "Rr", "eq", sep='.');

AC_Tot_data$Zone_f[which(AC_Tot_data$Dyn_state_f == paste("Cr", "Rr", "eq", sep='.') )] <- 2;
AC_Tot_data$Zone_f[which(AC_Tot_data$Dyn_state_f == paste("Cr", "Rr", "cycles", sep='.') )] <- 3;
# Case 3: only resources remain (either Rr alone, Ri alone, or both resource)
# A: From Tot == 2
dataRA <- AC_Tot_data[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Cr_fin ==0),];
#dataRB <- AC_Tot_data[which(AC_Tot_data$Tot == 1 & AC_Tot_data$Ri_fin ==0 & AC_Tot_data$Cr_fin ==0),];
#dataRC <- AC_Tot_data[which(AC_Tot_data$Tot == 1 & AC_Tot_data$Rr_fin ==0 & AC_Tot_data$Cr_fin ==0),];
dataRA$Rr_eq <- K(Temp= dataRA$Temp, Mx= dataRA$M_Rr,Car=dataRA$Car)
dataRA$Ri_eq <- K(Temp= dataRA$Temp, Mx= dataRA$M_Ri,Car=dataRA$Car)
AC_Tot_data$Rr_eq[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Cr_eq < Tr)] <- dataRA$Rr_eq;
AC_Tot_data$Ri_eq[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Cr_eq < Tr)] <- dataRA$Ri_eq;
AC_Tot_data$Tot_eq[which(AC_Tot_data$Tot == 2 & AC_Tot_data$Cr_eq < Tr)] <- 2;

#AC_Tot_data$Dyn_state_f[which(AC_Tot_data$Tot_eq == 1)] <- paste("R", "eq", sep=".");
AC_Tot_data$Dyn_state_f[which(AC_Tot_data$Tot_eq == 2 & AC_Tot_data$Cr_eq < Tr)] <- paste("R", "eq", sep=".");
AC_Tot_data$Zone_f[which(AC_Tot_data$Dyn_state_f == paste("R", "eq", sep="."))] <- 1;

# Case 0: no species remaining
AC_Tot_data$Dyn_state_f[which(AC_Tot_data$Tot_eq == 0)] <- paste(0, "species", sep=".");
AC_Tot_data$Zone_f[which(AC_Tot_data$Dyn_state_f == paste(0, "species", sep="."))] <- 0;

#### Biodiversity change ====
for(i in Alpha){
  for(j in Beta){
    gamma <- i*j;
    if(gamma<0.9){next}else{
      subinv <- AC_Tot_data[which(AC_Tot_data$Alpha == i & AC_Tot_data$Beta == j),];
      subres <- Res[which(Res$Gamma==gamma),];
      AC_Tot_data$BioD[which(AC_Tot_data$Alpha == i & AC_Tot_data$Beta == j)] <- (subinv$Tot_eq - subres$Tot_eq);
    }
  }
}
#### Invasion states ====
AC_Tot_data$State <- as.character(AC_Tot_data$State);
AC_Tot_data$State[which(AC_Tot_data$Tot_eq == 3)] <- "Integration";
AC_Tot_data$State[which(AC_Tot_data$BioD == 0 & AC_Tot_data$Tot_eq == 2 & AC_Tot_data$Rr_eq < Tr)] <- "Substitution";
AC_Tot_data$State[which(AC_Tot_data$BioD > 0 & AC_Tot_data$Tot_eq == 2 & AC_Tot_data$Rr_eq < Tr)] <- "Substitution";
AC_Tot_data$State[which(AC_Tot_data$BioD == 0 & AC_Tot_data$Tot_eq == 0)]  <- "Resistance";
AC_Tot_data$State[which(AC_Tot_data$BioD == 0 & AC_Tot_data$Tot_eq == 1)] <- "Resistance";
AC_Tot_data$State[which(AC_Tot_data$BioD == 0 & AC_Tot_data$Tot_eq == 2 & AC_Tot_data$Ri_eq < Tr)] <- "Resistance";

AC_Tot_data$State[which(AC_Tot_data$BioD < 0)] <- "Vulnerability";
AC_Tot_data$State[which(AC_Tot_data$BioD > 0 & AC_Tot_data$Ri_eq > Tr & AC_Tot_data$Cr_eq < Tr & AC_Tot_data$Tot_eq != 3)] <- "Occupancy";
AC_Tot_data$State[which(AC_Tot_data$BioD > 0 & AC_Tot_data$Ri_eq < Tr & AC_Tot_data$Tot_eq != 3)] <- "Rescue";
#### Changes in system stability in the community (initially vs. after 5000 years) ====
# 0. Initial-Final configurations
conf_name <- c(); initeq <- c("Eq", "Cycles"); fineq <- c("Null", "Eq", "Cycles")
for(i in initeq){conf_name <- append(conf_name, paste(i, fineq, sep='.'))};
conf_code <- c(0:5);
Conf_tab <- cbind.data.frame(conf_name, conf_code);

summary(as.factor(AC_Tot_data$Zone_i))
summary(as.factor(AC_Tot_data$Zone_f))
# A. Prepare the different combinations of stability changes
init_eq_code <- c(1, 2, 3);
init_eq <- c(rep("Eq", 2), "Cycles");
final_eq_code <- c(seq(0,5,by=1), 8);
final_eq <- c("Null", "Eq", rep( c("Eq", "Cycles"), 2), "Eq");

# Data of observe combinations
I_code <- c(); for(i in init_eq_code){I_code <- append(I_code, rep(i, length(final_eq)))};
F_code <- rep(final_eq_code, length(init_eq_code));

IF_eq <- c(); for(i in init_eq){IF_eq <- append(IF_eq, paste(i, final_eq, sep='.'))}
IF_code <- rep(0, length(IF_eq));
IF <- cbind.data.frame(IF_eq, IF_code, I_code, F_code);

IF$IF_eq <- as.character(IF$IF_eq);
for(i in 1:length(IF_eq)){ IF$IF_code[i] <- Conf_tab$conf_code[which(Conf_tab$conf_name==IF_eq[i])]};
# B. Attribute the AC_Tot_data's outcomes to the corresponding configurations EQ code and name
AC_Tot_data$Zone_i <- as.numeric(AC_Tot_data$Zone_i);
AC_Tot_data$Zone_f <- as.numeric(AC_Tot_data$Zone_f);
AC_Tot_data$EQ_transit_InitFin <- as.character(AC_Tot_data$EQ_transit_InitFin)
AC_Tot_data$Zone_EQtransit_InitFin <- as.numeric(AC_Tot_data$Zone_EQtransit_InitFin);
for(i in 1:dim(IF)[1]){
  AC_Tot_data$EQ_transit_InitFin[which(AC_Tot_data$Zone_i == IF$I_code[i] & AC_Tot_data$Zone_f==IF$F_code[i])] <- IF$IF_eq[i];
  AC_Tot_data$Zone_EQtransit_InitFin[which(AC_Tot_data$Zone_i == IF$I_code[i] & AC_Tot_data$Zone_f==IF$F_code[i])] <- IF$IF_code[i]
}
#### Changes in system stability between resident and invaded communities (both after 5000 years) ====
Res_zonef <- rep(0, dim(AC_Tot_data)[1]);
AC_Tot_data <- cbind.data.frame(AC_Tot_data, Res_zonef);
# 0. Resident-Invaded configurations
conf_name <- c(); initeq <- c("Null", "Eq", "Cycles");
for(i in initeq){conf_name <- append(conf_name, paste(i, initeq, sep='.'))};
conf_code <- c(0:8);
Conf_tab <- cbind.data.frame(conf_name, conf_code);
# A. Merge Res$Zone_f to invaded data frame
for(i in Alpha){
  for(j in Beta){
    gamma <- i*j;
    if(gamma<0.9){next}else{
      subinv <- subset(AC_Tot_data,AC_Tot_data$Alpha==i & AC_Tot_data$Beta==j);
      subres <- subset(Res,Res$Gamma==gamma);
      AC_Tot_data$Res_zonef[which(AC_Tot_data$Alpha == i & AC_Tot_data$Beta == j)] <- subres$Zone_f;
    }
  }
}
# B. Prepare the different combinations of stability changes
res_eq_code <- c(0, 1, 2, 3);
res_eq <- c("Null", rep("Eq", 2), "Cycles");
# EQ & Code for invaded community
inv_eq_code <- c(seq(0,5,by=1), 8);
inv_eq <- c("Null", "Eq", rep( c("Eq", "Cycles"), 2), "Eq");
# Data of observe combinations
Res_code <- c(); for(i in res_eq_code){Res_code <- append(Res_code, rep(i, length(inv_eq)))};
Inv_code <- rep(inv_eq_code, length(res_eq_code));

Inv_Res_eq <- c(); for(i in res_eq){Inv_Res_eq <- append(Inv_Res_eq, paste(i, inv_eq, sep='.'))}
Inv_Res_code <- rep(0, length(Inv_Res_eq));
Inv_Res <- cbind.data.frame(Inv_Res_eq, Inv_Res_code, Res_code, Inv_code);
Inv_Res$Inv_Res_eq <- as.character(Inv_Res$Inv_Res_eq);

for(i in 1:length(Inv_Res_eq)){ Inv_Res$Inv_Res_code[i] <- Conf_tab$conf_code[which(Conf_tab$conf_name==Inv_Res_eq[i])]};
# C. Attribute the AC_Tot_data's outcomes to the corresponding configurations EQ code and name
AC_Tot_data$EQ_transit_ResInv <- as.character(AC_Tot_data$EQ_transit_ResInv)
AC_Tot_data$Zone_f <- as.numeric(AC_Tot_data$Zone_f);
AC_Tot_data$Res_zonef <- as.numeric(AC_Tot_data$Res_zonef);
AC_Tot_data$Zone_EQtransit_ResInv <- as.numeric(AC_Tot_data$Zone_EQtransit_ResInv);
for(i in 1:dim(Inv_Res)[1]){
  AC_Tot_data$EQ_transit_ResInv[which(AC_Tot_data$Res_zonef == Inv_Res$Res_code[i] & AC_Tot_data$Zone_f==Inv_Res$Inv_code[i])] <- Inv_Res$Inv_Res_eq[i];
  AC_Tot_data$Zone_EQtransit_ResInv[which(AC_Tot_data$Res_zonef == Inv_Res$Res_code[i] & AC_Tot_data$Zone_f==Inv_Res$Inv_code[i])] <- Inv_Res$Inv_Res_code[i]
}
#### Savig the dataset ====
write.table(AC_Tot_data, "./data/AC_Tot_data.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
############# EC; Exploitative competition ====
Inv <- read.table("./data/EC.data.txt", h=T,  sep = "\t");     # Data from EC module simulation
Res <- read.table("./data/CR_Tot_Gamma.txt", h=T,  sep = "\t");# Complete Data from 2sp C-R (A, C) including simulation outcomes and equilibrium analyses

Alpha <- c(1, 2, 5, 10);                          # BMR B/A
Beta <- c(0.1, 0.2, 0.5, 0.75, 1, 1.5, 2, 5, 10);            # BMR C/B
Gamma <- c(1, 1.5, 2, 2.5, 3, 3.75, 4, 5, 7.5, 10, 15, 20, 25, 50, 100); 
Tr <- 1e-12;      # Extinction threshold

#### Initiating the dataset ====
InvTP <- rep(2, dim(Inv)[1]); # Trophic position of the invading species (here 2 = consumer)
Fw <- rep("EC", dim(Inv)[1]); # Food web module in which invasion occurs (here EC = Exploitative competition)
Res_BMR <- rep('Gamma', length=dim(Inv)[1]);# species bmr of resident system

Dyn_state_i <- rep('NA', length=dim(Inv)[1]);# Dynamic state initial of resident species
Dyn_state_f <- rep('NA', length=dim(Inv)[1]);# Dynamic state final, after invasion
Zone_i <- rep(0, length=dim(Inv)[1]);     # initial, identical to 2spAC
Zone_f <- rep(0, length=dim(Inv)[1]);     # final
BioD <- rep(0, length=dim(Inv)[1]);       # Biodiversity change after invasion
State <- rep(0, length=dim(Inv)[1]);      # Invasion states
Rr_i <- rep(0, length=dim(Inv)[1]);       # Initial biomass density of resident resource (starting at C-R eq)
Cr_i <- rep(0, length=dim(Inv)[1]);       # Initial biomass density of resident consumer (starting at C-R eq)
Ci_i <- rep(1E-6, length=dim(Inv)[1]);    # Initial biomass density of invading consumer

Rr_eq <- rep(0, length=dim(Inv)[1]);      # Biomass density of resident resource at equilibrium
Ci_eq <- rep(0, length=dim(Inv)[1]);      # Biomass density of invading consumer at equilibrium
Cr_eq <- rep(0, length=dim(Inv)[1]);      # Biomass density of resident consumer at equilibrium
Tot_eq <- rep(0, length=dim(Inv)[1]);     # Total diversity at equilibrium
Dom_Eigenvalue <- rep(0, length=dim(Inv)[1]); # Dominant eigenvalues for the different community structure

EQ_transit_InitFin <- rep(0, length=dim(Inv)[1]); # Change in the system stability between the initial and final outcomes
Zone_EQtransit_InitFin <- rep(0, length=dim(Inv)[1]); # Corresponding codes of EQ_transit_if
EQ_transit_ResInv <- rep(0, length=dim(Inv)[1]); # Change in the system stability between the final outcomes of resident and invaded communities
Zone_EQtransit_ResInv <- rep(0, length=dim(Inv)[1]); # Corresponding codes of EQ_transit_invres

Inv <- cbind.data.frame(Inv, Rr_i, Ci_i, Cr_i, Dyn_state_i, Dyn_state_f, Zone_i, Zone_f, BioD, State);
for(i in Alpha){
  for(j in Beta){
    gamma <- i*j;
    if(gamma<0.9){next}else{
      subinv <- Inv[which(Inv$Alpha==i & Inv$Beta==j),];
      subres <- Res[which(Res$Gamma==gamma),];
      Inv$Rr_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Rr_i;
      Inv$Cr_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Cr_i;
      Inv$Zone_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Zone_i;
      Inv$Dyn_state_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Dyn_state_i;
    }
  }
}

EC_Tot_data <- cbind.data.frame(
  Inv$Temperature, Inv$Carr, Fw, InvTP, Res_BMR, Inv$Alpha, Inv$Beta, Inv$Gamma, Inv$M_A, Inv$M_B, Inv$M_C,
  Inv$Rr_i, Inv$Ci_i, Inv$Cr_i, Inv$A, Inv$B, Inv$C, Inv$Afin, Inv$Bfin, Inv$Cfin, Inv$Tot, 
  Inv$ExtTimeA, Inv$ExtTimeB, Inv$ExtTimeC, Inv$Time,
  Rr_eq, Ci_eq, Cr_eq, Dom_Eigenvalue, Tot_eq, Inv$Dyn_state_i, Inv$Dyn_state_f, Inv$Zone_i, Inv$Zone_f,
  EQ_transit_InitFin, Zone_EQtransit_InitFin, EQ_transit_ResInv, Zone_EQtransit_ResInv, Inv$BioD, Inv$State);

colnames(EC_Tot_data)<-c("Temp", "Car", "Fw", "Inv.position", "Res_BMR", "Alpha", "Beta", "Gamma", "M_Rr", "M_Ci", "M_Cr",
                         "Rr_i", "Ci_i", "Cr_i", "Rr_f", "Ci_f", "Cr_f", "Rr_fin", "Ci_fin", "Cr_fin", "Tot",
                         "Rr_extime", "Ci_extime", "Cr_extime", "Sim_duration",
                         "Rr_eq", "Ci_eq", "Cr_eq", "Dom_Eigenvalue", "Tot_eq",
                         "Dyn_state_i", "Dyn_state_f", "Zone_i", "Zone_f",
                         "EQ_transit_InitFin", "Zone_EQtransit_InitFin", "EQ_transit_ResInv", "Zone_EQtransit_ResInv", "BioD", "State");

#### Species composition -> Rrule, Equilibrium analyses, Calculus of Eigenvalues ====
## 1. Cross-check of species diversity between transient analyses (after 5000 years) and the equilibrium analyses ====
data3sp <- EC_Tot_data[which(EC_Tot_data$Tot==3),];
PRdata <- as.data.frame(R.Rule.Expcomp(Temp=data3sp$Temp[1:dim(data3sp)[1]], Car=data3sp$Car[1:dim(data3sp)[1]], 
                                   MRr=data3sp$M_Rr[1:dim(data3sp)[1]], MCi=data3sp$M_Ci[1:dim(data3sp)[1]], MCr=data3sp$M_Cr[1:dim(data3sp)[1]]))
# Calculus of biomass ratio in R eq in Ci-R and Cr-R
R_eqir_ratio  <- PRdata$R_star_Ci/PRdata$R_star_Cr;
R_eqir_diff   <- PRdata$R_star_Ci-PRdata$R_star_Cr;
# species presence/absence
Rr_feq <- ifelse(PRdata$R_star_Cr > 1e-12 | PRdata$R_star_Ci > 1e-12, 1, 0);
Ci_feq <- ifelse(PRdata$Ci_star > 1e-12, 1, 0);
Cr_feq <- ifelse(PRdata$Cr_star > 1e-12, 1, 0);
Tot_eq <- rep(0, length=dim(PRdata)[1]);
Rrstar <- rep(0, length=dim(PRdata)[1]);
Cistar <- rep(0, length=dim(PRdata)[1]);
Crstar <- rep(0, length=dim(PRdata)[1]);
PRdata <- cbind.data.frame(PRdata[,1:9], R_eqir_ratio, R_eqir_diff, Rr_feq, Ci_feq, Cr_feq, Tot_eq, Rrstar, Cistar, Crstar);
# outcomes of P rules
PRdata$Ci_feq[which(PRdata$R_eqir_ratio>1)] <- 0;
PRdata$Cr_feq[which(PRdata$R_eqir_ratio<1)] <- 0;
PRdata$Tot_eq <- apply(PRdata[,c(12:14)], 1, sum);
# Biomass at equilibrium for 3sp system
PRdata$Rrstar[which(PRdata$Tot_eq==3)] <- PRdata$R_star_Cr[which(PRdata$Tot_eq==3)];
PRdata$Cistar[which(PRdata$Tot_eq==3)] <- PRdata$Ci_star[which(PRdata$Tot_eq==3)];
PRdata$Crstar[which(PRdata$Tot_eq==3)] <- PRdata$Cr_star[which(PRdata$Tot_eq==3)];
# correcting 2-sp biomass at equilibrium (Ci-Rr & Cr-Rr)
PRdata$Rrstar[which(PRdata$Tot_eq==2 & PRdata$Cr_feq == 0)] <- PRdata$R_star_Ci[which(PRdata$Tot_eq==2 & PRdata$Cr_feq == 0)];
PRdata$Cistar[which(PRdata$Tot_eq==2 & PRdata$Cr_feq == 0)] <- PRdata$Ci_star[which(PRdata$Tot_eq==2 & PRdata$Cr_feq == 0)];

PRdata$Rrstar[which(PRdata$Tot_eq==2 & PRdata$Ci_feq == 0)] <- PRdata$R_star_Cr[which(PRdata$Tot_eq==2 & PRdata$Ci_feq == 0)];
PRdata$Crstar[which(PRdata$Tot_eq==2 & PRdata$Ci_feq == 0)] <- PRdata$Crstar[which(PRdata$Tot_eq==2 & PRdata$Ci_feq == 0)];

EC_Tot_data$Tot_eq[which(EC_Tot_data$Tot==3)] <- PRdata$Tot_eq;
EC_Tot_data$Rr_eq[which(EC_Tot_data$Tot==3)] <- PRdata$Rrstar;
EC_Tot_data$Ci_eq[which(EC_Tot_data$Tot==3)] <- PRdata$Cistar;
EC_Tot_data$Cr_eq[which(EC_Tot_data$Tot==3)] <- PRdata$Crstar;
## 2. Calculus of Eigenvalues & Defining species composition ====
EC_Tot_data$Dyn_state_f <- as.character(EC_Tot_data$Dyn_state_f);
## Case 0: All species are present
data_3sp <- EC_Tot_data[which(EC_Tot_data$Tot_eq == 3),]
EC_Tot_data$Dyn_state_f[which(EC_Tot_data$Tot_eq == 3)] <- paste(3, "species", sep='.');
EC_Tot_data$Zone_f[which(EC_Tot_data$Dyn_state_f == paste(3, "species", sep='.') )] <- 8;
# Case 1: Invading consumer replaced the resident species (Ci & Rr)
# A: From Tot == 3
data_CiRrA <- EC_Tot_data[which(EC_Tot_data$Tot == 3 & EC_Tot_data$Tot_eq == 2 & EC_Tot_data$Cr_eq < Tr),];
lambda_CiRrA <- mapply(lambda_2sp_AB, Temp=data_CiRrA$Temp[1:dim(data_CiRrA)[1]], Car=data_CiRrA$Car[1:dim(data_CiRrA)[1]],
                       Bx=data_CiRrA$Rr_eq[1:dim(data_CiRrA)[1]], By=data_CiRrA$Ci_eq[1:dim(data_CiRrA)[1]],
                       Mx=data_CiRrA$M_Rr[1:dim(data_CiRrA)[1]], My=data_CiRrA$M_Ci[1:dim(data_CiRrA)[1]],
                       Alpha=data_CiRrA$Alpha[1:dim(data_CiRrA)[1]]);
data_CiRrA <- cbind.data.frame(data_CiRrA, lambda_CiRrA);
EC_Tot_data$Dom_Eigenvalue[which(EC_Tot_data$Tot == 3 & EC_Tot_data$Tot_eq == 2 & EC_Tot_data$Cr_eq < Tr)] <- data_CiRrA$lambda_CiRrA;
# B: From Tot == 2
data_CiRrB <- EC_Tot_data[which(EC_Tot_data$Tot == 2 & EC_Tot_data$Ci_fin == 1),];
data_CiRrBeq <- as.data.frame(EQUI_2spAB(Temp=data_CiRrB$Temp[1:dim(data_CiRrB)[1]], Car=data_CiRrB$Car[1:dim(data_CiRrB)[1]],
                                         Alpha=data_CiRrB$Alpha[1:dim(data_CiRrB)[1]]));
data_CiRrB$Rr_eq <- data_CiRrBeq$Bx;
data_CiRrB$Ci_eq <- data_CiRrBeq$By;

lambda_CiRrB <- mapply(lambda_2sp_AB, Temp=data_CiRrB$Temp[1:dim(data_CiRrB)[1]], Car=data_CiRrB$Car[1:dim(data_CiRrB)[1]],
                      Bx=data_CiRrB$Rr_eq[1:dim(data_CiRrB)[1]], By=data_CiRrB$Ci_eq[1:dim(data_CiRrB)[1]],
                      Mx=data_CiRrB$M_Rr[1:dim(data_CiRrB)[1]], My=data_CiRrB$M_Ci[1:dim(data_CiRrB)[1]],
                      Alpha=data_CiRrB$Alpha[1:dim(data_CiRrB)[1]]);
data_CiRrB <- cbind.data.frame(data_CiRrB, lambda_CiRrB);
EC_Tot_data$Tot_eq[which(EC_Tot_data$Tot == 2 & EC_Tot_data$Ci_fin == 1)] <- 2;
EC_Tot_data$Rr_eq[which(EC_Tot_data$Tot == 2 & EC_Tot_data$Ci_fin == 1)] <- data_CiRrB$Rr_eq;
EC_Tot_data$Ci_eq[which(EC_Tot_data$Tot == 2 & EC_Tot_data$Ci_fin == 1)] <- data_CiRrB$Ci_eq;
EC_Tot_data$Dom_Eigenvalue[which(EC_Tot_data$Tot == 2 & EC_Tot_data$Ci_fin == 1)] <- data_CiRrB$lambda_CiRrB;
# C: Defining system stability in Ci-Rr
EC_Tot_data$Dyn_state_f[which(EC_Tot_data$Tot_eq ==2 & EC_Tot_data$Ci_eq > Tr & EC_Tot_data$Dom_Eigenvalue > 0)] <- paste("Ci", "Rr", "cycles", sep='.');
EC_Tot_data$Dyn_state_f[which(EC_Tot_data$Tot_eq ==2 & EC_Tot_data$Ci_eq > Tr & (EC_Tot_data$Dom_Eigenvalue == 0 | EC_Tot_data$Dom_Eigenvalue < 0) )] <- paste("Ci", "Rr", "eq", sep='.');
EC_Tot_data$Zone_f[which(EC_Tot_data$Dyn_state_f == paste("Ci", "Rr", "cycles", sep='.') )] <- 7;
EC_Tot_data$Zone_f[which(EC_Tot_data$Dyn_state_f == paste("Ci", "Rr", "eq", sep='.') )] <- 6;
# Case 2: Only resident species remain (Cr & Rr)
# A: From Tot == 3
data_CrRrA <- EC_Tot_data[which(EC_Tot_data$Tot == 3 & EC_Tot_data$Tot_eq == 2 & EC_Tot_data$Ci_eq < Tr),];
lambda_CrRrA <- mapply(lambda_2sp_AC, Temp=data_CrRrA$Temp[1:dim(data_CrRrA)[1]], Car=data_CrRrA$Car[1:dim(data_CrRrA)[1]],
                       Bx=data_CrRrA$Rr_eq[1:dim(data_CrRrA)[1]], Bz=data_CrRrA$Cr_eq[1:dim(data_CrRrA)[1]],
                       Mx=data_CrRrA$M_Rr[1:dim(data_CrRrA)[1]], Mz=data_CrRrA$M_Cr[1:dim(data_CrRrA)[1]],
                       Gamma=data_CrRrA$Gamma[1:dim(data_CrRrA)[1]]);
data_CrRrA <- cbind.data.frame(data_CrRrA, lambda_CrRrA);
EC_Tot_data$Dom_Eigenvalue[which(EC_Tot_data$Tot == 3 & EC_Tot_data$Tot_eq == 2 & EC_Tot_data$Ci_eq < Tr)] <- data_CrRrA$lambda_CrRrA;
# B: From Tot == 2
data_CrRrB <- EC_Tot_data[which(EC_Tot_data$Tot == 2 & EC_Tot_data$Cr_fin == 1),];
data_CrRrBeq <- as.data.frame(EQUI_2spAC(Temp=data_CrRrB$Temp[1:dim(data_CrRrB)[1]], Car=data_CrRrB$Car[1:dim(data_CrRrB)[1]],
                                         Gamma=data_CrRrB$Gamma[1:dim(data_CrRrB)[1]]));
data_CrRrB$Rr_eq <- data_CrRrBeq$Bx;
data_CrRrB$Cr_eq <- data_CrRrBeq$Bz;

lambda_CrRrB <- mapply(lambda_2sp_AC, Temp=data_CrRrB$Temp[1:dim(data_CrRrB)[1]], Car=data_CrRrB$Car[1:dim(data_CrRrB)[1]],
                       Bx=data_CrRrB$Rr_eq[1:dim(data_CrRrB)[1]], Bz=data_CrRrB$Cr_eq[1:dim(data_CrRrB)[1]],
                       Mx=data_CrRrB$M_Rr[1:dim(data_CrRrB)[1]], Mz=data_CrRrB$M_Cr[1:dim(data_CrRrB)[1]],
                       Gamma=data_CrRrB$Gamma[1:dim(data_CrRrB)[1]]);
data_CrRrB <- cbind.data.frame(data_CrRrB, lambda_CrRrB);
EC_Tot_data$Tot_eq[which(EC_Tot_data$Tot == 2 & EC_Tot_data$Cr_fin == 1)] <- 2;
EC_Tot_data$Rr_eq[which(EC_Tot_data$Tot == 2 & EC_Tot_data$Cr_fin == 1)] <- data_CrRrB$Rr_eq;
EC_Tot_data$Cr_eq[which(EC_Tot_data$Tot == 2 & EC_Tot_data$Cr_fin == 1)] <- data_CrRrB$Cr_eq;
EC_Tot_data$Dom_Eigenvalue[which(EC_Tot_data$Tot == 2 & EC_Tot_data$Cr_fin == 1)] <- data_CrRrB$lambda_CrRrB;
# C: Defining system stability in Cr-Rr
EC_Tot_data$Dyn_state_f[which(EC_Tot_data$Tot_eq ==2 & EC_Tot_data$Cr_eq > Tr & EC_Tot_data$Dom_Eigenvalue > 0)] <- paste("Cr", "Rr", "cycles", sep='.');
EC_Tot_data$Dyn_state_f[which(EC_Tot_data$Tot_eq ==2 & EC_Tot_data$Cr_eq > Tr & (EC_Tot_data$Dom_Eigenvalue == 0 | EC_Tot_data$Dom_Eigenvalue < 0) )] <- paste("Cr", "Rr", "eq", sep='.');
EC_Tot_data$Zone_f[which(EC_Tot_data$Dyn_state_f == paste("Cr", "Rr", "cycles", sep='.') )] <- 3;
EC_Tot_data$Zone_f[which(EC_Tot_data$Dyn_state_f == paste("Cr", "Rr", "eq", sep='.') )] <- 2;
# Case 3: only resources remain
# A: From Tot == 3
dataRa <- EC_Tot_data[which(EC_Tot_data$Tot == 3 & EC_Tot_data$Tot_eq == 2 & EC_Tot_data$Ci_eq < Tr & EC_Tot_data$Cr_eq < Tr),];
dataRa$Rr_eq <- K(Temp= dataRa$Temp, Mx= dataRa$M_Rr,Car=dataRa$Car)
dataRa$Tot_eq <- 1;
EC_Tot_data$Rr_eq[which(EC_Tot_data$Tot == 3 & EC_Tot_data$Tot_eq == 2 & EC_Tot_data$Ci_eq < Tr & EC_Tot_data$Cr_eq < Tr)] <- dataRa$Rr_eq;
EC_Tot_data$Tot_eq[which(EC_Tot_data$Tot == 3 & EC_Tot_data$Tot_eq == 2 & EC_Tot_data$Ci_eq < Tr & EC_Tot_data$Cr_eq < Tr)] <- dataRa$Tot_eq;

# B: From Tot == 1
dataRb <- EC_Tot_data[which(EC_Tot_data$Tot == 1),];
dataRb$Rr_eq <- K(Temp= dataRb$Temp, Mx= dataRb$M_Rr,Car=dataRb$Car)
EC_Tot_data$Rr_eq[which(EC_Tot_data$Tot == 1)] <- dataRb$Rr_eq;
EC_Tot_data$Tot_eq[which(EC_Tot_data$Tot == 1)] <- 1;

# C: Defining system stability in Rr
EC_Tot_data$Dyn_state_f[which(EC_Tot_data$Tot_eq == 1)] <- paste("R", "eq", sep=".");
EC_Tot_data$Zone_f[which(EC_Tot_data$Dyn_state_f == paste("R", "eq", sep="."))] <- 1;

# Case 4: no species remaining
EC_Tot_data$Dyn_state_f[which(EC_Tot_data$Tot_eq == 0)] <- paste(0, "species", sep=".");
EC_Tot_data$Zone_f[which(EC_Tot_data$Dyn_state_f == paste(0, "species", sep="."))] <- 0;

#### Biodiversity change ====
for(i in Alpha){
  for(j in Beta){
    gamma <- i*j;
    if(gamma<0.9){next}else{
      subinv <- EC_Tot_data[which(EC_Tot_data$Alpha==i & EC_Tot_data$Beta==j),];
      subres <- Res[which(Res$Gamma==gamma),];
      EC_Tot_data$BioD[which(EC_Tot_data$Alpha == i & EC_Tot_data$Beta == j)] <- (subinv$Tot_eq - subres$Tot_eq);
    }
  }
}
#### Invasion states ====
EC_Tot_data$State <- as.character(EC_Tot_data$State);
EC_Tot_data$State[which(EC_Tot_data$Tot_eq == 3)] <- "Integration";
EC_Tot_data$State[which(EC_Tot_data$Tot_eq == 2 & EC_Tot_data$BioD == 0 & EC_Tot_data$Cr_eq < Tr)] <- "Substitution";
EC_Tot_data$State[which(EC_Tot_data$Tot_eq == 2 & EC_Tot_data$BioD == 0 & EC_Tot_data$Ci_eq < Tr)] <- "Resistance";
EC_Tot_data$State[which(EC_Tot_data$Tot_eq == 1 & EC_Tot_data$BioD == 0 )] <- "Resistance";
EC_Tot_data$State[which(EC_Tot_data$Tot_eq == 0 & EC_Tot_data$BioD == 0 )] <- "Resistance";

EC_Tot_data$State[which(EC_Tot_data$BioD < 0)] <- "Vulnerability";
EC_Tot_data$State[which(EC_Tot_data$BioD > 0 & EC_Tot_data$Tot_eq == 2 & EC_Tot_data$Cr_eq < Tr)] <- "Occupancy";
EC_Tot_data$State[which(EC_Tot_data$BioD > 0 & EC_Tot_data$Tot_eq != 3 & EC_Tot_data$Ci_eq < Tr)] <- "Rescue";

#### Changes in system stability in the community (initially vs. after 5000 years) ====
# 0. Initial-Final configurations
conf_name <- c(); initeq <- c("Eq", "Cycles"); fineq <- c("Null", "Eq", "Cycles")
for(i in initeq){conf_name <- append(conf_name, paste(i, fineq, sep='.'))};
conf_code <- c(0:5);
Conf_tab <- cbind.data.frame(conf_name, conf_code);

summary(as.factor(EC_Tot_data$Zone_i))
summary(as.factor(EC_Tot_data$Zone_f))
# A. Prepare the different combinations of stability changes
init_eq_code <- c(1, 2, 3);
init_eq <- c(rep("Eq", 2), "Cycles");
final_eq_code <- c(0:3, 6:8);
final_eq <- c("Null", "Eq", rep( c("Eq", "Cycles"), 2), "Eq");

# Data of observe combinations
I_code <- c(); for(i in init_eq_code){I_code <- append(I_code, rep(i, length(final_eq)))};
F_code <- rep(final_eq_code, length(init_eq_code));

IF_eq <- c(); for(i in init_eq){IF_eq <- append(IF_eq, paste(i, final_eq, sep='.'))}
IF_code <- rep(0, length(IF_eq));
IF <- cbind.data.frame(IF_eq, IF_code, I_code, F_code);

IF$IF_eq <- as.character(IF$IF_eq);
for(i in 1:length(IF_eq)){ IF$IF_code[i] <- Conf_tab$conf_code[which(Conf_tab$conf_name==IF_eq[i])]};
# B. Attribute the AC_Tot_data's outcomes to the corresponding configurations EQ code and name
EC_Tot_data$Zone_i <- as.numeric(EC_Tot_data$Zone_i);
EC_Tot_data$Zone_f <- as.numeric(EC_Tot_data$Zone_f);
EC_Tot_data$EQ_transit_InitFin <- as.character(EC_Tot_data$EQ_transit_InitFin)
EC_Tot_data$Zone_EQtransit_InitFin <- as.numeric(EC_Tot_data$Zone_EQtransit_InitFin);
for(i in 1:dim(IF)[1]){
  EC_Tot_data$EQ_transit_InitFin[which(EC_Tot_data$Zone_i == IF$I_code[i] & EC_Tot_data$Zone_f==IF$F_code[i])] <- IF$IF_eq[i];
  EC_Tot_data$Zone_EQtransit_InitFin[which(EC_Tot_data$Zone_i == IF$I_code[i] & EC_Tot_data$Zone_f==IF$F_code[i])] <- IF$IF_code[i]
}
#### Changes in system stability between resident and invaded communities (both after 5000 years) ====
Res_zonef <- rep(0, dim(EC_Tot_data)[1]);
EC_Tot_data <- cbind.data.frame(EC_Tot_data, Res_zonef);
# 0. Resident-Invaded configurations
conf_name <- c(); initeq <- c("Null", "Eq", "Cycles");
for(i in initeq){conf_name <- append(conf_name, paste(i, initeq, sep='.'))};
conf_code <- c(0:8);
Conf_tab <- cbind.data.frame(conf_name, conf_code);
# A. Merge Res$Zone_f to invaded data frame
for(i in Alpha){
  for(j in Beta){
    gamma <- i*j;
    if(gamma<0.9){next}else{
      subinv <- subset(EC_Tot_data,EC_Tot_data$Alpha==i & EC_Tot_data$Beta==j);
      subres <- subset(Res,Res$Gamma==gamma);
      EC_Tot_data$Res_zonef[which(EC_Tot_data$Alpha == i & EC_Tot_data$Beta == j)] <- subres$Zone_f;
    }
  }
}
# B. Prepare the different combinations of stability changes
res_eq_code <- c(0, 1, 2, 3);
res_eq <- c("Null", rep("Eq", 2), "Cycles");
# EQ & Code for invaded community
inv_eq_code <- c(seq(0,3,by=1), seq(6,8,by=1));
inv_eq <- c("Null", "Eq", rep( c("Eq", "Cycles"), 2), "Eq");
# Data of observe combinations
Res_code <- c(); for(i in res_eq_code){Res_code <- append(Res_code, rep(i, length(inv_eq)))};
Inv_code <- rep(inv_eq_code, length(res_eq_code));

Inv_Res_eq <- c(); for(i in res_eq){Inv_Res_eq <- append(Inv_Res_eq, paste(i, inv_eq, sep='.'))}
Inv_Res_code <- rep(0, length(Inv_Res_eq));
Inv_Res <- cbind.data.frame(Inv_Res_eq, Inv_Res_code, Res_code, Inv_code);

for(i in 1:length(Inv_Res_eq)){ Inv_Res$Inv_Res_code[i] <- Conf_tab$conf_code[which(Conf_tab$conf_name==Inv_Res_eq[i])]};
# C. Attribute the EC_Tot_data's outcomes to the corresponding configurations EQ code and name
EC_Tot_data$EQ_transit_ResInv <- as.character(EC_Tot_data$EQ_transit_ResInv)
EC_Tot_data$Zone_f <- as.numeric(EC_Tot_data$Zone_f);
EC_Tot_data$Res_zonef <- as.numeric(EC_Tot_data$Res_zonef);
EC_Tot_data$Zone_EQtransit_ResInv <- as.numeric(EC_Tot_data$Zone_EQtransit_ResInv);

for(i in 1:dim(Inv_Res)[1]){
  EC_Tot_data$EQ_transit_ResInv[which(EC_Tot_data$Res_zonef == Inv_Res$Res_code[i] & EC_Tot_data$Zone_f==Inv_Res$Inv_code[i])] <- Inv_Res$Inv_Res_eq[i];
  EC_Tot_data$Zone_EQtransit_ResInv[which(EC_Tot_data$Res_zonef == Inv_Res$Res_code[i] & EC_Tot_data$Zone_f==Inv_Res$Inv_code[i])] <- Inv_Res$Inv_Res_code[i]
}
#### Corrections on data made on 25/08/2021 Due to observed artifacts from transient analyses ====
# EC_Tot_data <- read.table("EC_Tot_data.txt", h=T,  sep = "\t");
# EC_Tot_data<-EC_Tot_data[,-c(42,43)]; #to remove code_state column
## Correction section 1: Beta < 1 & IS=Substitution + Zone_InvRes ====
dataobs <- EC_Tot_data[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Substitution' & EC_Tot_data$Dyn_state_f=="Ci.Rr.cycles"),];
PRdata <- as.data.frame(R.Rule.Expcomp(Temp=dataobs$Temp[1:dim(dataobs)[1]], Car=dataobs$Car[1:dim(dataobs)[1]], 
                                       MRr=dataobs$M_Rr[1:dim(dataobs)[1]], MCi=dataobs$M_Ci[1:dim(dataobs)[1]], MCr=dataobs$M_Cr[1:dim(dataobs)[1]]))
R_eqir_ratio  <- PRdata$R_star_Ci/PRdata$R_star_Cr;
R_eqir_diff   <- PRdata$R_star_Ci-PRdata$R_star_Cr;
Rr_feq <- ifelse(PRdata$R_star_Cr > 1e-12 | PRdata$R_star_Ci > 1e-12, 1, 0);
Ci_feq <- ifelse(PRdata$Ci_star > 1e-12, 1, 0);
Cr_feq <- ifelse(PRdata$Cr_star > 1e-12, 1, 0);
Tot_eq <- rep(0, length=dim(PRdata)[1]);
Rrstar <- rep(0, length=dim(PRdata)[1]);
Cistar <- rep(0, length=dim(PRdata)[1]);
Crstar <- rep(0, length=dim(PRdata)[1]);
PRdata <- cbind.data.frame(PRdata[,1:9], R_eqir_ratio, R_eqir_diff, Rr_feq, Ci_feq, Cr_feq, Tot_eq, Rrstar, Cistar, Crstar);

PRdata$Ci_feq[which(PRdata$R_eqir_ratio>1)] <- 0;
PRdata$Cr_feq[which(PRdata$R_eqir_ratio<1)] <- 0;
PRdata$Tot_eq <- apply(PRdata[,c(12:14)], 1, sum);

lambda_CrRr <- mapply(lambda_2sp_AC, Temp=PRdata$Temp[1:dim(PRdata)[1]], Car=PRdata$Car[1:dim(PRdata)[1]],
                       Bx=PRdata$R_star_Cr[1:dim(PRdata)[1]], Bz=PRdata$Cr_star[1:dim(PRdata)[1]],
                       Mx=PRdata$MRr[1:dim(PRdata)[1]], Mz=PRdata$MCr[1:dim(PRdata)[1]],
                       Gamma=1);

# Implemtenting the correction in the data now
EC_Tot_data$Rr_eq[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Substitution' & EC_Tot_data$Dyn_state_f=="Ci.Rr.cycles")] <- PRdata$R_star_Cr;
EC_Tot_data$Ci_eq[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Substitution' & EC_Tot_data$Dyn_state_f=="Ci.Rr.cycles")] <- 0;
EC_Tot_data$Cr_eq[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Substitution' & EC_Tot_data$Dyn_state_f=="Ci.Rr.cycles")] <- PRdata$Cr_star;
EC_Tot_data$Tot_eq[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Substitution' & EC_Tot_data$Dyn_state_f=="Ci.Rr.cycles")] <-PRdata$Tot_eq;
EC_Tot_data$Dom_Eigenvalue[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Substitution' & EC_Tot_data$Dyn_state_f=="Ci.Rr.cycles")] <- lambda_CrRr;
EC_Tot_data$Zone_f[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Substitution' & EC_Tot_data$Dyn_state_f=="Ci.Rr.cycles")] <- 3;
EC_Tot_data$Dyn_state_f[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Substitution' & EC_Tot_data$Dyn_state_f=="Ci.Rr.cycles")] <- paste("Cr", "Rr", "cycles", sep='.');
EC_Tot_data$State[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Substitution' & EC_Tot_data$Dyn_state_f=="Cr.Rr.cycles")] <- "Resistance";

## Correction section 2: Beta < 1 & Observed IS=Vulnerability, Bio=-1, InitFin==C.E ====
dataobs2 <- EC_Tot_data[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq"),];
PRdata <- as.data.frame(R.Rule.Expcomp(Temp=dataobs2$Temp[1:dim(dataobs2)[1]], Car=dataobs2$Car[1:dim(dataobs2)[1]], 
                                       MRr=dataobs2$M_Rr[1:dim(dataobs2)[1]], MCi=dataobs2$M_Ci[1:dim(dataobs2)[1]], MCr=dataobs2$M_Cr[1:dim(dataobs2)[1]]))

R_eqir_ratio  <- PRdata$R_star_Ci/PRdata$R_star_Cr;
R_eqir_diff   <- PRdata$R_star_Ci-PRdata$R_star_Cr;
Rr_feq <- ifelse(PRdata$R_star_Cr > 1e-12 | PRdata$R_star_Ci > 1e-12, 1, 0);
Ci_feq <- ifelse(PRdata$Ci_star > 1e-12, 1, 0);
Cr_feq <- ifelse(PRdata$Cr_star > 1e-12, 1, 0);
Tot_eq <- rep(0, length=dim(PRdata)[1]);
Rrstar <- rep(0, length=dim(PRdata)[1]);
Cistar <- rep(0, length=dim(PRdata)[1]);
Crstar <- rep(0, length=dim(PRdata)[1]);
PRdata <- cbind.data.frame(PRdata[,1:9], R_eqir_ratio, R_eqir_diff, Rr_feq, Ci_feq, Cr_feq, Tot_eq, Rrstar, Cistar, Crstar);

PRdata$Ci_feq[which(PRdata$R_eqir_ratio>1)] <- 0;
PRdata$Cr_feq[which(PRdata$R_eqir_ratio<1)] <- 0;
PRdata$Tot_eq <- apply(PRdata[,c(12:14)], 1, sum);

lambda_CrRr <- mapply(lambda_2sp_AC, Temp=PRdata$Temp[1:dim(PRdata)[1]], Car=PRdata$Car[1:dim(PRdata)[1]],
                      Bx=PRdata$R_star_Cr[1:dim(PRdata)[1]], Bz=PRdata$Cr_star[1:dim(PRdata)[1]],
                      Mx=PRdata$MRr[1:dim(PRdata)[1]], Mz=PRdata$MCr[1:dim(PRdata)[1]], Gamma=1);
# Implemtenting the correction in the data now
EC_Tot_data$Rr_eq[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq")] <- PRdata$R_star_Cr;
EC_Tot_data$Ci_eq[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq")] <- 0;
EC_Tot_data$Cr_eq[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq")] <- PRdata$Cr_star;
EC_Tot_data$Tot_eq[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq")] <-PRdata$Tot_eq;
EC_Tot_data$Dom_Eigenvalue[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq")] <- lambda_CrRr;

EC_Tot_data$BioD[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq")] <- 0;
EC_Tot_data$Zone_f[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq")] <- 3;

EC_Tot_data$EQ_transit_ResInv[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq")] <- paste("Cycles", "Cycles", sep='.');
EC_Tot_data$Zone_EQtransit_ResInv[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq")] <- 8
EC_Tot_data$EQ_transit_InitFin[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq")] <- paste("Cycles", "Cycles", sep='.');
EC_Tot_data$Zone_EQtransit_InitFin[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq")] <- 5;

EC_Tot_data$Dyn_state_f[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f=="R.eq")] <- paste("Cr", "Rr", "cycles", sep='.');
EC_Tot_data$State[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Vulnerability' & EC_Tot_data$Dyn_state_f==paste("Cr", "Rr", "cycles", sep='.'))] <- "Resistance";

## Correction section 3: Beta > 1 & Observed IS=Resistance, Cr.Rr.cycles ====
dataobs3 <- EC_Tot_data[which(EC_Tot_data$Beta>1 & EC_Tot_data$State == "Resistance" & EC_Tot_data$Dyn_state_f=='Cr.Rr.cycles'),];
PRdata <- as.data.frame(R.Rule.Expcomp(Temp=dataobs3$Temp[1:dim(dataobs3)[1]], Car=dataobs3$Car[1:dim(dataobs3)[1]], 
                                       MRr=dataobs3$M_Rr[1:dim(dataobs3)[1]], MCi=dataobs3$M_Ci[1:dim(dataobs3)[1]], MCr=dataobs3$M_Cr[1:dim(dataobs3)[1]]))

R_eqir_ratio  <- PRdata$R_star_Ci/PRdata$R_star_Cr;
R_eqir_diff   <- PRdata$R_star_Ci-PRdata$R_star_Cr;
Rr_feq <- ifelse(PRdata$R_star_Cr > 1e-12 | PRdata$R_star_Ci > 1e-12, 1, 0);
Ci_feq <- ifelse(PRdata$Ci_star > 1e-12, 1, 0);
Cr_feq <- ifelse(PRdata$Cr_star > 1e-12, 1, 0);
Tot_eq <- rep(0, length=dim(PRdata)[1]);
Rrstar <- rep(0, length=dim(PRdata)[1]);
Cistar <- rep(0, length=dim(PRdata)[1]);
Crstar <- rep(0, length=dim(PRdata)[1]);
PRdata <- cbind.data.frame(PRdata[,1:9], R_eqir_ratio, R_eqir_diff, Rr_feq, Ci_feq, Cr_feq, Tot_eq, Rrstar, Cistar, Crstar);

PRdata$Ci_feq[which(PRdata$R_eqir_ratio>1)] <- 0;
PRdata$Cr_feq[which(PRdata$R_eqir_ratio<1)] <- 0;
PRdata$Tot_eq <- apply(PRdata[,c(12:14)], 1, sum);

lambda_CiRr <- mapply(lambda_2sp_AB, Temp=PRdata$Temp[1:dim(PRdata)[1]], Car=PRdata$Car[1:dim(PRdata)[1]],
                      Bx=PRdata$R_star_Ci[1:dim(PRdata)[1]], By=PRdata$Ci_star[1:dim(PRdata)[1]],
                      Mx=PRdata$MRr[1:dim(PRdata)[1]], My=PRdata$MCi[1:dim(PRdata)[1]], Alpha=1);
# Implemtenting the correction in the data now
EC_Tot_data$Rr_eq[which(EC_Tot_data$Beta>1 & EC_Tot_data$State == "Resistance" & EC_Tot_data$Dyn_state_f=='Cr.Rr.cycles')] <- PRdata$R_star_Ci;
EC_Tot_data$Ci_eq[which(EC_Tot_data$Beta>1 & EC_Tot_data$State == "Resistance" & EC_Tot_data$Dyn_state_f=='Cr.Rr.cycles')] <- PRdata$Ci_star;
EC_Tot_data$Cr_eq[which(EC_Tot_data$Beta>1 & EC_Tot_data$State == "Resistance" & EC_Tot_data$Dyn_state_f=='Cr.Rr.cycles')] <- 0;
EC_Tot_data$Tot_eq[which(EC_Tot_data$Beta>1 & EC_Tot_data$State == "Resistance" & EC_Tot_data$Dyn_state_f=='Cr.Rr.cycles')] <-PRdata$Tot_eq;
EC_Tot_data$Dom_Eigenvalue[which(EC_Tot_data$Beta>1 & EC_Tot_data$State == "Resistance" & EC_Tot_data$Dyn_state_f=='Cr.Rr.cycles')] <- lambda_CiRr;
EC_Tot_data$Zone_f[which(EC_Tot_data$Beta>1 & EC_Tot_data$State == "Resistance" & EC_Tot_data$Dyn_state_f=='Cr.Rr.cycles')] <- 7;
EC_Tot_data$Dyn_state_f[which(EC_Tot_data$Beta>1 & EC_Tot_data$State == "Resistance" & EC_Tot_data$Dyn_state_f=='Cr.Rr.cycles')] <- 'Ci.Rr.cycles';
EC_Tot_data$State[which(EC_Tot_data$Beta>1 & EC_Tot_data$State == "Resistance" & EC_Tot_data$Dyn_state_f=='Ci.Rr.cycles')] <- "Substitution";

## Correction section 4: Beta < 1 & Observed IS=Resistance, Bio=0, InitFin==C.C ====
dataobs4 <- EC_Tot_data[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Resistance' & EC_Tot_data$BioD==0 & EC_Tot_data$EQ_transit_InitFin=="Cycles.Cycles"),];

for(i in 1:dim(dataobs4)[1]){
  if(dataobs4$Dom_Eigenvalue[i]>0){next
  }else{
    dataobs4$Dyn_state_f[i] <- "Cr.Rr.eq";
    dataobs4$Zone_f[i] <- 2;
    dataobs4$EQ_transit_InitFin[i] <- paste("Eq", "Eq", sep='.'); 
    dataobs4$Zone_EQtransit_InitFin[i] <- 2;
    dataobs4$EQ_transit_ResInv[i] <- paste("Eq", "Eq", sep='.'); 
    dataobs4$Zone_EQtransit_ResInv[i] <- 5;
    }
}
dataobs4 -> EC_Tot_data[which(EC_Tot_data$Beta<1 & EC_Tot_data$State=='Resistance' & EC_Tot_data$BioD==0 & EC_Tot_data$EQ_transit_InitFin=="Cycles.Cycles"),];

#### Savig the dataset ====
write.table(EC_Tot_data, "./data/EC_Tot_data.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
############# TC; Trophic chain ====
Inv <- read.table("./data/TC.data.txt", h=T,  sep = "\t");     # Data from TC module simulation
Res <- read.table("./data/CR_Tot_Alpha.txt", h=T,  sep = "\t");# Complete Data from 2sp C-R (A, B) including simulation outcomes and equilibrium analyses

Alpha <- c(1, 2, 5, 10);                          # BMR B/A
Beta <- c(1, 2, 5, 10);                           # BMR C/B
Gamma <- c(1, 2, 2.5, 4, 5, 10, 20, 25, 50, 100); # BMR B/A
Tr <- 1e-12;      # Extinction threshold

#### Initiating the dataset ====
InvTP <- rep(3, dim(Inv)[1]); # Trophic position of the invading species (here 3 = top predator)
Fw <- rep("TC", dim(Inv)[1]); # Food web module in which invasion occurs (here TC = Trophic chain)
Res_BMR <- rep('Alpha', length=dim(Inv)[1]);# species bmr of resident system

Dyn_state_i <- rep('NA', length=dim(Inv)[1]);# Dynamic state initial of resident species
Dyn_state_f <- rep('NA', length=dim(Inv)[1]);# Dynamic state final, after invasion
Zone_i <- rep(0, length=dim(Inv)[1]);     # initial, identical to 2spAB
Zone_f <- rep(0, length=dim(Inv)[1]);     # final
BioD <- rep(0, length=dim(Inv)[1]);       # Biodiversity change after invasion
State <- rep(0, length=dim(Inv)[1]);      # Invasion states
Rr_i <- rep(0, length=dim(Inv)[1]);       # Initial biomass density of resident resource (starting at C-R eq)
Cr_i <- rep(0, length=dim(Inv)[1]);       # Initial biomass density of resident consumer (starting at C-R eq)
Pi_i <- rep(1E-6, length=dim(Inv)[1]);    # Initial biomass density of invading predator

Rr_eq <- rep(0, length=dim(Inv)[1]);      # Resource biomass at eq (Equilibrium analyses)
Cr_eq <- rep(0, length=dim(Inv)[1]);      # Consumer biomass at eq (Equilibrium analyses)
Pi_eq <- rep(0, length=dim(Inv)[1]);      # Predator biomass at eq (Equilibrium analyses)
Tot_eq <- rep(0, length=dim(Inv)[1]);     # Diversity at eq (Equilibrium analyses)
Dom_Eigenvalue <- rep(0, length=dim(Inv)[1]); # Dominant eigenvalues for the different community structure
Inv <- cbind.data.frame(Inv, Rr_i, Cr_i, Pi_i, Dyn_state_i, Dyn_state_f, Zone_i, Zone_f, BioD, State);

EQ_transit_InitFin <- rep(0, length=dim(Inv)[1]); # Change in the system stability between the initial and final outcomes
Zone_EQtransit_InitFin <- rep(0, length=dim(Inv)[1]); # Corresponding codes of EQ_transit_if
EQ_transit_ResInv <- rep(0, length=dim(Inv)[1]); # Change in the system stability between the final outcomes of resident and invaded communities
Zone_EQtransit_ResInv <- rep(0, length=dim(Inv)[1]); # Corresponding codes of EQ_transit_invres

dynstati<-as.factor(c());
for(i in Alpha){
  for(j in Beta){
      subinv <- subset(Inv,Inv$Alpha==i & Inv$Beta==j);
      subres <- subset(Res,Res$Alpha==i);
      dynstati <- c(dynstati, subres$Dyn_state_i);
      Inv$Rr_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Rr_i;
      Inv$Cr_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Cr_i;
      Inv$Zone_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Zone_i;
  }
}
Inv$Dyn_state_i <- dynstati;

TC_Tot_data <- cbind.data.frame(
  Inv$Temperature, Inv$Carr, Fw, InvTP, Res_BMR, Inv$Alpha, Inv$Beta, Inv$Gamma, Inv$M_A, Inv$M_B, Inv$M_C,
  Inv$Rr_i, Inv$Cr_i, Inv$Pi_i, Inv$A, Inv$B, Inv$C, Inv$Afin, Inv$Bfin, Inv$Cfin, Inv$Tot, 
  Inv$ExtTimeA, Inv$ExtTimeB, Inv$ExtTimeC, Inv$Time, 
  Rr_eq, Cr_eq, Pi_eq, Dom_Eigenvalue, Tot_eq,
  Inv$Dyn_state_i, Inv$Dyn_state_f, Inv$Zone_i, Inv$Zone_f,
  EQ_transit_InitFin, Zone_EQtransit_InitFin, EQ_transit_ResInv, Zone_EQtransit_ResInv, Inv$BioD, Inv$State);

colnames(TC_Tot_data)<-c("Temp", "Car", "Fw", "Inv.position", "Res_BMR", "Alpha", "Beta", "Gamma", "M_Rr", "M_Cr", "M_Pi",
                         "Rr_i", "Cr_i", "Pi_i", "Rr_f", "Cr_f", "Pi_f", "Rr_fin", "Cr_fin", "Pi_fin", "Tot",
                         "Rr_extime", "Cr_extime", "Pi_extime", "Sim_duration",
                         "Rr_eq", "Cr_eq", "Pi_eq", "Dom_Eigenvalue", "Tot_eq",
                         "Dyn_state_i", "Dyn_state_f", "Zone_i", "Zone_f", 
                         "EQ_transit_InitFin", "Zone_EQtransit_InitFin", "EQ_transit_ResInv", "Zone_EQtransit_ResInv", "BioD", "State");

#### Species composition -> Equilibrium analyses & Calculus of Eigenvalues ====
## 1. Cross-check of species diversity between transient analyses (after 5000 years) and the equilibrium analyses ====
# a: Three species remaining at the end of the Transient analysis
data_3sp <- TC_Tot_data[which(TC_Tot_data$Tot == 3),]
data_3speq <- as.data.frame(EQUI_3sp(data_3sp$Temp, data_3sp$Car, data_3sp$Alpha, data_3sp$Beta, data_3sp$Gamma))
# Correction of species biomass at eq (if negative or < Tr)
data_3speq$Bx[which(data_3speq$Bx<Tr)] <- 0;
data_3speq$By[which(data_3speq$By<Tr)] <- 0;
data_3speq$Bz[which(data_3speq$Bz<Tr)] <- 0;
data_3speq$By[which(data_3speq$Bx<Tr)] <- 0; data_3speq$Bz[which(data_3speq$Bx<Tr)] <- 0;
data_3speq$Bz[which(data_3speq$By<Tr)] <- 0;

R <- ifelse(data_3speq$Bx>Tr,1,0); C <- ifelse(data_3speq$By>Tr,1,0); P <- ifelse(data_3speq$Bz>Tr,1,0);
Toteq <- as.vector(apply(cbind.data.frame(R, C, P),1,sum));
data_3speq <- cbind.data.frame(data_3speq, R, C, P, Toteq);
TC_Tot_data$Rr_eq[which(TC_Tot_data$Tot == 3)] <- data_3sp$Rr_eq <- data_3speq$Bx;
TC_Tot_data$Cr_eq[which(TC_Tot_data$Tot == 3)] <- data_3sp$Cr_eq <- data_3speq$By;
TC_Tot_data$Pi_eq[which(TC_Tot_data$Tot == 3)] <- data_3sp$Pi_eq <- data_3speq$Bz;
TC_Tot_data$Tot_eq[which(TC_Tot_data$Tot == 3)] <- data_3sp$Tot_eq <- data_3speq$Toteq;
# b: Only resident species remain (Cr & Rr)
data_2sp <- TC_Tot_data[which(TC_Tot_data$Tot == 2|TC_Tot_data$Tot_eq == 2),];
data_2speq <- as.data.frame(EQUI_2spAB(data_2sp$Temp, data_2sp$Car, data_2sp$Alpha));
data_2speq$Bx[which(data_2speq$Bx<Tr)] <- 0;
data_2speq$By[which(data_2speq$Bx<Tr)] <- 0; data_2speq$By[which(data_2speq$By<Tr)] <- 0;

R <- ifelse(data_2speq$Bx>Tr,1,0); C <- ifelse(data_2speq$By>Tr,1,0);
Toteq <- as.vector(apply(cbind.data.frame(R, C),1,sum));
data_2speq <- cbind.data.frame(data_2speq, R, C, Toteq);
TC_Tot_data$Rr_eq[which(TC_Tot_data$Tot == 2|TC_Tot_data$Tot_eq == 2)] <- data_2sp$Rr_eq <- data_2speq$Bx;
TC_Tot_data$Cr_eq[which(TC_Tot_data$Tot == 2|TC_Tot_data$Tot_eq == 2)] <- data_2sp$Cr_eq <- data_2speq$By;
TC_Tot_data$Tot_eq[which(TC_Tot_data$Tot == 2|TC_Tot_data$Tot_eq == 2)] <- data_2sp$Tot_eq <- data_2speq$Toteq;
# c: only resource remain, considered as R.eq
data_1sp <- TC_Tot_data[which(TC_Tot_data$Tot == 1|TC_Tot_data$Tot_eq == 1),];
data_1sp$Rr_eq <- K(Temp=data_1sp$Temp, Mx=data_1sp$M_Rr, Car=data_1sp$Car);

TC_Tot_data$Rr_eq[which(TC_Tot_data$Tot == 1|TC_Tot_data$Tot_eq == 1)] <- data_1sp$Rr_eq;
TC_Tot_data$Tot_eq[which(TC_Tot_data$Tot == 1|TC_Tot_data$Tot_eq == 1)] <- 1;

## 2. Calculus of Eigenvalues & Defining species composition ====
TC_Tot_data$Dyn_state_f <- as.character(TC_Tot_data$Dyn_state_f);
# a: P-C-R
data_3sp <- TC_Tot_data[which(TC_Tot_data$Tot_eq == 3),]
lambda_3sp <- mapply(lambda_3sp_TC, Temp=data_3sp$Temp[1:dim(data_3sp)[1]], Car=data_3sp$Car[1:dim(data_3sp)[1]],
                     Bx=data_3sp$Rr_eq[1:dim(data_3sp)[1]], By=data_3sp$Cr_eq[1:dim(data_3sp)[1]], Bz=data_3sp$Pi_eq[1:dim(data_3sp)[1]],
                     Mx=data_3sp$M_Rr[1:dim(data_3sp)[1]], My=data_3sp$M_Cr[1:dim(data_3sp)[1]], Mz=data_3sp$M_Pi[1:dim(data_3sp)[1]],
                     Alpha=data_3sp$Alpha[1:dim(data_3sp)[1]], Beta=data_3sp$Beta[1:dim(data_3sp)[1]], Gamma=data_3sp$Gamma[1:dim(data_3sp)[1]]);
data_3sp$Dom_Eigenvalue <- lambda_3sp;
data_3sp$Dyn_state_f[which(data_3sp$Dom_Eigenvalue > 0)] <- paste(3, "sp", "cycles", sep='.');
data_3sp$Dyn_state_f[which(data_3sp$Dom_Eigenvalue == 0 | data_3sp$Dom_Eigenvalue < 0)] <- paste(3, "sp", "eq", sep='.');
data_3sp$Zone_f[which(data_3sp$Dyn_state_f == paste(3, "sp", "eq", sep='.') )] <- 8;
data_3sp$Zone_f[which(data_3sp$Dyn_state_f == paste(3, "sp", "cycles", sep='.') )] <- 9;

TC_Tot_data$Dom_Eigenvalue[which(TC_Tot_data$Tot_eq == 3)] <- data_3sp$Dom_Eigenvalue;
TC_Tot_data$Dyn_state_f[which(TC_Tot_data$Tot_eq == 3)] <- data_3sp$Dyn_state_f;
TC_Tot_data$Zone_f[which(TC_Tot_data$Tot_eq == 3)] <- data_3sp$Zone_f;
# b: C-R
data_CrRr <- TC_Tot_data[which(TC_Tot_data$Tot_eq == 2),];
lambda_CrRr <- mapply(lambda_2sp_AB, Temp=data_CrRr$Temp[1:dim(data_CrRr)[1]], Car=data_CrRr$Car[1:dim(data_CrRr)[1]],
                      Bx=data_CrRr$Rr_eq[1:dim(data_CrRr)[1]], By=data_CrRr$Cr_eq[1:dim(data_CrRr)[1]],
                      Mx=data_CrRr$M_Rr[1:dim(data_CrRr)[1]], My=data_CrRr$M_Cr[1:dim(data_CrRr)[1]],
                      Alpha=data_CrRr$Alpha[1:dim(data_CrRr)[1]]);
data_CrRr$Dom_Eigenvalue <- lambda_CrRr;
data_CrRr$Dyn_state_f[which(data_CrRr$Dom_Eigenvalue > 0)] <- paste("Cr", "Rr", "cycles", sep='.'); 
data_CrRr$Dyn_state_f[which(data_CrRr$Dom_Eigenvalue == 0 | data_CrRr$Dom_Eigenvalue < 0)] <- paste("Cr", "Rr", "eq", sep='.');
data_CrRr$Zone_f[which(data_CrRr$Dyn_state_f == paste("Cr", "Rr", "eq", sep='.') )] <- 2;
data_CrRr$Zone_f[which(data_CrRr$Dyn_state_f == paste("Cr", "Rr", "cycles", sep='.') )] <- 3; 

TC_Tot_data$Dom_Eigenvalue[which(TC_Tot_data$Tot_eq == 2)] <- data_CrRr$Dom_Eigenvalue;
TC_Tot_data$Dyn_state_f[which(TC_Tot_data$Tot_eq == 2)] <- data_CrRr$Dyn_state_f;
TC_Tot_data$Zone_f[which(TC_Tot_data$Tot_eq == 2)] <- data_CrRr$Zone_f;
# c: R alone
TC_Tot_data$Dyn_state_f[which(TC_Tot_data$Tot_eq == 1)] <- paste("R", "eq", sep=".");
TC_Tot_data$Zone_f[which(TC_Tot_data$Dyn_state_f == paste("R", "eq", sep="."))] <- 1;

# cbis: correcting bugs
TC_Tot_data$Rr_eq[which(TC_Tot_data$Rr_eq=='NaN')] <- 0;
TC_Tot_data$Cr_eq[which(TC_Tot_data$Rr_eq=='NaN')] <- 0;
TC_Tot_data$Pi_eq[which(TC_Tot_data$Rr_eq=='NaN')] <- 0;
TC_Tot_data$Tot_eq[which(TC_Tot_data$Rr_eq=='NaN')] <- 0;
# d: 0 remaining species
TC_Tot_data$Dyn_state_f[which(TC_Tot_data$Tot_eq == 0)] <- paste(0, "species", sep=".");
TC_Tot_data$Zone_f[which(TC_Tot_data$Dyn_state_f == paste(0, "species", sep=".") )] <- 0;

#### Biodiversity change ====
for(i in Alpha){
  for(j in Beta){
    subinv <- TC_Tot_data[which(TC_Tot_data$Alpha==i & TC_Tot_data$Beta==j),];
    subres <- Res[which(Res$Alpha==i),];
    TC_Tot_data$BioD[which(TC_Tot_data$Alpha == i & TC_Tot_data$Beta == j)] <- (subinv$Tot_eq - subres$Tot_eq);
  }
}
#### Invasion states ====
TC_Tot_data$State <- as.character(TC_Tot_data$State);
TC_Tot_data$State[which(TC_Tot_data$Tot_eq == 3)] <- "Integration";
TC_Tot_data$State[which(TC_Tot_data$BioD < 0)] <- "Vulnerability";
TC_Tot_data$State[which(TC_Tot_data$BioD == 0 & TC_Tot_data$Pi_eq < Tr)] <- "Resistance";
TC_Tot_data$State[which(TC_Tot_data$BioD > 0 & TC_Tot_data$Pi_eq < Tr & TC_Tot_data$Tot_eq != 3)] <- "Rescue";
#### Changes in system stability in the community (initially vs. after 5000 years) ====
# 0. Initial-Final configurations
conf_name <- c(); initeq <- c("Eq", "Cycles"); fineq <- c("Null", "Eq", "Cycles")
for(i in initeq){conf_name <- append(conf_name, paste(i, fineq, sep='.'))};
conf_code <- c(0:5);
Conf_tab <- cbind.data.frame(conf_name, conf_code);

summary(as.factor(TC_Tot_data$Zone_i))
summary(as.factor(TC_Tot_data$Zone_f))
# A. Prepare the different combinations of stability changes
init_eq_code <- c(1, 2, 3);
init_eq <- c(rep("Eq", 2), "Cycles");
final_eq_code <- c(0:3, 8:9);
final_eq <- c("Null", "Eq", rep( c("Eq", "Cycles"), 2));

# Data of observe combinations
I_code <- c(); for(i in init_eq_code){I_code <- append(I_code, rep(i, length(final_eq)))};
F_code <- rep(final_eq_code, length(init_eq_code));

IF_eq <- c(); for(i in init_eq){IF_eq <- append(IF_eq, paste(i, final_eq, sep='.'))}
IF_code <- rep(0, length(IF_eq));
IF <- cbind.data.frame(IF_eq, IF_code, I_code, F_code);

IF$IF_eq <- as.character(IF$IF_eq);
for(i in 1:length(IF_eq)){ IF$IF_code[i] <- Conf_tab$conf_code[which(Conf_tab$conf_name==IF_eq[i])]};
# B. Attribute the AC_Tot_data's outcomes to the corresponding configurations EQ code and name
TC_Tot_data$Zone_i <- as.numeric(TC_Tot_data$Zone_i);
TC_Tot_data$Zone_f <- as.numeric(TC_Tot_data$Zone_f);
TC_Tot_data$EQ_transit_InitFin <- as.character(TC_Tot_data$EQ_transit_InitFin)
TC_Tot_data$Zone_EQtransit_InitFin <- as.numeric(TC_Tot_data$Zone_EQtransit_InitFin);
for(i in 1:dim(IF)[1]){
  TC_Tot_data$EQ_transit_InitFin[which(TC_Tot_data$Zone_i == IF$I_code[i] & TC_Tot_data$Zone_f==IF$F_code[i])] <- IF$IF_eq[i];
  TC_Tot_data$Zone_EQtransit_InitFin[which(TC_Tot_data$Zone_i == IF$I_code[i] & TC_Tot_data$Zone_f==IF$F_code[i])] <- IF$IF_code[i]
}
#### Changes in system stability between resident and invaded communities (both after 5000 years) ====
Res_zonef <- rep(0, dim(TC_Tot_data)[1]);
TC_Tot_data <- cbind.data.frame(TC_Tot_data, Res_zonef);
# 0. Resident-Invaded configurations
conf_name <- c(); initeq <- c("Null", "Eq", "Cycles");
for(i in initeq){conf_name <- append(conf_name, paste(i, initeq, sep='.'))};
conf_code <- c(0:8);
Conf_tab <- cbind.data.frame(conf_name, conf_code);

# A. Merge Res$Zone_f to invaded data frame
for(i in Alpha){
  for(j in Beta){
    subinv <- TC_Tot_data[which(TC_Tot_data$Alpha==i & TC_Tot_data$Beta==j),];
    subres <- Res[which(Res$Alpha==i),];
    TC_Tot_data$Res_zonef[which(TC_Tot_data$Alpha == i & TC_Tot_data$Beta == j)] <- subres$Zone_f;
  }
}
# B. Prepare the different combinations of stability changes
res_eq_code <- c(0, 1, 2, 3);
res_eq <- c("Null", rep("Eq", 2), "Cycles");
# EQ & Code for invaded community
inv_eq_code <- c(seq(0,3,by=1), 8:9);
inv_eq <- c("Null", "Eq", rep( c("Eq", "Cycles"), 2));
# Data of observe combinations
Res_code <- c(); for(i in res_eq_code){Res_code <- append(Res_code, rep(i, length(inv_eq)))};
Inv_code <- rep(inv_eq_code, length(res_eq_code));

Inv_Res_eq <- c(); for(i in res_eq){Inv_Res_eq <- append(Inv_Res_eq, paste(i, inv_eq, sep='.'))}
Inv_Res_code <- rep(0, length(Inv_Res_eq));
Inv_Res <- cbind.data.frame(Inv_Res_eq, Inv_Res_code, Res_code, Inv_code);

for(i in 1:length(Inv_Res_eq)){ Inv_Res$Inv_Res_code[i] <- Conf_tab$conf_code[which(Conf_tab$conf_name==Inv_Res_eq[i])]};
# C. Attribute the AC_Tot_data's outcomes to the corresponding configurations EQ code and name
TC_Tot_data$EQ_transit_ResInv <- as.character(TC_Tot_data$EQ_transit_ResInv)
TC_Tot_data$Zone_f <- as.numeric(TC_Tot_data$Zone_f);
TC_Tot_data$Res_zonef <- as.numeric(TC_Tot_data$Res_zonef);
TC_Tot_data$Zone_EQtransit_ResInv <- as.numeric(TC_Tot_data$Zone_EQtransit_ResInv);

for(i in 1:dim(Inv_Res)[1]){
  TC_Tot_data$EQ_transit_ResInv[which(TC_Tot_data$Res_zonef == Inv_Res$Res_code[i] & TC_Tot_data$Zone_f==Inv_Res$Inv_code[i])] <- Inv_Res$Inv_Res_eq[i];
  TC_Tot_data$Zone_EQtransit_ResInv[which(TC_Tot_data$Res_zonef == Inv_Res$Res_code[i] & TC_Tot_data$Zone_f==Inv_Res$Inv_code[i])] <- Inv_Res$Inv_Res_code[i]
}
#### Savig the dataset ====
write.table(TC_Tot_data, "./data/TC_Tot_data.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

############# IGPB; Intraguild Predation B (Intraguild-prey invasion) ====
Inv <- read.table("./data/IGPB.data.txt", h=T,  sep = "\t");     # Data from IGPB module simulation
Res <- read.table("./data/CR_Tot_Gamma.txt", h=T,  sep = "\t");  # Complete Data from 2sp C-R (A, C) including simulation outcomes and equilibrium analyses

Alpha <- c(1, 2, 5, 10);                          # BMR B/A
Beta <- c(1, 2, 5, 10);                           # BMR C/B
Gamma <- c(1, 2, 2.5, 4, 5, 10, 20, 25, 50, 100); # BMR B/A
Tr <- 1e-12;      # Extinction threshold

#### Initiating the dataset ====
InvTP <- rep(2, dim(Inv)[1]);   # Trophic position of the invading species (here 2 = intraguild-prey, a consumer)
Fw <- rep("IGPB", dim(Inv)[1]); # Food web module in which invasion occurs (here IGP = Intraguild predation, B for invading consumer)
Res_BMR <- rep('Gamma', length=dim(Inv)[1]);# species bmr of resident system

Dyn_state_i <- rep('NA', length=dim(Inv)[1]);# Dynamic state initial of resident species
Dyn_state_f <- rep('NA', length=dim(Inv)[1]);# Dynamic state final, after invasion
Zone_i <- rep(0, length=dim(Inv)[1]);     # initial, identical to 2spAC
Zone_f <- rep(0, length=dim(Inv)[1]);     # final
BioD <- rep(0, length=dim(Inv)[1]);       # Biodiversity change after invasion
State <- rep(0, length=dim(Inv)[1]);      # Invasion states
Rr_i <- rep(0, length=dim(Inv)[1]);       # Initial biomass density of resident resource (starting at C-R eq)
Pr_i <- rep(0, length=dim(Inv)[1]);       # Initial biomass density of resident predator (starting at C-R eq)
Ci_i <- rep(1E-6, length=dim(Inv)[1]);    # Initial biomass density of invading consumer
Rr_eq <- rep(0, length=dim(Inv)[1]);      # Resource biomass at eq (Equilibrium analyses)
Ci_eq <- rep(0, length=dim(Inv)[1]);      # Consumer biomass at eq (Equilibrium analyses)
Pr_eq <- rep(0, length=dim(Inv)[1]);      # Predator biomass at eq (Equilibrium analyses)
Tot_eq <- rep(0, length=dim(Inv)[1]);     # Diversity at eq (Equilibrium analyses)
Dom_Eigenvalue <- rep(0, length=dim(Inv)[1]); # Dominant eigenvalues for the different community structure
Inv <- cbind.data.frame(Inv, Rr_i, Ci_i, Pr_i, Dyn_state_i, Dyn_state_f, Zone_i, Zone_f, BioD, State);

EQ_transit_InitFin <- rep(0, length=dim(Inv)[1]); # Change in the system stability between the initial and final outcomes
Zone_EQtransit_InitFin <- rep(0, length=dim(Inv)[1]); # Corresponding codes of EQ_transit_if
EQ_transit_ResInv <- rep(0, length=dim(Inv)[1]); # Change in the system stability between the final outcomes of resident and invaded communities
Zone_EQtransit_ResInv <- rep(0, length=dim(Inv)[1]); # Corresponding codes of EQ_transit_invres

dynstati<-as.factor(c());
for(i in Alpha){
  for(j in Beta){
    gamma=i*j;
    subinv <- subset(Inv,Inv$Alpha==i & Inv$Beta==j);
    subres <- subset(Res,Res$Gamma==gamma);
    dynstati <- c(dynstati, subres$Dyn_state_i);
    Inv$Rr_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Rr_i;
    Inv$Cr_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Cr_i;
    Inv$Zone_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Zone_i;
  }
}
Inv$Dyn_state_i <- dynstati;

IGPB_Tot_data <- cbind.data.frame(
  Inv$Temperature, Inv$Carr, Fw, InvTP, Res_BMR, Inv$Alpha, Inv$Beta, Inv$Gamma, Inv$M_A, Inv$M_B, Inv$M_C,
  Inv$Rr_i, Inv$Ci_i, Inv$Pr_i, Inv$A, Inv$B, Inv$C, Inv$Afin, Inv$Bfin, Inv$Cfin, Inv$Tot, 
  Inv$ExtTimeA, Inv$ExtTimeB, Inv$ExtTimeC, Inv$Time, 
  Rr_eq, Ci_eq, Pr_eq, Dom_Eigenvalue, Tot_eq,
  Inv$Dyn_state_i, Inv$Dyn_state_f, Inv$Zone_i, Inv$Zone_f,
  EQ_transit_InitFin, Zone_EQtransit_InitFin, EQ_transit_ResInv, Zone_EQtransit_ResInv, Inv$BioD, Inv$State);

colnames(IGPB_Tot_data)<-c("Temp", "Car", "Fw", "Inv.position", "Res_BMR", "Alpha", "Beta", "Gamma", "M_Rr", "M_Ci", "M_Pr",
                         "Rr_i", "Ci_i", "Pr_i", "Rr_f", "Ci_f", "Pr_f", "Rr_fin", "Ci_fin", "Pr_fin", "Tot",
                         "Rr_extime", "Ci_extime", "Pr_extime", "Sim_duration", 
                         "Rr_eq", "Ci_eq", "Pr_eq", "Dom_Eigenvalue", "Tot_eq",
                         "Dyn_state_i", "Dyn_state_f", "Zone_i", "Zone_f",
                         "EQ_transit_InitFin", "Zone_EQtransit_InitFin", "EQ_transit_ResInv", "Zone_EQtransit_ResInv", "BioD", "State");

#### Species composition -> Equilibrium analyses & Calculus of Eigenvalues ====
## 1. Cross-check of species diversity between transient analyses (after 5000 years) and the equilibrium analyses
## 2. Calculus of Eigenvalues & Defining species composition
# Case 1: All species are present ###### To correct once we have the isoclines for 3sp
data_3sp <- IGPB_Tot_data[which(IGPB_Tot_data$Tot == 3),];

IGPB_Tot_data$Rr_eq[which(IGPB_Tot_data$Tot == 3)] <- data_3sp$Rr_eq <- data_3sp$Rr_f;
IGPB_Tot_data$Ci_eq[which(IGPB_Tot_data$Tot == 3)] <- data_3sp$Ci_eq <- data_3sp$Ci_f;
IGPB_Tot_data$Pr_eq[which(IGPB_Tot_data$Tot == 3)] <- data_3sp$Pr_eq <- data_3sp$Pr_f;
IGPB_Tot_data$Tot_eq[which(IGPB_Tot_data$Tot == 3)]<- data_3sp$Tot_eq <- 3;
IGPB_Tot_data$Dom_Eigenvalue[which(IGPB_Tot_data$Tot == 3)] <- data_3sp$Tot_eq <- NA;

#data_3sp$Dyn_state_f[which(data_3sp$lambda_3sp == 0 | data_3sp$lambda_3sp < 0)] <- paste(3, "sp", "eq", sep='.');
IGPB_Tot_data$Dyn_state_f[which(IGPB_Tot_data$Tot_eq == 3)] <- paste(3, "sp", "eq", sep='.');
IGPB_Tot_data$Zone_f[which(IGPB_Tot_data$Dyn_state_f == paste(3, "sp", "eq", sep='.') )] <- 8;
# Case 2: Invading consumer replaced the resident top predator species (Ci & Rr)
data_CiRr <- IGPB_Tot_data[which(IGPB_Tot_data$Tot == 2 & IGPB_Tot_data$Pr_fin == 0),];
data_CiRreq <- as.data.frame(EQUI_2spAB(Temp=data_CiRr$Temp, Car=data_CiRr$Car, Alpha=data_CiRr$Alpha));

data_CiRreq$By <- ifelse(data_CiRreq$By>Tr, data_CiRreq$By, 0);
data_CiRreq$Bx <- ifelse(data_CiRreq$By>Tr, data_CiRreq$Bx, K(Temp=data_CiRreq$Temp, Mx=data_CiRreq$Mx, Car=data_CiRreq$Car));

R <- ifelse(data_CiRreq$Bx > Tr, 1, 0); C <- ifelse(data_CiRreq$By > Tr, 1, 0);
Toteq <- as.vector(apply(cbind.data.frame(R, C),1,sum));
data_CiRreq <- cbind.data.frame(data_CiRreq, Toteq);

IGPB_Tot_data$Rr_eq[which(IGPB_Tot_data$Tot == 2 & IGPB_Tot_data$Pr_fin == 0 )] <- data_CiRreq$Bx;
IGPB_Tot_data$Ci_eq[which(IGPB_Tot_data$Tot == 2 & IGPB_Tot_data$Pr_fin == 0 )] <- data_CiRreq$By;
IGPB_Tot_data$Tot_eq[which(IGPB_Tot_data$Tot == 2 & IGPB_Tot_data$Pr_fin == 0 )] <-data_CiRreq$Toteq;
# Only taking 2 species at equilibria
data_CiRrB <- IGPB_Tot_data[which(IGPB_Tot_data$Tot_eq == 2 & IGPB_Tot_data$Pr_fin == 0),];
lambda_CiRrB <- mapply(lambda_2sp_AB, Temp=data_CiRrB$Temp[1:dim(data_CiRrB)[1]], Car=data_CiRrB$Car[1:dim(data_CiRrB)[1]],
                      Bx=data_CiRrB$Rr_eq[1:dim(data_CiRrB)[1]], By=data_CiRrB$Ci_eq[1:dim(data_CiRrB)[1]],
                      Mx=data_CiRrB$M_Rr[1:dim(data_CiRrB)[1]], My=data_CiRrB$M_Ci[1:dim(data_CiRrB)[1]],
                      Alpha=data_CiRrB$Alpha[1:dim(data_CiRrB)[1]]);
data_CiRrB <- cbind.data.frame(data_CiRrB, lambda_CiRrB);
data_CiRrB$Dyn_state_f[which(data_CiRrB$lambda_CiRrB > 0)] <- paste("Ci", "Rr", "cycles", sep='.');
data_CiRrB$Dyn_state_f[which(data_CiRrB$lambda_CiRrB == 0 | data_CiRrB$lambda_CiRrB < 0)] <- paste("Ci", "Rr", "eq", sep='.');

IGPB_Tot_data$Dyn_state_f[which(IGPB_Tot_data$Tot_eq == 2 & IGPB_Tot_data$Pr_fin == 0 )] <- data_CiRrB$Dyn_state_f;
IGPB_Tot_data$Zone_f[which(IGPB_Tot_data$Dyn_state_f == paste("Ci", "Rr", "eq", sep='.') )] <- 6;
IGPB_Tot_data$Zone_f[which(IGPB_Tot_data$Dyn_state_f == paste("Ci", "Rr", "cycles", sep='.') )] <- 7;
# Case 3: Only resident species remain (Pr & Rr)
data_CrRr <- IGPB_Tot_data[which(IGPB_Tot_data$Tot == 2 & IGPB_Tot_data$Ci_fin == 0 ),];
data_CrRreq <- as.data.frame(EQUI_2spAC(Temp=data_CrRr$Temp, Car=data_CrRr$Car, Gamma=data_CrRr$Gamma));

lambda_CrRr <- mapply(lambda_2sp_AC, Temp=data_CrRreq$Temp[1:dim(data_CrRreq)[1]], Car=data_CrRreq$Car[1:dim(data_CrRreq)[1]],
                      Bx=data_CrRreq$Bx[1:dim(data_CrRreq)[1]], Bz=data_CrRreq$Bz[1:dim(data_CrRreq)[1]],
                      Mx=data_CrRreq$Mx[1:dim(data_CrRreq)[1]], Mz=data_CrRreq$Mz[1:dim(data_CrRreq)[1]],
                      Gamma=data_CrRreq$Gamma[1:dim(data_CrRreq)[1]]);
data_CrRr <- cbind.data.frame(data_CrRr, lambda_CrRr);

IGPB_Tot_data$Rr_eq[which(IGPB_Tot_data$Tot == 2 & IGPB_Tot_data$Ci_fin == 0 )] <- data_CrRreq$Bx;
IGPB_Tot_data$Pr_eq[which(IGPB_Tot_data$Tot == 2 & IGPB_Tot_data$Ci_fin == 0 )] <- data_CrRreq$Bz;
IGPB_Tot_data$Tot_eq[which(IGPB_Tot_data$Tot == 2 & IGPB_Tot_data$Ci_fin == 0 )] <- 2;
IGPB_Tot_data$Dom_Eigenvalue[which(IGPB_Tot_data$Tot == 2 & IGPB_Tot_data$Ci_fin == 0 )] <- lambda_CrRr;

IGPB_Tot_data$Dyn_state_f[which(IGPB_Tot_data$Tot == 2 & IGPB_Tot_data$Ci_fin == 0 & IGPB_Tot_data$Dom_Eigenvalue > 0)] <- paste("Cr", "Rr", "cycles", sep='.');
IGPB_Tot_data$Dyn_state_f[which(IGPB_Tot_data$Tot == 2 & IGPB_Tot_data$Ci_fin == 0 & (IGPB_Tot_data$Dom_Eigenvalue == 0 | IGPB_Tot_data$Dom_Eigenvalue < 0) )] <- paste("Cr", "Rr", "eq", sep='.');
IGPB_Tot_data$Zone_f[which(IGPB_Tot_data$Dyn_state_f == paste("Cr", "Rr", "eq", sep='.') )] <- 2;
IGPB_Tot_data$Zone_f[which(IGPB_Tot_data$Dyn_state_f == paste("Cr", "Rr", "cycles", sep='.') )] <- 3;
# Case 4: only resource remain, considered as R.eq
dataR <- IGPB_Tot_data[which(IGPB_Tot_data$Tot == 1),];
IGPB_Tot_data$Tot_eq[which(IGPB_Tot_data$Tot == 1)] <- 1;
IGPB_Tot_data$Rr_eq[which(IGPB_Tot_data$Tot == 1)] <- K(Temp=dataR$Temp, Mx=dataR$M_Rr, Car=dataR$Car);
IGPB_Tot_data$Dyn_state_f[which(IGPB_Tot_data$Tot_eq == 1)] <- paste("R", "eq", sep=".");
IGPB_Tot_data$Zone_f[which(IGPB_Tot_data$Tot_eq == 1)] <- 1;
# Case 0: no species remaining
IGPB_Tot_data$Dyn_state_f[which(IGPB_Tot_data$Tot_eq == 0)] <- paste(0, "species", sep=".");
IGPB_Tot_data$Zone_f[which(IGPB_Tot_data$Tot_eq == 0)] <- 0;
#### Biodiversity change ====
for(i in Alpha){
  for(j in Beta){
    gamma=i*j;
    subinv <- IGPB_Tot_data[which(IGPB_Tot_data$Alpha==i & IGPB_Tot_data$Beta==j),];
    subres <- Res[which(Res$Gamma==gamma),];
    IGPB_Tot_data$BioD[which(IGPB_Tot_data$Alpha == i & IGPB_Tot_data$Beta == j)] <- (subinv$Tot_eq - subres$Tot_eq);
  }
}
#### Invasion states ====
IGPB_Tot_data$State[which(IGPB_Tot_data$Tot_eq == 3)] <- "Integration";
IGPB_Tot_data$State[which(IGPB_Tot_data$BioD < 0)] <- "Vulnerability";
IGPB_Tot_data$State[which(IGPB_Tot_data$BioD == 0 & IGPB_Tot_data$Tot_eq == 0)]  <- "Resistance";
IGPB_Tot_data$State[which(IGPB_Tot_data$BioD == 0 & IGPB_Tot_data$Tot_eq == 1)]  <- "Resistance";
IGPB_Tot_data$State[which(IGPB_Tot_data$BioD == 0 & IGPB_Tot_data$Tot_eq==2 & IGPB_Tot_data$Pr_eq > Tr)] <- "Resistance";
IGPB_Tot_data$State[which(IGPB_Tot_data$BioD == 0 & IGPB_Tot_data$Tot_eq==2 & IGPB_Tot_data$Ci_eq > Tr)] <- "Substitution";
IGPB_Tot_data$State[which(IGPB_Tot_data$BioD > 0 & IGPB_Tot_data$Ci_eq ==0 & IGPB_Tot_data$Tot_eq != 3)] <- "Rescue";
IGPB_Tot_data$State[which(IGPB_Tot_data$BioD > 0 & IGPB_Tot_data$Ci_eq > Tr & IGPB_Tot_data$Tot_eq != 3)] <- "Occupancy";
#### Changes in system stability in the community (initially vs. after 5000 years) ====
# 0. Initial-Final configurations
conf_name <- c(); initeq <- c("Eq", "Cycles"); fineq <- c("Null", "Eq", "Cycles")
for(i in initeq){conf_name <- append(conf_name, paste(i, fineq, sep='.'))};
conf_code <- c(0:5);
Conf_tab <- cbind.data.frame(conf_name, conf_code);

summary(as.factor(IGPB_Tot_data$Zone_i))
summary(as.factor(IGPB_Tot_data$Zone_f))
# A. Prepare the different combinations of stability changes
init_eq_code <- c(1, 2, 3);
init_eq <- c(rep("Eq", 2), "Cycles");
final_eq_code <- c(0:3, 6:8);
final_eq <- c("Null", "Eq", rep( c("Eq", "Cycles"), 2), "Eq");

# Data of observe combinations
I_code <- c(); for(i in init_eq_code){I_code <- append(I_code, rep(i, length(final_eq)))};
F_code <- rep(final_eq_code, length(init_eq_code));

IF_eq <- c(); for(i in init_eq){IF_eq <- append(IF_eq, paste(i, final_eq, sep='.'))}
IF_code <- rep(0, length(IF_eq));
IF <- cbind.data.frame(IF_eq, IF_code, I_code, F_code);

IF$IF_eq <- as.character(IF$IF_eq);
for(i in 1:length(IF_eq)){ IF$IF_code[i] <- Conf_tab$conf_code[which(Conf_tab$conf_name==IF_eq[i])]};
# B. Attribute the AC_Tot_data's outcomes to the corresponding configurations EQ code and name
IGPB_Tot_data$Zone_i <- as.numeric(IGPB_Tot_data$Zone_i);
IGPB_Tot_data$Zone_f <- as.numeric(IGPB_Tot_data$Zone_f);
IGPB_Tot_data$EQ_transit_InitFin <- as.character(IGPB_Tot_data$EQ_transit_InitFin)
IGPB_Tot_data$Zone_EQtransit_InitFin <- as.numeric(IGPB_Tot_data$Zone_EQtransit_InitFin);
for(i in 1:dim(IF)[1]){
  IGPB_Tot_data$EQ_transit_InitFin[which(IGPB_Tot_data$Zone_i == IF$I_code[i] & IGPB_Tot_data$Zone_f==IF$F_code[i])] <- IF$IF_eq[i];
  IGPB_Tot_data$Zone_EQtransit_InitFin[which(IGPB_Tot_data$Zone_i == IF$I_code[i] & IGPB_Tot_data$Zone_f==IF$F_code[i])] <- IF$IF_code[i]
}
#### Changes in system stability between resident and invaded communities (both after 5000 years) ====
Res_zonef <- rep(0, dim(IGPB_Tot_data)[1]);
IGPB_Tot_data <- cbind.data.frame(IGPB_Tot_data, Res_zonef);
# 0. Resident-Invaded configurations
conf_name <- c(); initeq <- c("Null", "Eq", "Cycles");
for(i in initeq){conf_name <- append(conf_name, paste(i, initeq, sep='.'))};
conf_code <- c(0:8);
Conf_tab <- cbind.data.frame(conf_name, conf_code);

# A. Merge Res$Zone_f to invaded data frame
for(i in Alpha){
  for(j in Beta){
    gamma <- i*j;
    subinv <- IGPB_Tot_data[which(IGPB_Tot_data$Alpha==i & IGPB_Tot_data$Beta==j),];
    subres <- Res[which(Res$Gamma==gamma),];
    IGPB_Tot_data$Res_zonef[which(IGPB_Tot_data$Alpha == i & IGPB_Tot_data$Beta == j)] <- subres$Zone_f;
  }
}
# B. Prepare the different combinations of stability changes
res_eq_code <- c(0, 1, 2, 3);
res_eq <- c("Null", rep("Eq", 2), "Cycles");
# EQ & Code for invaded community
inv_eq_code <- c(0:3, 6:8);
inv_eq <- c("Null", "Eq", rep( c("Eq", "Cycles"), 2), "Eq");
# Data of observe combinations
Res_code <- c(); for(i in res_eq_code){Res_code <- append(Res_code, rep(i, length(inv_eq)))};
Inv_code <- rep(inv_eq_code, length(res_eq_code));

Inv_Res_eq <- c(); for(i in res_eq){Inv_Res_eq <- append(Inv_Res_eq, paste(i, inv_eq, sep='.'))}
Inv_Res_code <- rep(0, length(Inv_Res_eq));
Inv_Res <- cbind.data.frame(Inv_Res_eq, Inv_Res_code, Res_code, Inv_code);

for(i in 1:length(Inv_Res_eq)){ Inv_Res$Inv_Res_code[i] <- Conf_tab$conf_code[which(Conf_tab$conf_name==Inv_Res_eq[i])]};
# C. Attribute the AC_Tot_data's outcomes to the corresponding configurations EQ code and name
IGPB_Tot_data$EQ_transit_ResInv <- as.character(IGPB_Tot_data$EQ_transit_ResInv)
IGPB_Tot_data$Zone_f <- as.numeric(IGPB_Tot_data$Zone_f);
IGPB_Tot_data$Res_zonef <- as.numeric(IGPB_Tot_data$Res_zonef);
IGPB_Tot_data$Zone_EQtransit_ResInv <- as.numeric(IGPB_Tot_data$Zone_EQtransit_ResInv);

for(i in 1:dim(Inv_Res)[1]){
  IGPB_Tot_data$EQ_transit_ResInv[which(IGPB_Tot_data$Res_zonef == Inv_Res$Res_code[i] & IGPB_Tot_data$Zone_f==Inv_Res$Inv_code[i])] <- Inv_Res$Inv_Res_eq[i];
  IGPB_Tot_data$Zone_EQtransit_ResInv[which(IGPB_Tot_data$Res_zonef == Inv_Res$Res_code[i] & IGPB_Tot_data$Zone_f==Inv_Res$Inv_code[i])] <- Inv_Res$Inv_Res_code[i]
}
#### Savig the dataset ====
write.table(IGPB_Tot_data, "./data/IGPB_Tot_data.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

############# IGPC; Intraguild Predation C (Intraguild predator invasion) ====
Inv <- read.table("./data/IGPC.data.txt", h=T,  sep = "\t");     # Data from IGPC module simulation
Res <- read.table("./data/CR_Tot_Alpha.txt", h=T,  sep = "\t");  # Complete Data from 2sp C-R (A, B) including simulation outcomes and equilibrium analyses

Alpha <- c(1, 2, 5, 10);                          # BMR B/A
Beta <- c(1, 2, 5, 10);                           # BMR C/B
Gamma <- c(1, 2, 2.5, 4, 5, 10, 20, 25, 50, 100); # BMR B/A
Tr <- 1e-12;      # Extinction threshold

#### Initiating the dataset ====
InvTP <- rep(3, dim(Inv)[1]);   # Trophic position of the invading species (here 3 = intraguild-predator, a top predator)
Fw <- rep("IGPC", dim(Inv)[1]); # Food web module in which invasion occurs (here IGP = Intraguild predation, C for invading top predator)
Res_BMR <- rep('Alpha', length=dim(Inv)[1]);# species bmr of resident system

Dyn_state_i <- rep('NA', length=dim(Inv)[1]);# Dynamic state initial of resident species
Dyn_state_f <- rep('NA', length=dim(Inv)[1]);# Dynamic state final, after invasion
Zone_i <- rep(0, length=dim(Inv)[1]);     # initial, identical to 2spAC
Zone_f <- rep(0, length=dim(Inv)[1]);     # final
BioD <- rep(0, length=dim(Inv)[1]);       # Biodiversity change after invasion
State <- rep(0, length=dim(Inv)[1]);      # Invasion states
Rr_i <- rep(0, length=dim(Inv)[1]);       # Initial biomass density of resident resource (starting at C-R eq)
Cr_i <- rep(0, length=dim(Inv)[1]);       # Initial biomass density of resident consumer (starting at C-R eq)
Pi_i <- rep(1E-6, length=dim(Inv)[1]);    # Initial biomass density of invading top predator
Rr_eq <- rep(0, length=dim(Inv)[1]);      # Resource biomass at eq (Equilibrium analyses)
Cr_eq <- rep(0, length=dim(Inv)[1]);      # Consumer biomass at eq (Equilibrium analyses)
Pi_eq <- rep(0, length=dim(Inv)[1]);      # Predator biomass at eq (Equilibrium analyses)
Tot_eq <- rep(0, length=dim(Inv)[1]);     # Diversity at eq (Equilibrium analyses)
Dom_Eigenvalue <- rep(0, length=dim(Inv)[1]); # Dominant eigenvalues for the different community structure
Inv <- cbind.data.frame(Inv, Rr_i, Cr_i, Pi_i, Dyn_state_i, Dyn_state_f, Zone_i, Zone_f, BioD, State);

EQ_transit_InitFin <- rep(0, length=dim(Inv)[1]); # Change in the system stability between the initial and final outcomes
Zone_EQtransit_InitFin <- rep(0, length=dim(Inv)[1]); # Corresponding codes of EQ_transit_if
EQ_transit_ResInv <- rep(0, length=dim(Inv)[1]); # Change in the system stability between the final outcomes of resident and invaded communities
Zone_EQtransit_ResInv <- rep(0, length=dim(Inv)[1]); # Corresponding codes of EQ_transit_invres

dynstati<-as.factor(c());
for(i in Alpha){
  for(j in Beta){
    subinv <- subset(Inv,Inv$Alpha==i & Inv$Beta==j);
    subres <- subset(Res,Res$Alpha==i);
    dynstati <- c(dynstati, subres$Dyn_state_i);
    Inv$Rr_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Rr_i;
    Inv$Cr_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Cr_i;
    Inv$Zone_i[which(Inv$Alpha == i & Inv$Beta == j)] <- subres$Zone_i;
  }
}
Inv$Dyn_state_i <- dynstati;

IGPC_Tot_data <- cbind.data.frame(
  Inv$Temperature, Inv$Carr, Fw, InvTP, Res_BMR, Inv$Alpha, Inv$Beta, Inv$Gamma, Inv$M_A, Inv$M_B, Inv$M_C,
  Inv$Rr_i, Inv$Cr_i, Inv$Pi_i, Inv$A, Inv$B, Inv$C, Inv$Afin, Inv$Bfin, Inv$Cfin, Inv$Tot, 
  Inv$ExtTimeA, Inv$ExtTimeB, Inv$ExtTimeC, Inv$Time, 
  Rr_eq, Cr_eq, Pi_eq, Dom_Eigenvalue, Tot_eq,
  Inv$Dyn_state_i, Inv$Dyn_state_f, Inv$Zone_i, Inv$Zone_f,
  EQ_transit_InitFin, Zone_EQtransit_InitFin, EQ_transit_ResInv, Zone_EQtransit_ResInv, Inv$BioD, Inv$State);

colnames(IGPC_Tot_data)<-c("Temp", "Car", "Fw", "Inv.position", "Res_BMR", "Alpha", "Beta", "Gamma", "M_Rr", "M_Cr", "M_Pi",
                           "Rr_i", "Cr_i", "Pi_i", "Rr_f", "Cr_f", "Pi_f", "Rr_fin", "Cr_fin", "Pi_fin", "Tot",
                           "Rr_extime", "Cr_extime", "Pi_extime", "Sim_duration",
                           "Rr_eq", "Cr_eq", "Pi_eq", "Dom_Eigenvalue", "Tot_eq", "Dyn_state_i", "Dyn_state_f", "Zone_i", "Zone_f",
                           "EQ_transit_InitFin", "Zone_EQtransit_InitFin", "EQ_transit_ResInv", "Zone_EQtransit_ResInv", "BioD", "State");

#### Species composition -> Equilibrium analyses & Calculus of Eigenvalues ====
## 1. Cross-check of species diversity between transient analyses (after 5000 years) and the equilibrium analyses
## 2. Calculus of Eigenvalues & Defining species composition
IGPC_Tot_data$Dyn_state_f <- as.character(IGPC_Tot_data$Dyn_state_f);
# Case 1: All species are present ###### To correct once we have the isoclines for 3sp
data_3sp <- IGPC_Tot_data[which(IGPC_Tot_data$Tot == 3),];

IGPC_Tot_data$Rr_eq[which(IGPC_Tot_data$Tot == 3)] <- data_3sp$Rr_eq <- data_3sp$Rr_f;
IGPC_Tot_data$Cr_eq[which(IGPC_Tot_data$Tot == 3)] <- data_3sp$Cr_eq <- data_3sp$Cr_f;
IGPC_Tot_data$Pi_eq[which(IGPC_Tot_data$Tot == 3)] <- data_3sp$Pi_eq <- data_3sp$Pi_f;
IGPC_Tot_data$Tot_eq[which(IGPC_Tot_data$Tot == 3)]<- data_3sp$Tot_eq <- 3;
IGPC_Tot_data$Dom_Eigenvalue[which(IGPC_Tot_data$Tot == 3)] <- data_3sp$Tot_eq <- NA;

IGPC_Tot_data$Dyn_state_f[which(IGPC_Tot_data$Tot == 3)] <- paste(3, "sp", "eq", sep='.');
IGPC_Tot_data$Zone_f[which(IGPC_Tot_data$Dyn_state_f == paste(3, "sp", "eq", sep='.') )] <- 8;

# Case 2: Invading predator replaced the resident consumer species (Pi & Rr)
data_PiRr <- IGPC_Tot_data[which(IGPC_Tot_data$Tot == 2 & IGPC_Tot_data$Cr_fin == 0 ),];
data_PiRreq <- as.data.frame(EQUI_2spAC(Temp=data_PiRr$Temp, Car=data_PiRr$Car, Gamma=data_PiRr$Gamma));

R <- ifelse(data_PiRreq$Bx > Tr, 1, 0); P <- ifelse(data_PiRreq$Bz > Tr, 1, 0);
Toteq <- as.vector(apply(cbind.data.frame(R, P),1,sum));
data_PiRreq <- cbind.data.frame(data_PiRreq, Toteq);

lambda_PiRr <- mapply(lambda_2sp_AC, Temp=data_PiRreq$Temp[1:dim(data_PiRreq)[1]], Car=data_PiRreq$Car[1:dim(data_PiRreq)[1]],
                      Bx=data_PiRreq$Bx[1:dim(data_PiRreq)[1]], Bz=data_PiRreq$Bz[1:dim(data_PiRreq)[1]],
                      Mx=data_PiRreq$Mx[1:dim(data_PiRreq)[1]], Mz=data_PiRreq$Mz[1:dim(data_PiRreq)[1]],
                      Gamma=data_PiRreq$Gamma[1:dim(data_PiRreq)[1]]);
data_PiRr$Dom_Eigenvalue <- lambda_PiRr;
data_PiRr$Dyn_state_f[which(data_PiRr$Dom_Eigenvalue > 0)] <- paste("Ci", "Rr", "cycles", sep='.');
data_PiRr$Dyn_state_f[which(data_PiRr$Dom_Eigenvalue == 0 | data_PiRr$Dom_Eigenvalue < 0)] <- paste("Ci", "Rr", "eq", sep='.');

IGPC_Tot_data$Rr_eq[which(IGPC_Tot_data$Tot == 2 & IGPC_Tot_data$Cr_fin == 0 )] <- data_PiRreq$Bx
IGPC_Tot_data$Pi_eq[which(IGPC_Tot_data$Tot == 2 & IGPC_Tot_data$Cr_fin == 0 )] <- data_PiRreq$Bz
IGPC_Tot_data$Tot_eq[which(IGPC_Tot_data$Tot == 2 & IGPC_Tot_data$Cr_fin == 0 )] <- data_PiRreq$Toteq;
IGPC_Tot_data$Dom_Eigenvalue[which(IGPC_Tot_data$Tot == 2 & IGPC_Tot_data$Cr_fin == 0 )] <- lambda_PiRr;
IGPC_Tot_data$Dyn_state_f[which(IGPC_Tot_data$Tot == 2 & IGPC_Tot_data$Cr_fin == 0 )] <- data_PiRr$Dyn_state_f;

IGPC_Tot_data$Zone_f[which(IGPC_Tot_data$Dyn_state_f == paste("Ci", "Rr", "eq", sep='.') )] <- 6;
IGPC_Tot_data$Zone_f[which(IGPC_Tot_data$Dyn_state_f == paste("Ci", "Rr", "cycles", sep='.') )] <- 7;

# Case 3: Only resident species remain (Cr & Rr)
data_CrRr <- IGPC_Tot_data[which(IGPC_Tot_data$Tot == 2 & IGPC_Tot_data$Pi_fin == 0 ),];
data_CrRreq <- as.data.frame(EQUI_2spAB(Temp=data_CrRr$Temp, Car=data_CrRr$Car, Alpha=data_CrRr$Alpha));

R <- ifelse(data_CrRreq$Bx > Tr, 1, 0); C <- ifelse(data_CrRreq$By > Tr, 1, 0);
Toteq <- as.vector(apply(cbind.data.frame(R, C),1,sum));
data_CrRreq <- cbind.data.frame(data_CrRreq, Toteq);

lambda_CrRr <- mapply(lambda_2sp_AB, Temp=data_CrRreq$Temp[1:dim(data_CrRreq)[1]], Car=data_CrRreq$Car[1:dim(data_CrRreq)[1]],
                      Bx=data_CrRreq$Bx[1:dim(data_CrRreq)[1]], By=data_CrRreq$By[1:dim(data_CrRreq)[1]],
                      Mx=data_CrRreq$Mx[1:dim(data_CrRr)[1]], My=data_CrRreq$My[1:dim(data_CrRreq)[1]],
                      Alpha=data_CrRreq$Alpha[1:dim(data_CrRreq)[1]]);
data_CrRr$Dom_Eigenvalue <- lambda_CrRr;
data_CrRr$Dyn_state_f[which(data_CrRr$Dom_Eigenvalue > 0)] <- paste("Cr", "Rr", "cycles", sep='.');
data_CrRr$Dyn_state_f[which(data_CrRr$Dom_Eigenvalue == 0 | data_CrRr$Dom_Eigenvalue < 0)] <- paste("Cr", "Rr", "eq", sep='.');

IGPC_Tot_data$Rr_eq[which(IGPC_Tot_data$Tot == 2 & IGPC_Tot_data$Pi_fin == 0 )] <- data_CrRreq$Bx;
IGPC_Tot_data$Cr_eq[which(IGPC_Tot_data$Tot == 2 & IGPC_Tot_data$Pi_fin == 0 )] <- data_CrRreq$By;
IGPC_Tot_data$Tot_eq[which(IGPC_Tot_data$Tot == 2 & IGPC_Tot_data$Pi_fin == 0 )] <- data_CrRreq$Toteq;
IGPC_Tot_data$Dom_Eigenvalue[which(IGPC_Tot_data$Tot == 2 & IGPC_Tot_data$Pi_fin == 0 )] <- lambda_CrRr;
IGPC_Tot_data$Dyn_state_f[which(IGPC_Tot_data$Tot == 2 & IGPC_Tot_data$Pi_fin == 0 )] <- data_CrRr$Dyn_state_f;

IGPC_Tot_data$Zone_f[which(IGPC_Tot_data$Dyn_state_f == paste("Cr", "Rr", "eq", sep='.') )] <- 2;
IGPC_Tot_data$Zone_f[which(IGPC_Tot_data$Dyn_state_f == paste("Cr", "Rr", "cycles", sep='.') )] <- 3;

# Case 4: only resource remain, considered as R.eq
data1sp <- IGPC_Tot_data[which(IGPC_Tot_data$Tot == 1),];
data1sp$Rr_eq <- K(Temp=data1sp$Temp, Mx=data1sp$M_Rr, Car = data1sp$Car);
IGPC_Tot_data$Rr_eq[which(IGPC_Tot_data$Tot == 1)] <- data1sp$Rr_eq;
IGPC_Tot_data$Tot_eq[which(IGPC_Tot_data$Tot == 1)] <- 1;
IGPC_Tot_data$Dyn_state_f[which(IGPC_Tot_data$Tot == 1)] <- paste("R", "eq", sep=".");
IGPC_Tot_data$Zone_f[which(IGPC_Tot_data$Dyn_state_f == paste("R", "eq", sep="."))] <- 1;
# Case 0: no species remaining
IGPC_Tot_data$Dyn_state_f[which(IGPC_Tot_data$Tot_eq == 0)] <- paste(0, "species", sep=".");
IGPC_Tot_data$Zone_f[which(IGPC_Tot_data$Tot_eq == 0)] <- 0;
#### Biodiversity change ====
for(i in Alpha){
  subres <- Res[which(Res$Alpha==i),];
  for(j in Beta){
    subinv <- IGPC_Tot_data[which(IGPC_Tot_data$Alpha==i & IGPC_Tot_data$Beta==j),];
    IGPC_Tot_data$BioD[which(IGPC_Tot_data$Alpha == i & IGPC_Tot_data$Beta == j)] <- (subinv$Tot_eq - subres$Tot_eq);
  }
}
#### Invasion states ====
IGPC_Tot_data$State[which(IGPC_Tot_data$Tot_eq == 3)] <- "Integration";
IGPC_Tot_data$State[which(IGPC_Tot_data$BioD < 0)] <- "Vulnerability";
IGPC_Tot_data$State[which(IGPC_Tot_data$BioD == 0 & IGPC_Tot_data$Tot_eq == 0)] <- "Resistance";
IGPC_Tot_data$State[which(IGPC_Tot_data$BioD == 0 & IGPC_Tot_data$Tot_eq == 1)] <- "Resistance";
IGPC_Tot_data$State[which(IGPC_Tot_data$BioD == 0 & IGPC_Tot_data$Tot_eq == 2 & IGPC_Tot_data$Cr_eq > Tr)] <- "Resistance";
IGPC_Tot_data$State[which(IGPC_Tot_data$BioD == 0 & IGPC_Tot_data$Tot_eq == 2 & IGPC_Tot_data$Cr_eq < Tr)] <- "Substitution";
IGPC_Tot_data$State[which(IGPC_Tot_data$BioD > 0 & IGPC_Tot_data$Tot_eq != 3 & IGPC_Tot_data$Pi_eq < Tr)] <- "Rescue";
IGPC_Tot_data$State[which(IGPC_Tot_data$BioD > 0 & IGPC_Tot_data$Tot_eq == 2 & IGPC_Tot_data$Pi_eq > Tr)] <- "Occupancy";
#### Changes in system stability in the community (initially vs. after 5000 years) ====
# 0. Initial-Final configurations
conf_name <- c(); initeq <- c("Eq", "Cycles"); fineq <- c("Null", "Eq", "Cycles")
for(i in initeq){conf_name <- append(conf_name, paste(i, fineq, sep='.'))};
conf_code <- c(0:5);
Conf_tab <- cbind.data.frame(conf_name, conf_code);

summary(as.factor(IGPC_Tot_data$Zone_i))
summary(as.factor(IGPC_Tot_data$Zone_f))
# A. Prepare the different combinations of stability changes
init_eq_code <- c(1, 2, 3);
init_eq <- c(rep("Eq", 2), "Cycles");
final_eq_code <- c(0:3, 6:8);
final_eq <- c("Null", "Eq", rep( c("Eq", "Cycles"), 2), "Eq");

# Data of observe combinations
I_code <- c(); for(i in init_eq_code){I_code <- append(I_code, rep(i, length(final_eq)))};
F_code <- rep(final_eq_code, length(init_eq_code));

IF_eq <- c(); for(i in init_eq){IF_eq <- append(IF_eq, paste(i, final_eq, sep='.'))}
IF_code <- rep(0, length(IF_eq));
IF <- cbind.data.frame(IF_eq, IF_code, I_code, F_code);

IF$IF_eq <- as.character(IF$IF_eq);
for(i in 1:length(IF_eq)){ IF$IF_code[i] <- Conf_tab$conf_code[which(Conf_tab$conf_name==IF_eq[i])]};
# B. Attribute the AC_Tot_data's outcomes to the corresponding configurations EQ code and name
IGPC_Tot_data$Zone_i <- as.numeric(IGPC_Tot_data$Zone_i);
IGPC_Tot_data$Zone_f <- as.numeric(IGPC_Tot_data$Zone_f);
IGPC_Tot_data$EQ_transit_InitFin <- as.character(IGPC_Tot_data$EQ_transit_InitFin)
IGPC_Tot_data$Zone_EQtransit_InitFin <- as.numeric(IGPC_Tot_data$Zone_EQtransit_InitFin);
for(i in 1:dim(IF)[1]){
  IGPC_Tot_data$EQ_transit_InitFin[which(IGPC_Tot_data$Zone_i == IF$I_code[i] & IGPC_Tot_data$Zone_f==IF$F_code[i])] <- IF$IF_eq[i];
  IGPC_Tot_data$Zone_EQtransit_InitFin[which(IGPC_Tot_data$Zone_i == IF$I_code[i] & IGPC_Tot_data$Zone_f==IF$F_code[i])] <- IF$IF_code[i]
}
#### Changes in system stability between resident and invaded communities (both after 5000 years) ====
Res_zonef <- rep(0, dim(IGPC_Tot_data)[1]);
IGPC_Tot_data <- cbind.data.frame(IGPC_Tot_data, Res_zonef);
# 0. Resident-Invaded configurations
conf_name <- c(); initeq <- c("Null", "Eq", "Cycles");
for(i in initeq){conf_name <- append(conf_name, paste(i, initeq, sep='.'))};
conf_code <- c(0:8);
Conf_tab <- cbind.data.frame(conf_name, conf_code);
# A. Merge Res$Zone_f to invaded data frame
for(i in Alpha){
  for(j in Beta){
    subinv <- IGPC_Tot_data[which(IGPC_Tot_data$Alpha==i & IGPC_Tot_data$Beta==j),];
    subres <- Res[which(Res$Alpha==i),];
    IGPC_Tot_data$Res_zonef[which(IGPC_Tot_data$Alpha == i & IGPC_Tot_data$Beta == j)] <- subres$Zone_f;
  }
}
# B. Prepare the different combinations of stability changes
res_eq_code <- c(0, 1, 2, 3);
res_eq <- c("Null", rep("Eq", 2), "Cycles");
# EQ & Code for invaded community
inv_eq_code <- c(0:3, 6:8);
inv_eq <- c("Null", "Eq", rep( c("Eq", "Cycles"), 2), "Eq");
# Data of observe combinations
Res_code <- c(); for(i in res_eq_code){Res_code <- append(Res_code, rep(i, length(inv_eq)))};
Inv_code <- rep(inv_eq_code, length(res_eq_code));

Inv_Res_eq <- c(); for(i in res_eq){Inv_Res_eq <- append(Inv_Res_eq, paste(i, inv_eq, sep='.'))}
Inv_Res_code <- rep(0, length(Inv_Res_eq));
Inv_Res <- cbind.data.frame(Inv_Res_eq, Inv_Res_code, Res_code, Inv_code);

for(i in 1:length(Inv_Res_eq)){ Inv_Res$Inv_Res_code[i] <- Conf_tab$conf_code[which(Conf_tab$conf_name==Inv_Res_eq[i])]};
# C. Attribute the AC_Tot_data's outcomes to the corresponding configurations EQ code and name
IGPC_Tot_data$EQ_transit_ResInv <- as.character(IGPC_Tot_data$EQ_transit_ResInv)
IGPC_Tot_data$Zone_f <- as.numeric(IGPC_Tot_data$Zone_f);
IGPC_Tot_data$Res_zonef <- as.numeric(IGPC_Tot_data$Res_zonef);
IGPC_Tot_data$Zone_EQtransit_ResInv <- as.numeric(IGPC_Tot_data$Zone_EQtransit_ResInv);

for(i in 1:dim(Inv_Res)[1]){
  IGPC_Tot_data$EQ_transit_ResInv[which(IGPC_Tot_data$Res_zonef == Inv_Res$Res_code[i] & IGPC_Tot_data$Zone_f==Inv_Res$Inv_code[i])] <- as.character(Inv_Res$Inv_Res_eq)[i];
  IGPC_Tot_data$Zone_EQtransit_ResInv[which(IGPC_Tot_data$Res_zonef == Inv_Res$Res_code[i] & IGPC_Tot_data$Zone_f==Inv_Res$Inv_code[i])] <- Inv_Res$Inv_Res_code[i]
}
#### Savig the dataset ====
write.table(IGPC_Tot_data, "./data/IGPC_Tot_data.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);