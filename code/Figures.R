########################################################################
### R Script                     Figures.R
###
###   Codes used to obtain the selected figures
########################################################################
######## Required packages ====
library(deSolve)
library(rootSolve)
library(ggplot2)
library(lattice)
library(gridExtra)

#library(devtools)
#devtools::install_github("ropenscilabs/ochRe")
library(RColorBrewer)
library(ochRe)

######## Required functions
source("./code/Functions.R"); # Load all required functions

######## Color selections used across the Figures ====
##### Color palette for Biodiversity change
col.bio <- c(brewer.pal(n=11, name='RdYlBu')[c(10,8)],  # -2, -1
             brewer.pal(n=9, name='Greys')[4],          # 0
             brewer.pal(n=11, name='RdYlBu')[c(5,3,1)]);# +1, +2, +3
tab_bio <-  data.frame(matrix(0, ncol=2, nrow=length(col.bio)));
colnames(tab_bio) <- c("Bio", "Colours");
tab_bio$Bio <- c(-2, -1, 0, 1, 2, 3);
tab_bio$Colours <- col.bio;
atpar_bio = seq(-2.5, 3.5, by=1); seq_bio = seq(-2, 3, by=1);

##### Color palette for Community composition
col.compo <- c(brewer.pal(n=9, name='Greens')[5],       # 0 species
               brewer.pal(n=9, name='Purples')[c(5,7)], # R eq, Cr-Rr eq
               brewer.pal(n=9, name='Oranges')[6],      # Cr-Rr cycles
               brewer.pal(n=11, name='PuOr')[c(9,3)],   # Cr-Ri
               brewer.pal(n=11, name='PuOr')[c(10,2)],  # Ci-Rr
               brewer.pal(n=11, name='PuOr')[c(11,1)] ) # 3sp

tab_compo <- data.frame(matrix(0, ncol=3, nrow=length(col.compo)));
colnames(tab_compo) <- c("Compo", "Colours", "Zone");
tab_compo$Compo <- c("0.sp", "R eq", "Cr.Rr eq", "Cr.Rr cycles", "Cr.Ri eq", "Cr.Ri cycles", "Ci.Rr eq", "Ci.Rr cycles", "3.sp eq", "3.sp cycles");
tab_compo$Colours <- col.compo;
tab_compo$Zone <- seq(0,9)
atpar_eq = seq(-0.5, 9.5, by=1);

#2spAC
tab_compo_2AC_i <- tab_compo[c(2:4),];
tab_compo_2AC_f <- tab_compo[c(1:4),];
#AC
tab_compo_AC <- tab_compo[c(1:6),];
tab_compo_ACb <- tab_compo[c(1:6, 9),];
#EC
tab_compo_EC <- tab_compo[c(1:4,7:8),];
tab_compo_ECb <- tab_compo[c(1:4,7:9),];
#TC
tab_compo_TC <- tab_compo[c(1:4, 9:10),];
#IGP
tab_compo_IGPB <- tab_compo; tab_compo_IGPC <- tab_compo;
tab_compo_IGPB[c(3:4), 1] <- c("Pr.Rr eq", "Pr.Rr cycles");
tab_compo_IGPC[c(5:6), 1] <- c("Pi.Rr eq", "Pi.Rr cycles");

##### Color palette for Invasion states
# using the ochRe palette; olsen qual
pal <- as.vector(unlist(ochre_palettes[8]));
col.state <- c(brewer.pal(n=9, name='Greys')[9],   # Vulnerability
               pal[5],   # Resistance
               pal[6],   # Substitution
               pal[2],   # Rescue
               pal[4],   # Occupancy
               pal[1]);  # Integration

tab_state <- data.frame(matrix(0, ncol=4, nrow=6));
colnames(tab_state) <- c("State", "Code", "Colours", "Symbols");
tab_state$State <- c("Vulnerability", "Resistance", "Substitution", "Rescue", "Occupancy", "Integration");
tab_state$Code <- seq(0,5, by=1);
tab_state$Colours <- col.state;
tab_state$Symbols <- c(17,16,16,15,15,15);

##### Color palette for Changes in stability regime
# N.N. O.N. E.N.:
colN <- c(brewer.pal(n=9, name='Greens')[c(5,7,9)]);
# N.O. O.O. E.O.:
colC <- c(brewer.pal(n=9, name='Oranges')[c(4,6,9)]);
# N.E. O.E. E.E.:
colE <- c(brewer.pal(n=9, name='Purples')[c(5,7,9)]);

# Tab for change in stability regime Resident.Invasion
tab_stab_ResInv<- data.frame(matrix(0, ncol=6, nrow=9));
colnames(tab_stab_ResInv) <- c("Stab", "State", "Code", "Colours", "Symbols", "Linetypes");

stab <- c(); vecstab <- c("Null", "Eq", "Cycles");
for(i in vecstab){stab <- append(stab, paste(i, vecstab, sep='.'))};
tab_stab_ResInv$Stab <- stab;
colRI <- NULL; for(i in c(1,3,2)){colRI <- append(colRI, c(colN[i], colE[i], colC[i]))}
tab_stab_ResInv$Colours <- colRI;
tab_stab_ResInv$Code <- seq(0,8);
tab_stab_ResInv$State <- c("N.N", "N.E", "N.O", "E.N", "E.E", "E.O", "O.N", "O.E", "O.O");
tab_stab_ResInv$Symbols <- c(16,15,15,17,16,17,17,15,16);
tab_stab_ResInv$Linetypes <- c(3,1,1,1,3,1,1,1,3);
atpar_stab_ResInv = seq(-0.5, 8.5, by=1);

# Tab for change in stability regime Initial.Final (t0 vs t=5000y)
tab_stab_IniFin <- tab_stab_ResInv[c(4:dim(tab_stab_ResInv)[1]),];
atpar_stab_InitFin = seq(-0.5, 5.5, by=1);

# Tab of Regime state before invasion 
tab_stab_Res <- tab_stab_ResInv[c(1,5,9),c(1:2,4:5)];
tab_stab_Res$Stab <- c("Null", "Eq", "Cycles");
tab_stab_Res$State <- c("N", "E", "O");

##### Color palette To illustrate the outcomes of R* and P* rules
col.Rrule <- c(brewer.pal(n=11, name='PuOr')[3],
               brewer.pal(n=9, name='Greys')[6],
               brewer.pal(n=11, name='PuOr')[9]);

col.Prule <- c(brewer.pal(n=11, name='PuOr')[9],
               brewer.pal(n=9, name='Greys')[6],
               brewer.pal(n=11, name='PuOr')[3]);
##### Palette for the Populations Biomass extremes (Figs S6-7)
pal <- as.vector(unlist(ochre_palettes[8]));
col.bmr <- c(pal[2], # bmr=1
             pal[4], # bmr=4
             pal[1], # bmr=25
             pal[3]);# bmr=100

##### Palette for Environmental regions
col.reg <- c(brewer.pal(n=3, name='Paired')[1],
             brewer.pal(n=8, name='Set2')[6],
             brewer.pal(n=8, name='Dark2')[c(5,2)]);
########################################################################
####################### MAIN FIGURES =====
######## Fig. 1 => see Fig. S1a to reproduce Fig.1b
######## Fig. 2: Ggplots of Invasion outcomes and Regime states along environmental gradients across trophic modules and fixed body mass ratio (Gamma=100, Alpha= Beta= 10) ====
Fw <- c("AC", "EC", "TC", "IGPB", "IGPC");

# AC ====
name <- paste(paste(Fw[1],"Tot", "data", sep='_'),"txt", sep='.');
F1 <- read.table( paste(".","data",name,sep='/'), h=T,  sep = "\t");

sub = F1[which(F1$Gamma == 100),];
code_state <- rep(0, dim(sub)[1]);
sub <- cbind.data.frame(sub, code_state);
for(j in 1:dim(tab_state)[1]){
  if(length(sub$code_state[which(sub$State == tab_state$State[j])]) != 0){
    sub$code_state[which(sub$State == tab_state$State[j]) ] <- tab_state$Code[j];
  }else{next}
}

p1 <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(code_state)))+theme_bw()+
  geom_tile()+
  labs(x="", y=expression("Temperature ("*degree*C*')'))+
  scale_fill_manual(values = tab_state$Colour, breaks=tab_state$Code,
                    labels=tab_state$State, name=NULL, guide="none")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))

p2 <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(Zone_EQtransit_ResInv)))+theme_bw()+
  geom_tile()+
  labs(x=expression("Nutrient levels ("*g.L^-1*')'), y=expression("Temperature ("*degree*C*')'))+
  scale_fill_manual(values = tab_stab_ResInv$Colour, breaks=tab_stab_ResInv$Code,
                    labels=tab_stab_ResInv$Stab, name=NULL, guide="none")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))

# EC ====
name <- paste(paste(Fw[2],"Tot", "data", sep='_'),"txt", sep='.');
F1 <- read.table( paste(".","data",name,sep='/'), h=T,  sep = "\t");

sub = F1[which(F1$Gamma == 100),];
code_state <- rep(0, dim(sub)[1]);
sub <- cbind.data.frame(sub, code_state);
for(j in 1:dim(tab_state)[1]){
  if(length(sub$code_state[which(sub$State == tab_state$State[j])]) != 0){
    sub$code_state[which(sub$State == tab_state$State[j]) ] <- tab_state$Code[j];
  }else{next}
}

p3 <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(code_state)))+theme_bw()+
  geom_tile()+
  labs(x="", y="")+
  scale_fill_manual(values = tab_state$Colour, breaks=tab_state$Code,
                    labels=tab_state$State, name=NULL, guide="none")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))

p4 <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(Zone_EQtransit_ResInv)))+theme_bw()+
  geom_tile()+
  labs(x=expression("Nutrient levels ("*g.L^-1*')'), y="")+
  scale_fill_manual(values = tab_stab_ResInv$Colour, breaks=tab_stab_ResInv$Code,
                    labels=tab_stab_ResInv$Stab, name=NULL, guide="none")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))

# TC ====
name <- paste(paste(Fw[3],"Tot", "data", sep='_'),"txt", sep='.');
F1 <- read.table( paste(".","data",name,sep='/'), h=T,  sep = "\t");

sub = F1[which(F1$Gamma == 100),];
code_state <- rep(0, dim(sub)[1]);
sub <- cbind.data.frame(sub, code_state);
for(j in 1:dim(tab_state)[1]){
  if(length(sub$code_state[which(sub$State == tab_state$State[j])]) != 0){
    sub$code_state[which(sub$State == tab_state$State[j]) ] <- tab_state$Code[j];
  }else{next}
}

p5 <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(code_state)))+theme_bw()+
  geom_tile()+
  labs(x="", y="")+
  scale_fill_manual(values = tab_state$Colour, breaks=tab_state$Code,
                    labels=tab_state$State, name=NULL, guide="none")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))

p6 <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(Zone_EQtransit_ResInv)))+theme_bw()+
  geom_tile()+
  labs(x=expression("Nutrient levels ("*g.L^-1*')'), y="")+
  scale_fill_manual(values = tab_stab_ResInv$Colour, breaks=tab_stab_ResInv$Code,
                    labels=tab_stab_ResInv$Stab, name=NULL, guide="none")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))

# IGPB ====
name <- paste(paste(Fw[4],"Tot", "data", sep='_'),"txt", sep='.');
F1 <- read.table( paste(".","data",name,sep='/'), h=T,  sep = "\t");

sub = F1[which(F1$Gamma == 100),];
code_state <- rep(0, dim(sub)[1]);
sub <- cbind.data.frame(sub, code_state);
for(j in 1:dim(tab_state)[1]){
  if(length(sub$code_state[which(sub$State == tab_state$State[j])]) != 0){
    sub$code_state[which(sub$State == tab_state$State[j]) ] <- tab_state$Code[j];
  }else{next}
}

p7 <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(code_state)))+theme_bw()+
  geom_tile()+
  labs(x="", y="")+
  scale_fill_manual(values = tab_state$Colour, breaks=tab_state$Code,
                    labels=tab_state$State, name=NULL, guide="none")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))

p8 <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(Zone_EQtransit_ResInv)))+theme_bw()+
  geom_tile()+
  labs(x=expression("Nutrient levels ("*g.L^-1*')'), y="")+
  scale_fill_manual(values = tab_stab_ResInv$Colour, breaks=tab_stab_ResInv$Code,
                    labels=tab_stab_ResInv$Stab, name=NULL, guide="none")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))

# IGPC ====
name <- paste(paste(Fw[5],"Tot", "data", sep='_'),"txt", sep='.');
F1 <- read.table( paste(".","data",name,sep='/'), h=T,  sep = "\t");

sub = F1[which(F1$Gamma == 100),];
code_state <- rep(0, dim(sub)[1]);
sub <- cbind.data.frame(sub, code_state);
for(j in 1:dim(tab_state)[1]){
  if(length(sub$code_state[which(sub$State == tab_state$State[j])]) != 0){
    sub$code_state[which(sub$State == tab_state$State[j]) ] <- tab_state$Code[j];
  }else{next}
}

p9 <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(code_state)))+theme_bw()+
  geom_tile()+
  labs(x="", y="")+
  scale_fill_manual(values = tab_state$Colour, breaks=tab_state$Code,
                    labels=tab_state$State, name=NULL, guide="none")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))

p10 <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(Zone_EQtransit_ResInv)))+theme_bw()+
  geom_tile()+
  labs(x=expression("Nutrient levels ("*g.L^-1*')'), y="")+
  scale_fill_manual(values = tab_stab_ResInv$Colour, breaks=tab_stab_ResInv$Code,
                    labels=tab_stab_ResInv$Stab, name=NULL, guide="none")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))

# Arrange & Save ====
pdf(file="Fig_1.pdf"); # here without dimension adjustment for the figure
grid.arrange(p1, p3, p5, p7, p9, p2, p4, p6, p8, p10, ncol=5, nrow=2, newpage = T)
dev.off()
######## Fig. 3: Averaged proportions (Invasion success, diversity change and stability change) observed within 4 environmental regions ====
# along gradients of body mass ratio between competing species (AC, EC, IGP) and trophic levels (TC and IPGP)
# data ====
IS_Reg_AC <- read.table("./data/BMR_AC_IS_RegionS.txt", h=T,  sep = "\t");
RS_Reg_AC <- read.table("./data/BMR_AC_RS-ResInv_Regions.txt", h=T,  sep = "\t");
IS_Reg_AC <- IS_Reg_AC[-which(IS_Reg_AC$Alpha==1),]; RS_Reg_AC <- RS_Reg_AC[-which(RS_Reg_AC$Alpha==1),];
colnames(IS_Reg_AC)[13] <- colnames(RS_Reg_AC)[11] <- "BMR";

IS_Reg_EC <- read.table("./data/BMR_EC_IS_Regions.txt", h=T,  sep = "\t");
RS_EC <- read.table("./data/BMR_EC_RS-ResInv.txt", h=T,  sep = "\t");
RS_Reg_EC <- read.table("./data/BMR_EC_RS-ResInv_Regions.txt", h=T,  sep = "\t");
IS_Reg_EC <- IS_Reg_EC[-which(IS_Reg_EC$Beta==1),]; RS_Reg_EC <- RS_Reg_EC[-which(RS_Reg_EC$Beta==1),]; 

colnames(IS_Reg_EC)[13] <- colnames(RS_Reg_EC)[14] <- "BMR";
IS_Reg_EC[,13] <- (1/IS_Reg_EC$BMR); # Correct bmr as 1/Beta
RS_Reg_EC[,14] <- (1/RS_Reg_EC$BMR);

# IGP
IS_Reg_IGP <- read.table("./data/BMR_IGP_IS_Beta_Regions.txt", h=T,  sep = "\t");
RS_Reg_IGP <- read.table("./data/BMR_IGP_RS-ResInv_Beta_Regions.txt", h=T,  sep = "\t");
colnames(IS_Reg_IGP)[13] <- colnames(RS_Reg_IGP)[12] <- "BMR";

#TC
IS_Reg_TC <- read.table("./data/BMR_TC_IS_Delta_RegionS.txt", h=T,  sep = "\t");
RS_Reg_TC <- read.table("./data/BMR_TC_RS-ResInv_Delta_Regions.txt", h=T,  sep = "\t");
colnames(IS_Reg_TC)[11] <- colnames(RS_Reg_TC)[14] <- "BMR";

#IGPC
IS_Reg_IGPC <- read.table("./data/BMR_IGPC_IS_Delta_RegionS.txt", h=T,  sep = "\t");
RS_Reg_IGPC <- read.table("./data/BMR_IGPC_RS-ResInv_Delta_Regions.txt", h=T,  sep = "\t");
colnames(IS_Reg_IGPC)[12] <- colnames(RS_Reg_IGPC)[12] <- "BMR";

# figures ====
Fw <- c("AC", "EC", "IGP", "TC", "IGPC");
tab_name <- c("IS_Reg_","IS_Reg_", "RS_Reg_");
subdat <- c("IS", "IS", "RS");
subvec <- c("Bio","DeltaS");
col.reg <- c("lightblue", "yellow","green","orange");

pdf(file="Fig_3.pdf", width=15, height=9);

par(mfrow = c(3, length(Fw)));
for(r in 1:3){
  for(i in 1:length(Fw) ){
    par(mfg=c(r,i), mar = c(3, 4, 0.5, 0.5));
    plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(-1.1,1.1), ylim=c(-5,105));
    #if(r==1){axis(3, at=log10(1), labels=Fw[i], cex.axis=2.5) }else{ axis(3, at=log10(1), labels=F)}  ;
    if(r==3){axis(1, at=log10(c(0.1,1,10)), labels=c(expression(10^{-1}), 1, expression(10^{1})), cex.axis=2) };
    if(i==1){axis(2, at=c(0,25,50,75,100), label=c(0,"",50,"",100), cex.axis=2, las=2)};
    abline(v=log10(1), lty=3);
    
    if(r == 1){
      for(j in 1:4){
        YYY <- get(paste0(tab_name[r],Fw[i]))[which(get(paste0(tab_name[r],Fw[i]))$Regions==LETTERS[j]),]$Inv_succ;
        lines(get(paste0(tab_name[r],Fw[i]))[which(get(paste0(tab_name[r],Fw[i]))$Regions==LETTERS[j]),]$Inv_succ~log10(get(paste0(tab_name[r],Fw[i]))[which(get(paste0(tab_name[r],Fw[i]))$Regions==LETTERS[r]),]$BMR), type="l", lwd=2, col=col.reg[j]);
      }
    }
    if(r == 2){
      for(j in 1:4){
        lines(get(paste0(subdat[r],'_Reg_',Fw[i]))[which(get(paste0(subdat[r],'_Reg_',Fw[i]))$Regions==LETTERS[j]),][paste0(subvec[1],'_gain')][,1]~log10(get(paste0(subdat[r],'_Reg_',Fw[i]))[which(get(paste0(subdat[r],'_Reg_',Fw[i]))$Regions==LETTERS[j]),]$BMR), type="l", lwd=2, col=col.reg[j], lty=1); #gain
        lines(get(paste0(subdat[r],'_Reg_',Fw[i]))[which(get(paste0(subdat[r],'_Reg_',Fw[i]))$Regions==LETTERS[j]),][paste0(subvec[1],'_loss')][,1]~log10(get(paste0(subdat[r],'_Reg_',Fw[i]))[which(get(paste0(subdat[r],'_Reg_',Fw[i]))$Regions==LETTERS[j]),]$BMR), type="l", lwd=2, col=col.reg[j], lty=2); #loss 
      }
    }
    if(r == 3){
      for(j in 1:4){
        if(i==5){
          lines(get(paste0(subdat[r],'_Reg_',Fw[i]))[which(get(paste0(subdat[r],'_Reg_',Fw[i]))$Regions==LETTERS[j]),][paste0(subvec[2],'_gain')][,1]~log10(sort(get(paste0(subdat[r],'_Reg_',Fw[i]))[which(get(paste0(subdat[r],'_Reg_',Fw[i]))$Regions==LETTERS[j]),]$BMR,decreasing=F)), type="l", lwd=2, col=col.reg[j], lty=1); #gain
          lines(get(paste0(subdat[r],'_Reg_',Fw[i]))[which(get(paste0(subdat[r],'_Reg_',Fw[i]))$Regions==LETTERS[j]),][paste0(subvec[2],'_loss')][,1]~log10(sort(get(paste0(subdat[r],'_Reg_',Fw[i]))[which(get(paste0(subdat[r],'_Reg_',Fw[i]))$Regions==LETTERS[j]),]$BMR,decreasing=F)), type="l", lwd=2, col=col.reg[j], lty=2); #loss
        }else{
          lines(get(paste0(subdat[r],'_Reg_',Fw[i]))[which(get(paste0(subdat[r],'_Reg_',Fw[i]))$Regions==LETTERS[j]),][paste0(subvec[2],'_gain')][,1]~log10(get(paste0(subdat[r],'_Reg_',Fw[i]))[which(get(paste0(subdat[r],'_Reg_',Fw[i]))$Regions==LETTERS[j]),]$BMR), type="l", lwd=2, col=col.reg[j], lty=1); #gain
          lines(get(paste0(subdat[r],'_Reg_',Fw[i]))[which(get(paste0(subdat[r],'_Reg_',Fw[i]))$Regions==LETTERS[j]),][paste0(subvec[2],'_loss')][,1]~log10(get(paste0(subdat[r],'_Reg_',Fw[i]))[which(get(paste0(subdat[r],'_Reg_',Fw[i]))$Regions==LETTERS[j]),]$BMR), type="l", lwd=2, col=col.reg[j], lty=2); #loss
        }
      }
    }
  }
}
dev.off()
######## Fig. 4: Barplots of Regime states ~ Invading mechanisms ====
# data ====
Inv.AC <- read.table("./data/AC_Tot_data.txt", h=T,  sep = "\t");  # AC Complete Data
Inv.EC <- read.table("./data/EC_Tot_data.txt", h=T,  sep = "\t");  # EC Complete Data
Inv.TC <- read.table("./data/TC_Tot_data.txt", h=T,  sep = "\t");  # TC Complete Data
Inv.IGPB <- read.table("./data/IGPB_Tot_data.txt", h=T,  sep = "\t");  # IGPB Complete Data
Inv.IGPC <- read.table("./data/IGPC_Tot_data.txt", h=T,  sep = "\t");  # IGPC Complete Data

AC2 <-Inv.AC[which(Inv.AC$Alpha!=1),]; AC <- AC2;
EC2 <- Inv.EC[which(Inv.EC$Beta!=1),]; EC <- EC2;
TC <- Inv.TC
IGPB <- Inv.IGPB
IGPC <- Inv.IGPC

data <- rbind.data.frame(cbind(AC2$State, AC2$EQ_transit_ResInv),
                         cbind(EC2$State, EC2$EQ_transit_ResInv), 
                         cbind(TC$State, TC$EQ_transit_ResInv),
                         cbind(IGPB$State, IGPB$EQ_transit_ResInv),
                         cbind(IGPC$State, IGPC$EQ_transit_ResInv))
colnames(data) <- c("State", "Regimes");
State <- tab_state$State;
N_Obs <- dim(data)[1];

# Proportions of regime states for each invasion outcomes ====
Tab_obs_State_Regime <- as.data.frame(matrix(0, ncol=4));
colnames(Tab_obs_State_Regime) <- c("Mec", "Lev", "Obs", "Prob");

for(i in 1:length(State)){
  Lev <- levels(factor(data$Regimes[data$State==State[i]]));
  Obs <- as.numeric(summary(as.factor(data$Regimes[data$State==State[i]])));
  Mec <- rep(State[i], length(Lev));
  Prob <- Proportion(Obs,sum(Obs))
  sub <- cbind.data.frame(Mec, Lev, Obs, Prob)
  Tab_obs_State_Regime <- rbind(Tab_obs_State_Regime, sub )
}

Tab_obs_State_Regime <- Tab_obs_State_Regime[-1,];
Tab_obs_State_Regime <- cbind(Tab_obs_State_Regime, rep(NA, dim(Tab_obs_State_Regime)[1]), rep(NA, dim(Tab_obs_State_Regime)[1]), rep(NA, dim(Tab_obs_State_Regime)[1]), rep(NA, dim(Tab_obs_State_Regime)[1]))
colnames(Tab_obs_State_Regime) <- c("Mechanisms", "Regimes", "N_obs", "P_obs", "Colours", "Reg", "Reg.Lev", "Div.Lev")

for(i in 1:dim(tab_stab_ResInv)[1]){
  Tab_obs_State_Regime$Colours[Tab_obs_State_Regime$Regimes==tab_stab_ResInv$Stab[i]] <- tab_stab_ResInv$Colours[i];
  Tab_obs_State_Regime$Reg[Tab_obs_State_Regime$Regimes==tab_stab_ResInv$Stab[i]] <- tab_stab_ResInv$State[i];
}

Destabilizing.S <- c("O.N", "E.O", "E.N");
Neutral.S <- c("N.N", "O.O", "E.E");
Stabilizing.S <- c("N.E", "N.O", "O.E");
Div.loss <- c("Vulnerability");
Div.neut <- c("Resistance", "Substitution");
Div.gain <- c("Integration", "Rescue", "Occupancy");

Tab_obs_State_Regime$Div.Lev[Tab_obs_State_Regime$Mechanisms=="Vulnerability"] <- LETTERS[1];
for(i in 1:3){
  if(i<3){ Tab_obs_State_Regime$Div.Lev[Tab_obs_State_Regime$Mechanisms==Div.neut[i]] <- LETTERS[2];}
  Tab_obs_State_Regime$Reg.Lev[Tab_obs_State_Regime$Reg==Destabilizing.S[i]] <- LETTERS[1];
  Tab_obs_State_Regime$Reg.Lev[Tab_obs_State_Regime$Reg==Neutral.S[i]] <- LETTERS[2];
  Tab_obs_State_Regime$Reg.Lev[Tab_obs_State_Regime$Reg==Stabilizing.S[i]] <- LETTERS[3];
  Tab_obs_State_Regime$Div.Lev[Tab_obs_State_Regime$Mechanisms==Div.gain[i]] <- LETTERS[3];
}

# Removing Resistance from the outcome (for only for Fig 4):
Tab_obs_State_Regime <- Tab_obs_State_Regime[which(Tab_obs_State_Regime$Mechanisms!="Resistance"),];

# Delta S regarding to Invasion Mechanisms
Reg.Lev <- c(); Mechanisms <- c();
for(j in 1:length(State)){
  Lev <- levels(factor(Tab_obs_State_Regime$Reg.Lev[Tab_obs_State_Regime$Mechanisms==State[j]]));
  Mech <- rep(State[j], length(Lev));
  Reg.Lev <- append(Reg.Lev, Lev);
  Mechanisms <- append(Mechanisms, Mech);
}
tab_one <- cbind.data.frame(Mechanisms, Reg.Lev, rep(NA, length(Mechanisms)), rep(NA, length(Mechanisms)), rep(NA, length(Mechanisms)))
colnames(tab_one) <- c("Mechanisms", "Reg.Lev", "N_obs", "P_obs", "Colours")

for(i in 1:dim(tab_one)[1]){
  Totobs <- sum(Tab_obs_State_Regime$N_obs[Tab_obs_State_Regime$Mechanisms==tab_one$Mechanisms[i] ]);
  tab_one$N_obs[i] <- sum(Tab_obs_State_Regime$N_obs[Tab_obs_State_Regime$Mechanisms==tab_one$Mechanisms[i] & Tab_obs_State_Regime$Reg.Lev==tab_one$Reg.Lev[i]])
  tab_one$P_obs[i] <- Proportion(tab_one$N_obs[i],Totobs)
}
tab_one$Colours[which(tab_one$Reg.Lev=="A")] <- brewer.pal(n=9, name='Greys')[9]
tab_one$Colours[which(tab_one$Reg.Lev=="B")] <- brewer.pal(n=9, name='Greys')[1]
tab_one$Colours[which(tab_one$Reg.Lev=="C")] <- brewer.pal(n=9, name='Greys')[4]

## Delta S regarding to Delta D
Div.Lev <- c(); Reg.Lev <- c();
for(i in 1:3){
  Lev <- levels(factor(Tab_obs_State_Regime$Reg.Lev[Tab_obs_State_Regime$Div.Lev==LETTERS[i]]));
  Reg.Lev <- append(Reg.Lev, Lev);
  Div.Lev <- append(Div.Lev, rep(LETTERS[i], length(Lev) ));
}
tab_two <- cbind.data.frame(Div.Lev, Reg.Lev, rep(NA, length(Div.Lev)), rep(NA, length(Div.Lev)), rep(NA, length(Div.Lev)))
colnames(tab_two) <- c("Div.Lev", "Reg.Lev", "N_obs", "P_obs", "Colours")

for(i in 1:dim(tab_two)[1]){
  Totobs <- sum(Tab_obs_State_Regime$N_obs[Tab_obs_State_Regime$Div.Lev==tab_two$Div.Lev[i] ]);
  tab_two$N_obs[i] <- sum(Tab_obs_State_Regime$N_obs[Tab_obs_State_Regime$Div.Lev==tab_two$Div.Lev[i] & Tab_obs_State_Regime$Reg.Lev==tab_two$Reg.Lev[i]])
  tab_two$P_obs[i] <- Proportion(tab_two$N_obs[i],Totobs)
}
tab_two$Colours[which(tab_two$Reg.Lev=="A")] <- brewer.pal(n=9, name='Greys')[9]
tab_two$Colours[which(tab_two$Reg.Lev=="B")] <- brewer.pal(n=9, name='Greys')[1]
tab_two$Colours[which(tab_two$Reg.Lev=="C")] <- brewer.pal(n=9, name='Greys')[4]

# Barplots ====
q <-ggplot(data=tab_one, aes(x=Mechanisms, y=P_obs, color=Reg.Lev))+
  geom_bar(stat="identity", color="black", fill=tab_one$Colours)+
  labs(x="Invasion mechanisms", y= expression("Stability change"~Delta*S~"(%)"))+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15), legend.position="none");

r <-ggplot(data=tab_two, aes(x=Div.Lev, y=P_obs, color=Reg.Lev))+
  geom_bar(stat="identity", color="black", fill=tab_two$Colours)+
  scale_x_discrete(labels=c("A" = expression(Delta*D~"< 0"), "B" = expression(Delta*D~"= 0"), "C" = expression(Delta*D~"> 0") ))+
  labs(x="Biodiversity change", y= expression("Stability change"~Delta*S~"(%)"))+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15), legend.position="none");

# Arrange & Save ====
pdf(file="Fig_4.pdf", width=15, height=9);
grid.arrange(q, r, ncol=2, nrow=1);
dev.off()

####################### SUPPLEMENTARY FIGURES =====
######## Fig. S1: Structure of resident communities along gradients of abiotic condition and size structure ====
Res <- read.table("./data/CR_Tot_Gamma.txt", h=T, sep="\t", dec='.');
Res.bmr <- read.table("./data/BMR_2spAC_RS-F.txt", h=T, sep="\t", dec='.');
sub = Res[which(Res$Gamma == 100),];

p1 <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(Zone_f) ))+theme_bw()+
  geom_tile()+
  labs(x=expression("Nutrient levels ("*g.L^-1*')'), y=expression("Temperature ("*degree*C*')'))+
  scale_fill_manual(values = tab_compo_2AC_f$Colours, breaks=levels(as.factor(sub$Zone_f)),
                    labels=c("N", "E (R)", "E (C-R)", "C (C-R)"), name=NULL, guide="none")+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))

p2 <- ggplot(Res.bmr)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=Null), colour=tab_stab_Res$Colours[1], size=1.2)+
  geom_line(aes(x=Gamma, y=Eq), colour=tab_stab_Res$Colours[2], size=1.2, linetype=1)+
  geom_line(aes(x=Gamma, y=R.eq), colour=tab_compo_2AC_f$Colours[2], size=1.2, linetype=2)+
  geom_line(aes(x=Gamma, y=Cr.Rr.eq), colour=tab_compo_2AC_f$Colours[3], size=1.2, linetype=3)+
  geom_line(aes(x=Gamma, y=Cycle), colour=tab_stab_Res$Colours[3], size=1.2, linetype=1)+
  geom_point(aes(x=Gamma, y=Null), colour=tab_stab_Res$Colours[1], shape=tab_stab_Res$Symbols[1] ,size=3.2)+
  geom_point(aes(x=Gamma, y=Eq), colour=tab_stab_Res$Colours[2], shape=tab_stab_Res$Symbols[2] ,size=3.2)+
  geom_point(aes(x=Gamma, y=Cr.Rr.cycles), colour=tab_stab_Res$Colours[3], shape=tab_stab_Res$Symbols[3] ,size=3.2)+
  labs(x=expression(C['RES']*':'*R['RES']*" mass ratio"), y="Regime states of resident system (%)")+
  scale_x_log10(breaks = c(1,5,10,25,50,100), labels=c(1,5,10,25,50,100))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))

# Arrange & Save
pdf(file="Fig_S1.pdf");
grid.arrange(p1, p2, ncol=2, nrow=1, newpage = T)
dev.off()
######## Fig. S2-S3: Sensitivity analysis (Temperature and mass dependencies of biological rates and Rules) ====
# Parameter settings according to temperature- and mass- dependencies ====
#Fixe parameters used in all functions
Boltz = 8.617*10^(-5);  # Boltzmann constant (eV/K)
T0    = 293.15;         # normalization temperature (K)
T0K   = 273.15;         # 0 degres in Kalvin
e     = 0.85;           # Conversion efficiency constant of consumed biomass into biomass gain (see Yodzis and Innes 1992).

#### parameter values for temperature- and size-dependencies of prey growth rate #
Ir    = -15.68;         # rate specific constant
Sr    = -0.25;          # rate specific scaling coefficient
Ear   = -0.84;          # activation energy (eV)

#### parameter values for temperature- and size-dependencies of carrying capacity #
Sk    = 0.28;           # rate specific scaling coefficient
Eak   = 0.71;           # activation energy (eV)

#### parameter values for temperature- and size-dependencies of metabolic rate #
Ix    = -16.54;         # rate specific constant 
Sx    = -0.31;          # rate specific scaling coefficient
Eax   = -0.69;          # activation energy (eV)

#### parameter values for temperature- and size-dependencies of handling time #
#parameter values for the handling time basal relationship
Iy      = 9.66;         # rate specific constant
Sy_pred = 0.47;         # rate specific scaling coefficient
Sy_prey = -0.45;        # rate specific scaling coefficient
Eay     = 0.26;         # activation energy (eV)
#quadratic relationship between handling time and body mass
IThm   = 1.92;          # intercept
S1_Thm = -0.48;         # linear slope term 1
S2_Thm = 0.0256;        # quadratic slope term 2
#quadratic relationship between handling time and temperature
IThT   = 0.5;           # intercept
S1_ThT = -0.055;        # linear slope term 1
S2_ThT = 0.0013;        # quadratic slope term 2

#### parameter values for temperature- and size-dependencies of attack rate #
#parameter values for the basal relationship
I_B0     = -13.1;       # rate specific constant (g m^-2)
SB0_pred = -0.8;        # rate specific scaling coefficient
SB0_prey = 0.25;        # rate specific scaling coefficient
EaB0     = -0.38;       # activation energy (eV)
#quadratic relationship between search rate and body mass
Iam    = -1.81;         # intercept
S1_am  = 0.39;          # linear slope term 1
S2_am  = -0.017;        # quadratic slope term 2


# Functions for Temperature- and mass dependent rates ====
#### temperature- and mass-dependent growth rate of basal species #
r<-function(Temp,Mx){
  r= exp(Ir)*Mx^Sr*exp(Ear*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
}

#### temperature- and mass-dependent carrying capacity of basal species #
K<-function(Temp,Mx,Car){
  K = Car*Mx^Sk*exp(Eak*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
}

#### temperature- and mass-dependent metabolic rate of consumer and top predator species #
x2<-function(Temp,My){
  x2= exp(Ix)*My^Sx*exp(Eax*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
}

#### temperature- and mass-dependent handling time of consumer and top predator species #
th<-function(Temp,My,Mx){
  # handling time 
  Th_basalY = exp(Iy) * My^Sy_pred * Mx^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # size-dependent component of handling time (quadratic correction, eq 2.13 in Binzer et al. 2012)
  Th_massY = exp(IThm + S1_Thm * log(My/Mx) + S2_Thm * (log(My/Mx))^2);
  # temperature-dependent component of handling time (quadratic correction, eq 2.14 in Binzer et al. 2012)
  Th_TempY = exp(IThT + S1_ThT * (Temp) + S2_ThT * (Temp)^2);
  # combine functions to get the final resultat
  th = Th_basalY*Th_massY*Th_TempY;
}

#### temperature- and mass-dependent attack rate of consumer and top predator species #
a<-function(Temp,My,Mx){
  # size-dependent half-saturation constant (baseline allometry, eq 2.10 in Binzer et al. 2012)
  a_baseY = exp(I_B0) * My^(SB0_pred) * Mx^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # size-dependent attack rate (quadratic correction, eq 2.12 in Binzer et al. 2012)
  a_massY = exp(Iam + S1_am * log(My/Mx) + S2_am * (log(My/Mx))^2);
  ##search rate for the intermediary predator preying on the prey
  a = a_baseY*a_massY;
}
# Sentivity analysis ====
bmr_AC <- c(0.1, 0.2, 0.5, 1, 2, 5, 10); #body mass ratio between invading:resident basal species
bmr_EC <- c(0.1, 0.2, 0.5, 1, 2, 5, 10); #body mass ratio between invading:resident consumer species
Mass_R_res <- 0.001; #body mass of resident basal resource (g)

Car = 5; # Fixed Nutrient levels
Temperature <- c(1,10,20,30, 39); # Fixed temperature
cc<- c("darkblue", "lightblue", "yellow", "orange", "darkred");

#The following functions can be found in 'Functions.R' script to calculate the P.Rule.Appcomp function
# P* rule
Mass_R_inv <- Mass_R_res*bmr_AC;  # body mass of invading basal species (lowest to largest resource used)
Mass_C_species <- 0.01;           # body mass of resident consumer

r_1 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
K_1 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
X_1 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
H_1 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
A_1 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
R_eq1 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
C_eq1 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
P_ratio <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );

for(i in 1:5){
  r_1[i,] <- r(Mx= Mass_R_inv, Temp = Temperature[i]);
  K_1[i,] <- K(Mx= Mass_R_inv, Temp = Temperature[i], Car = 5);
  X_1[i,] <- x2(Temp = Temperature[i], My=Mass_C_species);
  H_1[i,] <- th(Temp = Temperature[i], Mx= Mass_R_inv, My=Mass_C_species);
  A_1[i,] <- a(Temp = Temperature[i], Mx= Mass_R_inv, My=Mass_C_species);
  R_eq1[i,] <- as.numeric(X_1[i,])/(as.numeric(A_1[i,])*(e - as.numeric(X_1[i,])*as.numeric(H_1[i,]) ));
  C_eq1[i,] <- (as.numeric(r_1[i,])/(as.numeric(K_1[i,])*as.numeric(A_1[i,])))*(1+ (as.numeric(A_1[i,])*as.numeric(H_1[i,])*as.numeric(R_eq1[i,])))*(as.numeric(K_1[i,]) - as.numeric(R_eq1[i,]));
  #C_eq1[i, which(C_eq1[i, ]< 0)] <- NA;
  if(C_eq1[i,4] >0){ P_ratio[i,] <- (as.numeric(C_eq1[i,])/as.numeric(C_eq1[i,4])) }else{ P_ratio[i,] <-NA} 
};

# for R* rule
Mass_C_res <- 0.01;               # body mass of resident consumer (fixed for C:R bmr = 10)
Mass_C_inv <- Mass_C_res*bmr_EC;  # body mass of invading consumer species

r_2 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
K_2 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
X_2 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
H_2 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
A_2 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
R_eq2 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
C_eq2 <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );
R_ratio <- data.frame(matrix(NA,nrow= length(Temperature), ncol=length(Mass_R_inv)) );

for(i in 1:5){
  r_2[i,] <- r(Mx= Mass_R_res, Temp = Temperature[i]);
  K_2[i,] <- K(Mx= Mass_R_res, Temp = Temperature[i], Car = 5);
  X_2[i,] <- x2(Temp = Temperature[i], My=Mass_C_inv);
  H_2[i,] <- th(Temp = Temperature[i], Mx= Mass_R_res, My=Mass_C_inv);
  A_2[i,] <- a(Temp = Temperature[i], Mx= Mass_R_res, My=Mass_C_inv);
  R_eq2[i,] <- as.numeric(X_2[i,])/(as.numeric(A_2[i,])*(e - as.numeric(X_2[i,])*as.numeric(H_2[i,]) ));
  C_eq2[i,] <- (as.numeric(r_2[i,])/(as.numeric(K_2[i,])*as.numeric(A_2[i,])))*(1+ (as.numeric(A_2[i,])*as.numeric(H_2[i,])*as.numeric(R_eq2[i,])))*(as.numeric(K_2[i,]) - as.numeric(R_eq2[i,]));
  R_ratio[i,] <- R_eq2[i,]/R_eq2[i,4];
};
## Fig. S2: Biological rate ~ species mass (basal and consumer species) for varying temperatures ====
par(mfrow = c(2, 3));
# Plots along gradient of species mass of basal resource
par(mfg=c(1,1), mar = c(1, 5, 2, 1));
# resource intrinsic growth
plot(log10(as.numeric(r_1[1,]))~log10(Mass_R_inv), ylim=c(-7.4,-4.8),
     type='l', lwd=2, col=cc[1],  xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n");
axis(2, at=log10(c(1E-7,1E-6,1E-5)), labels=c(expression(10^{-7}), expression(10^{-6}), expression(10^{-5})), las=2, cex.axis=1.5);
axis(1, at=log10(Mass_R_inv), labels=F);
mtext("Intrinsic growth rate", 2, line=3.5, cex=1.5)
#mtext("Body mass of basal resource (g)", 1, line=2.2, cex=1.2);
abline(v=log10(1E-3), lty=3);
for(i in 2:5){lines(log10(as.numeric(r_1[i,]))~log10(Mass_R_inv), type='l', lwd=2, col=cc[i])}

par(mfg=c(2,1), mar = c(4, 5, 0.5, 1));
# resource carrying capacity
plot(log10(as.numeric(K_1[1,]))~log10(Mass_R_inv), ylim=c(-1.5,1),
     type='l', lwd=2, col=cc[1],  xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n");
axis(2, at=log10(c(1E-2,1E-1,1E0, 10)), labels=c(expression(10^{-2}), expression(10^{-1}), expression(10^{0}), expression(10^{1})), las=2, cex.axis=1.5);
axis(1, at=log10(Mass_R_inv), labels=c(expression(10^{-4}),"","",expression(10^{-3}),"","",expression(10^{-2})), cex.axis=1.5);
mtext("Carrying capacity", 2, line=3.5, cex=1.5)
mtext("Body mass of basal resource (g)", 1, line=2.5, cex=1.5);
abline(v=log10(1E-3), lty=3);
for(i in 2:5){lines(log10(as.numeric(K_1[i,]))~log10(Mass_R_inv), type='l', lwd=2, col=cc[i])}

# Plots along gradient of species mass of consumer
par(mfg=c(1,2), mar = c(1, 5, 2, 1));
plot(log10(as.numeric(A_2[1,]))~log10(Mass_C_inv), ylim=c(-6.3,-3.9),
     type='l', lwd=2, col=cc[1],  xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n");
axis(1, at=log10(Mass_C_inv), labels=F);
axis(2, at=log10(c(1E-6,1E-5, 1E-4)), labels=c(expression(10^{-6}), expression(10^{-5}), expression(10^{-4})), las=2, cex.axis=1.5);
mtext("Attack rate", 2, line=3.5, cex=1.5)
#mtext("Body mass of consumer (g)", 1, line=2.2, cex=1.2);
for(i in 2:5){lines(log10(as.numeric(A_2[i,]))~log10(Mass_C_inv), type='l', lwd=2, col=cc[i])}
abline(v=log10(1E-2), lty=3)

par(mfg=c(2,2), mar = c(4, 5, 0.5, 1));
plot(log10(as.numeric(H_2[1,]))~log10(Mass_C_inv), ylim=c(3.9,6),
     type='l', lwd=2, col=cc[1],  xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n");
axis(1, at=log10(Mass_C_inv), labels=c(expression(10^{-3}),"","",expression(10^{-2}),"","",expression(10^{-1})), cex.axis=1.5);
axis(2, at=log10(c(1E4,1E5,1E6)), labels=c(expression(10^{4}),expression(10^{5}),expression(10^{6})), las=2, cex.axis=1.5);
mtext("Handling time", 2, line=3.5, cex=1.5)
mtext("Body mass of consumer (g)", 1, line=2.5, cex=1.5);
abline(v=log10(1E-2), lty=3)
for(i in 2:5){lines(log10(as.numeric(H_2[i,]))~log10(Mass_C_inv), type='l', lwd=2, col=cc[i])}

par(mfg=c(1,3), mar = c(1, 5, 2, 1));
plot(log10(as.numeric(X_2[1,]))~log10(Mass_C_inv), ylim=c(-8,-5),
     type='l', lwd=2, col=cc[1],  xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n");
axis(1, at=log10(Mass_C_inv), labels=c(expression(10^{-3}),"","",expression(10^{-2}),"","",expression(10^{-1})), cex.axis=1.5);
axis(2, at=log10(c(1E-8,1E-7,1E-6,1E-5)), labels=c(expression(10^{-8}), expression(10^{-7}),expression(10^{-6}),expression(10^{-5})), las=2, cex.axis=1.5);
mtext("Metabolic loss rate", 2, line=3.5, cex=1.5)
mtext("Body mass of consumer (g)", 1, line=2.5, cex=1.5);
abline(v=log10(1E-2), lty=3);
for(i in 2:5){lines(log10(as.numeric(X_2[i,]))~log10(Mass_C_inv), type='l', lwd=2, col=cc[i])}

## Fig. S3: P* & R* Rules   ~ species mass (basal and consumer species) for varying temperatures ====
par(mfrow = c(2, 2));
# C star (P* rule)
par(mfg=c(1,1), mar = c(0, 5, 4, 1));
plot(as.numeric(P_ratio[1,])~log10(Mass_R_inv), ylim=c(0.3,1.8),
     type='l', lwd=2, col=cc[1],  xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n");
axis(2, at=c(0.5,1,1.5), labels=c(0.5, 1, 1.5), las=2, cex.axis=1.5);
axis(1, at=log10(Mass_R_inv), labels=F, cex.axis=1.5);
mtext("P* ratio", 2, line=3.1, cex=1.5)
mtext("P* Rule (AC module)", 3, line=1.5, cex=1.5);
abline(v=log10(1E-3), lty=3); abline(h=1, lty=2);
for(i in 2:5){lines(as.numeric(P_ratio[i,])~log10(Mass_R_inv), type='l', lwd=2, col=cc[i])}

par(mfg=c(2,1), mar = c(4, 5, 0, 1));
plot(as.numeric(C_eq1[1,])~log10(Mass_R_inv), ylim=c(0,0.3),
     type='l', lwd=2, col=cc[1],  xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n");
axis(1, at=log10(Mass_R_inv), labels=c(expression(10^{-4}),"","",expression(10^{-3}),"","",expression(10^{-2})), cex.axis=1.5);
axis(2, at=seq(0,0.2, by=0.1), las=2, labels=T, cex.axis=1.2);
mtext("Consumer biomass at equilibrium", 2, line=3.1, cex=1.5)
mtext("Body mass of basal resource (g)", 1, line=2.5, cex=1.5);
abline(v=log10(1E-3), lty=3);
for(i in 2:5){lines(as.numeric(C_eq1[i,])~log10(Mass_R_inv), type='l', lwd=2, col=cc[i])}

# R star (R* rule)
par(mfg=c(1,2), mar = c(0, 5, 4, 1));

plot(as.numeric(R_ratio[1,])~log10(Mass_C_inv), ylim=c(0.3,1.8),
     type='l', lwd=2, col=cc[1],  xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n");
axis(2, at=c(0.5,1,1.5), labels=c(0.5, 1, 1.5), las=2, cex.axis=1.5);
axis(1, at=log10(Mass_C_inv), labels=F, cex.axis=1.5);
mtext("R* ratio", 2, line=3.1, cex=1.5)
mtext("R* Rule (EC module)", 3, line=1.5, cex=1.5);
abline(v=log10(1E-2), lty=3); abline(h=1, lty=2);
for(i in 2:5){lines(as.numeric(R_ratio[i,])~log10(Mass_C_inv), type='l', lwd=2, col=cc[i])}

par(mfg=c(2,2), mar = c(4, 5, 0, 1));
plot(as.numeric(R_eq2[1,])~log10(Mass_C_inv), ylim=c(0,0.3),
     type='l', lwd=2, col=cc[1],  xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n");
axis(1, at=log10(Mass_C_inv), labels=c(expression(10^{-3}),"","",expression(10^{-2}),"","",expression(10^{-1})), cex.axis=1.5);
axis(2, at=seq(0,0.2, by=0.1), las=2, labels=T, cex.axis=1.2);
mtext("Resource biomass at equilibrium", 2, line=3.1, cex=1.5)
mtext("Body mass of consumer (g)", 1, line=2.5, cex=1.5);
abline(v=log10(1E-2), lty=3);
for(i in 2:5){lines(as.numeric(R_eq2[i,])~log10(Mass_C_inv), type='l', lwd=2, col=cc[i])}

######## Fig. S4: P* & R* Rules ~ species mass (basal and consumer species) for varying temperatures ====
ratio <- c(0.1, 0.2, 0.5, 2, 5, 10);
# AC ====
Inv <- read.table("./data/AC_Tot_data.txt", h=T,  sep = "\t");  # AC Complete Data
AC <- Inv[which(Inv$Alpha!=1 & Inv$Beta==10 & Inv$Tot_eq>1 & Inv$Cr_fin==1),];

# Calculus of biomass ratio in P eq in C-Ri and C-Rr
PRdata <- as.data.frame(P.Rule.Appcomp(Temp=AC$Temp[1:dim(AC)[1]], Car=AC$Car[1:dim(AC)[1]], 
                                       MRr=AC$M_Rr[1:dim(AC)[1]], MRi=AC$M_Ri[1:dim(AC)[1]], MCr=AC$M_Cr[1:dim(AC)[1]]))

P_ratio  <- PRdata$C_star_Ri/PRdata$C_star_Rr;
P_code <- rep(0, dim(PRdata)[1]);
State_code <- rep(0, dim(PRdata)[1]);

PRdata <- cbind.data.frame(PRdata[,1:5], 
                           AC$Alpha, AC$Beta, AC$Gamma, 
                           PRdata[,6:9], P_ratio, P_code, AC$State, State_code);
colnames(PRdata)<- c(colnames(PRdata)[1:5], "Alpha", "Beta", "Gamma", colnames(PRdata[,9:12]), "P_ratio", "P_code", "State", "State_code")

PRdata$P_code[which(PRdata$P_ratio < 1 )] <- 1;
PRdata$P_code[which(PRdata$P_ratio == 1)] <- 2;
PRdata$P_code[which(PRdata$P_ratio > 1)] <- 3;

for(j in 1:dim(tab_state)[1]){
  if(length(PRdata$State_code[which(PRdata$State == tab_state$State[j])]) != 0){
    PRdata$State_code[which(PRdata$State == tab_state$State[j]) ] <- tab_state$Code[j];
  }else{next}
}

# Creating the ggplots
p <- list(); q <- list();
for(i in 1:length(ratio)){
  sub <- PRdata[which(PRdata$Alpha==ratio[i]),];
  
  ylabel <- ifelse(i==1, expression("Temperature ("*degree*C*')'), "");
  
  p[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(State_code)))+theme_bw()+
    geom_tile()+
    labs(x="", y= ylabel)+
    scale_fill_manual(values = tab_state$Colour, breaks=tab_state$Code,
                      labels=tab_state$State, name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
  
  q[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(P_code)))+theme_bw()+
    geom_tile()+
    labs(x="", y= ylabel)+
    scale_fill_manual(values = col.Prule, breaks=c(1:3),
                      labels=c("P<1", "P=1", "P>1"), name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
}
# EC ====
Inv <- read.table("./data/EC_Tot_data.txt", h=T,  sep = "\t");  # EC Complete Data
EC <- Inv[which(Inv$Beta!=1 & Inv$Alpha==10 & Inv$Tot_eq>1),];
EC$Beta <- 1/EC$Beta;
# Calculus of biomass ratio at equilibrium
PRdata <- as.data.frame(R.Rule.Expcomp(Temp=EC$Temp[1:dim(EC)[1]], Car=EC$Car[1:dim(EC)[1]], 
                                       MRr=EC$M_Rr[1:dim(EC)[1]], MCi=EC$M_Ci[1:dim(EC)[1]], MCr=EC$M_Cr[1:dim(EC)[1]]))

R_ratio  <- PRdata$R_star_Ci/PRdata$R_star_Cr;
R_code <- rep(0, dim(PRdata)[1]);
State_code <- rep(0, dim(PRdata)[1]);

PRdata <- cbind.data.frame(PRdata[,1:5], 
                           EC$Alpha, EC$Beta, EC$Gamma, 
                           PRdata[,6:9], R_ratio, R_code, EC$State, State_code);
colnames(PRdata)<- c(colnames(PRdata)[1:5], "Alpha", "Beta", "Gamma", colnames(PRdata[,9:12]), "R_ratio", "R_code", "State", "State_code")

PRdata$R_code[which(PRdata$R_ratio < 1 )] <- 1;
PRdata$R_code[which(PRdata$R_ratio == 1)] <- 2;
PRdata$R_code[which(PRdata$R_ratio > 1)] <- 3;

for(j in 1:dim(tab_state)[1]){
  if(length(PRdata$State_code[which(PRdata$State == tab_state$State[j])]) != 0){
    PRdata$State_code[which(PRdata$State == tab_state$State[j]) ] <- tab_state$Code[j];
  }else{next}
}

# Creating the ggplots
r <- list(); s <- list();
for(i in 1:length(ratio)){
  sub <- PRdata[which(PRdata$Beta==ratio[i]),];
  
  ylabel <- ifelse(i==1, expression("Temperature ("*degree*C*')'), "")
  r[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(State_code)))+theme_bw()+
    geom_tile()+
    labs(x="", y= ylabel)+
    scale_fill_manual(values = tab_state$Colour, breaks=tab_state$Code,
                      labels=tab_state$State, name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
  
  s[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(R_code)))+theme_bw()+
    geom_tile()+
    labs(x=expression("Nutrient levels ("*g.L^-1*')'), y= ylabel)+
    scale_fill_manual(values = col.Rrule, breaks=c(1:3),
                      labels=c("R<1", "R=1", "R>1"), name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
}

# Arrange & Save ====
pdf(file="Fig_S4.pdf", width=20, height=11);
grid.arrange(grobs=c(p,q, r, s), ncol=6, nrow=4)
dev.off()
######## Fig. S5: R* rule & Initial biomass growth of IG-prey in IGP ====
ratio <- c(0.1, 0.2, 0.5, 2, 5, 10);
# IGPB (invasion of IGP-prey) ====
Inv1 <- read.table("./data/IGPB_Tot_data.txt", h=T,  sep = "\t");  # IGPB Complete Data
IGPB <- Inv1[which(Inv1$Alpha==10 & Inv1$Tot_eq>1),];
IGPB$Beta <- 1/IGPB$Beta;

# R Rule and initial population density growth of invader
PR_data <- as.data.frame(PR.Rules.IGP(Temp=IGPB$Temp[1:dim(IGPB)[1]], Car=IGPB$Car[1:dim(IGPB)[1]], 
                                      MRr=IGPB$M_Rr[1:dim(IGPB)[1]], MCi=IGPB$M_Ci[1:dim(IGPB)[1]], MCr=IGPB$M_Pr[1:dim(IGPB)[1]]))

R_ratio  <- PR_data$R_star_Ci/PR_data$R_star_Cr;
R_code <- rep(0, dim(PR_data)[1]);
I_code <- rep(0, dim(PR_data)[1]);
State_code <- rep(0, dim(PR_data)[1]);

PR_data <- cbind.data.frame(PR_data[,1:5], 
                            IGPB$Alpha, IGPB$Beta, IGPB$Gamma,
                            PR_data[,6:10], R_ratio, R_code, I_code, IGPB$State, State_code);
colnames(PR_data)<- c(colnames(PR_data)[1:5], "Alpha", "Beta", "Gamma",
                      colnames(PR_data[,9:13]), "R_ratio", "R_code", "I_code", "State", "State_code")

# Attributing a color code for R levels, P levels and States
PR_data$R_code[which(PR_data$R_ratio < 1 )] <- 1;
PR_data$R_code[which(PR_data$R_ratio == 1)] <- 2;
PR_data$R_code[which(PR_data$R_ratio > 1)] <- 3;

PR_data$I_code[which(PR_data$dIi < 0 )] <- 1;
PR_data$I_code[which(PR_data$dIi == 0)] <- 2;
PR_data$I_code[which(PR_data$dIi > 0)] <- 3;

for(j in 1:dim(tab_state)[1]){
  if(length(PR_data$State_code[which(PR_data$State == tab_state$State[j])]) != 0){
    PR_data$State_code[which(PR_data$State == tab_state$State[j]) ] <- tab_state$Code[j];
  }else{next}
}

# Creating the ggplots 
s <- list(); t <- list(); u <- list();
ratio <- c(0.1, 0.2, 0.5)
for(i in 1:length(ratio)){
  sub <- PR_data[which(PR_data$Beta==ratio[i]),];
  
  ylabel <- ifelse(i==1, expression("Temperature ("*degree*C*')'), "")
  s[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(State_code)))+theme_bw()+
    geom_tile()+
    labs(x="", y= ylabel)+
    scale_fill_manual(values = tab_state$Colour, breaks=tab_state$Code,
                      labels=tab_state$State, name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
  
  t[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(R_code)))+theme_bw()+
    geom_tile()+
    labs(x="", y= ylabel)+
    scale_fill_manual(values = col.Rrule, breaks=c(1:3),
                      labels=c("R<1", "R=1", "R>1"), name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
  
  u[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(I_code)))+theme_bw()+
    geom_tile()+
    labs(x=expression("Nutrient levels ("*g.L^-1*')'), y= ylabel)+
    scale_fill_manual(values = c("black", "white", "grey"), breaks=c(1:3),
                      labels=c("dI<0", "dI=0", "dI>0"), name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
}

# IGPC (invasion of IGP-predator) ====
Inv2 <- read.table("./data/IGPC_Tot_data.txt", h=T,  sep = "\t");  # IGPC Complete Data
IGPC <- Inv2[which(Inv2$Alpha==10 & Inv2$Tot_eq>1 ),]; 

# R Rule and initial population density growth of invader 
PR_data2 <- as.data.frame(PR.Rules.IGP(Temp=IGPC$Temp[1:dim(IGPC)[1]], Car=IGPC$Car[1:dim(IGPC)[1]], 
                                       MRr=IGPC$M_Rr[1:dim(IGPC)[1]], MCr=IGPC$M_Cr[1:dim(IGPC)[1]], MCi=IGPC$M_Pi[1:dim(IGPC)[1]]))

R_ratio2  <- PR_data2$R_star_Ci/PR_data2$R_star_Cr;
R_code2 <- rep(0, dim(PR_data2)[1]);
I_code2 <- rep(0, dim(PR_data2)[1]);
State_code2 <- rep(0, dim(PR_data2)[1]);

PR_data2 <- cbind.data.frame(PR_data2[,1:5], 
                             IGPC$Alpha, IGPC$Beta, IGPC$Gamma,
                             PR_data2[,c(6:9,11:12)], R_ratio2, R_code2, I_code2, IGPC$State, State_code2);
colnames(PR_data2)<- c(colnames(PR_data2)[1:5], "Alpha", "Beta", "Gamma",
                       colnames(PR_data2[,9:14]), "R_ratio", "R_code", "I_code","State", "State_code")

# Attributing a color code for R levels, P levels and States
PR_data2$R_code[which(PR_data2$R_ratio < 1 )] <- 1;
PR_data2$R_code[which(PR_data2$R_ratio == 1)] <- 2;
PR_data2$R_code[which(PR_data2$R_ratio > 1)] <- 3;

PR_data2$I_code[which(PR_data2$dIr < 0 )] <- 1;
PR_data2$I_code[which(PR_data2$dIr == 0)] <- 2;
PR_data2$I_code[which(PR_data2$dIr > 0)] <- 3;

for(j in 1:dim(tab_state)[1]){
  if(length(PR_data2$State_code[which(PR_data2$State == tab_state$State[j])]) != 0){
    PR_data2$State_code[which(PR_data2$State == tab_state$State[j]) ] <- tab_state$Code[j];
  }else{next}
}

# Creating the ggplots 
v <- list(); w <- list(); x <- list();
ratio <- c(2, 5, 10)
for(i in 1:length(ratio)){
  sub <- PR_data2[which(PR_data2$Beta==ratio[i]),];
  
  ylabel <- ifelse(i==1, expression("Temperature ("*degree*C*')'), "")
  v[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(State_code)))+theme_bw()+
    geom_tile()+
    labs(x="", y= ylabel)+
    scale_fill_manual(values = tab_state$Colour, breaks=tab_state$Code,
                      labels=tab_state$State, name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
  
  w[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(R_code)))+theme_bw()+
    geom_tile()+
    labs(x="", y= ylabel)+
    scale_fill_manual(values = col.Rrule, breaks=c(1:3),
                      labels=c("R<1", "R=1", "R>1"), name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
  
  x[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(I_code)))+theme_bw()+
    geom_tile()+
    labs(x=expression("Nutrient levels ("*g.L^-1*')'), y= ylabel)+
    scale_fill_manual(values = c("black", "white", "grey"), breaks=c(1:3),
                      labels=c("dI<0", "dI=0", "dI>0"), name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
}

# Arrange & Save ====
pdf(file="Fig_S5.pdf", width=20, height=11);
grid.arrange(grobs=c(s,v,t,w,u,x), ncol=6, nrow=3)
dev.off()

######## Fig. S6-7: Population biomass extremes along abiotic gradients ====
# Creation of the subsets from the original dataframes: selection for gradients of N=5 and T=10 ====
# TC
TC <- read.table("./data/TC_Pop.txt", sep = "\t", h=T )

N <- data.frame()
for(bmr in c(1, 2, 5, 10)){
  sub <- TC[which(TC$Carr==5 & TC$Alpha==bmr & TC$Beta==bmr & TC$Gamma==bmr^2 ),]
  N <- rbind.data.frame(N,sub)
}
N_TC <- N

T <- data.frame()
for(bmr in c(1, 2, 5, 10)){
  sub <- TC[which(TC$Temperature==10 & TC$Alpha==bmr & TC$Beta==bmr & TC$Gamma==bmr^2 ),]
  T <- rbind.data.frame(T,sub)
}
T_TC <- T

write.table(N_TC, "./data/N_TC.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)
write.table(T_TC, "./data/T_TC.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)

# EC
EC <- read.table("./data/EC_Pop.txt", sep = "\t", h=T)

N <- data.frame()
for(bmr in c(1, 2, 5, 10)){
  sub <- EC[which(EC$Carr==5 & EC$Alpha==bmr & EC$Beta==bmr & EC$Gamma==bmr^2 ),]
  N <- rbind.data.frame(N,sub)
}
N_EC <- N

T <- data.frame()
for(bmr in c(1, 2, 5, 10)){
  sub <- EC[which(EC$Temperature==10 & EC$Alpha==bmr & EC$Beta==bmr & EC$Gamma==bmr^2 ),]
  T <- rbind.data.frame(T,sub)
}
T_EC <- T

write.table(N_EC, "./data/N_EC.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)
write.table(T_EC, "./data/T_EC.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)
# AC
AC <- read.table("./data/AC_Pop.txt", sep = "\t", h=T)

N <- data.frame()
for(bmr in c(1, 2, 5, 10)){
  sub <- AC[which(AC$Carr==5 & AC$Alpha==bmr & AC$Beta==bmr & AC$Gamma==bmr^2 ),]
  N <- rbind.data.frame(N,sub)
}
N_AC <- N

T <- data.frame()
for(bmr in c(1, 2, 5, 10)){
  sub <- AC[which(AC$Temperature==10 & AC$Alpha==bmr & AC$Beta==bmr & AC$Gamma==bmr^2 ),]
  T <- rbind.data.frame(T,sub)
}
T_AC <- T

write.table(N_AC, "./data/N_AC.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)
write.table(T_AC, "./data/T_AC.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)
# IGPB
IGPB <- read.table("./data/IGPB_Pop.txt", sep = "\t", h=T)

N <- data.frame()
for(bmr in c(1, 2, 5, 10)){
  sub <- IGPB[which(IGPB$Carr==5 & IGPB$Alpha==bmr & IGPB$Beta==bmr & IGPB$Gamma==bmr^2 ),]
  N <- rbind.data.frame(N,sub)
}
N_IGPB <- N

T <- data.frame()
for(bmr in c(1, 2, 5, 10)){
  sub <- IGPB[which(IGPB$Temperature==10 & IGPB$Alpha==bmr & IGPB$Beta==bmr & IGPB$Gamma==bmr^2 ),]
  T <- rbind.data.frame(T,sub)
}
T_IGPB <- T

write.table(N_IGPB, "./data/N_IGPB.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)
write.table(T_IGPB, "./data/T_IGPB.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)
# IGPC
IGPC <- read.table("./data/IGPC_Pop.txt", sep = "\t", h=T)

N <- data.frame()
for(bmr in c(1, 2, 5, 10)){
  sub <- IGPC[which(IGPC$Carr==5 & IGPC$Alpha==bmr & IGPC$Beta==bmr & IGPC$Gamma==bmr^2 ),]
  N <- rbind.data.frame(N,sub)
}
N_IGPC <- N

T <- data.frame()
for(bmr in c(1, 2, 5, 10)){
  sub <- IGPC[which(IGPC$Temperature==10 & IGPC$Alpha==bmr & IGPC$Beta==bmr & IGPC$Gamma==bmr^2 ),]
  T <- rbind.data.frame(T,sub)
}
T_IGPC <- T

write.table(N_IGPC, "./data/N_IGPC.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)
write.table(T_IGPC, "./data/T_IGPC.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE)

######## Fig. S6: Population biomass extremes along temperature gradient for a fixed N=5 (selection for Year=10) ====
# data ====
N_AC <- read.table("./data/N_AC.txt", sep = "\t", h=T )
N_AC <- N_AC[which(N_AC$Year==10),];
N_AC$Gamma<-as.factor(N_AC$Gamma);
N_EC <- read.table("./data/N_EC.txt", sep = "\t", h=T )
N_EC <- N_EC[which(N_EC$Year==10),];
N_EC$Gamma<-as.factor(N_EC$Gamma);
N_TC <- read.table("./data/N_TC.txt", sep = "\t", h=T )
N_TC <- N_TC[which(N_TC$Year==10),];
N_TC$Gamma<-as.factor(N_TC$Gamma);
N_IGPB <- read.table("./data/N_IGPB.txt", sep = "\t", h=T )
N_IGPB <-N_IGPB[which(N_IGPB$Year==10),]
N_IGPB$Gamma<-as.factor(N_IGPB$Gamma);
N_IGPC <- read.table("./data/N_IGPC.txt", sep = "\t", h=T )
N_IGPC <-N_IGPC[which(N_IGPC$Year==10),];
N_IGPC$Gamma<-as.factor(N_IGPC$Gamma);

bmr <- c(1, 4, 25, 100);
Nut <- 5;
p <- list(); q <- list(); r <- list();
#AC ====
Yax <- expression(R[RES]);
Bax <- expression(R[INV]);
Cax <- expression(C[RES]);

p[[1]] <- ggplot(N_AC)+theme_bw()+
  geom_line(aes(x=Carr, y=A_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=A_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), labels = c(0, "", 0.5, "", 1))+
  labs(x=expression("Temperature ("*degree*C*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[1]] <- ggplot(N_AC)+theme_bw()+
  geom_line(aes(x=Carr, y=B_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=B_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .65), breaks = seq(0, .6, .3))+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

r[[1]] <- ggplot(N_AC)+theme_bw()+
  geom_line(aes(x=Carr, y=C_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=C_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .85), breaks = seq(0, .8, .4))+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");
#EC ====
Yax <- expression(R[RES])
Bax <- expression(C[INV])
Cax <- expression(C[RES])

p[[2]] <- ggplot(N_EC)+theme_bw()+
  geom_line(aes(x=Temperature, y=A_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Temperature, y=A_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), labels = c(0, "", 0.5, "", 1))+
  labs(x=expression("Temperature ("*degree*C*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[2]] <- ggplot(N_EC)+theme_bw()+
  geom_line(aes(x=Temperature, y=B_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Temperature, y=B_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .65), breaks = seq(0, .6, .3))+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

r[[2]] <- ggplot(N_EC)+theme_bw()+
  geom_line(aes(x=Temperature, y=C_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Temperature, y=C_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .85), breaks = seq(0, .8, .4))+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");
#TC ====
Yax <- expression(R[RES])
Bax <- expression(C[RES])
Cax <- expression(P[INV])

p[[3]] <- ggplot(N_TC)+theme_bw()+
  geom_line(aes(x=Temperature, y=A_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Temperature, y=A_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 2.5), labels = c(0, 2.5, 5))+
  labs(x=expression("Temperature ("*degree*C*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[3]] <- ggplot(N_TC)+theme_bw()+
  geom_line(aes(x=Temperature, y=B_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Temperature, y=B_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .65), breaks = seq(0, .6, .3))+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

r[[3]] <- ggplot(N_TC)+theme_bw()+
  geom_line(aes(x=Temperature, y=C_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Temperature, y=C_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .85), breaks = seq(0, .8, .4))+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

#IGPB ====
Yax <- expression(R[RES])
Bax <- expression(C[INV])
Cax <- expression(P[RES])

p[[4]] <- ggplot(N_IGPB)+theme_bw()+
  geom_line(aes(x=Temperature, y=A_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Temperature, y=A_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), labels = c(0, "", 0.5, "", 1))+
  labs(x=expression("Temperature ("*degree*C*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[4]] <- ggplot(N_IGPB)+theme_bw()+
  geom_line(aes(x=Temperature, y=B_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Temperature, y=B_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .65), breaks = seq(0, .6, .3))+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

r[[4]] <- ggplot(N_IGPB)+theme_bw()+
  geom_line(aes(x=Temperature, y=C_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Temperature, y=C_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .85), breaks = seq(0, .8, .4))+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");
#IGPC ====
Yax <- expression(R[RES])
Bax <- expression(C[RES])
Cax <- expression(P[INV])

p[[5]] <- ggplot(N_IGPC)+theme_bw()+
  geom_line(aes(x=Temperature, y=A_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Temperature, y=A_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), labels = c(0, "", 0.5, "", 1))+
  labs(x=expression("Temperature ("*degree*C*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[5]] <- ggplot(N_IGPC)+theme_bw()+
  geom_line(aes(x=Temperature, y=B_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Temperature, y=B_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .65), breaks = seq(0, .6, .3))+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

r[[5]] <- ggplot(N_IGPC)+theme_bw()+
  geom_line(aes(x=Temperature, y=C_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Temperature, y=C_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .85), breaks = seq(0, .8, .4))+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

#Arrange & Save ====
pdf(file="Fig_S6.pdf", width=25, height=15);
grid.arrange(grobs=c(r,q,p), ncol=5);
dev.off()
######### Fig. S7: Population biomass extremes along nutrient gradient for a fixed T=10 (selection for Year=10) ====
# data ====
T_TC <- read.table("./data/T_TC.txt", sep = "\t", h=T )
T_TC <- T_TC[which(T_TC$Year==10),];
T_TC$Gamma<-as.factor(T_TC$Gamma);
T_EC <- read.table("./data/T_EC.txt", sep = "\t", h=T )
T_EC <- T_EC[which(T_EC$Year==10),];
T_EC$Gamma<-as.factor(T_EC$Gamma);
T_AC <- read.table("./data/T_AC.txt", sep = "\t", h=T )
T_AC <- T_AC[which(T_AC$Year==10),];
T_AC$Gamma<-as.factor(T_AC$Gamma);
T_IGPB <- read.table("./data/T_IGPB.txt", sep = "\t", h=T )
T_IGPB <- T_IGPB[which(T_IGPB$Year==10),];
T_IGPB$Gamma<-as.factor(T_IGPB$Gamma);
T_IGPC <- read.table("./data/T_IGPC.txt", sep = "\t", h=T )
T_IGPC <- T_IGPC[which(T_IGPC$Year==10),];
T_IGPC$Gamma<-as.factor(T_IGPC$Gamma);
#AC ====
Yax <- expression(R[RES]);
Bax <- expression(R[INV]);
Cax <- expression(C[RES]);

p[[1]] <- ggplot(T_AC)+theme_bw()+
  geom_line(aes(x=Carr, y=A_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=A_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), labels = c(0, "", 0.5, "", 1))+
  labs(x=expression("Nutrient levels ("*g.L^-1*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[1]] <- ggplot(T_AC)+theme_bw()+
  geom_line(aes(x=Carr, y=B_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=B_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .65), breaks = seq(0, .6, .3))+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

r[[1]] <- ggplot(T_AC)+theme_bw()+
  geom_line(aes(x=Carr, y=C_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=C_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .85), breaks = seq(0, .8, .4))+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

#EC ====
Yax <- expression(R[RES])
Bax <- expression(C[INV])
Cax <- expression(C[RES])

p[[2]] <- ggplot(T_EC)+theme_bw()+
  geom_line(aes(x=Carr, y=A_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=A_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), labels = c(0, "", 0.5, "", 1))+
  labs(x=expression("Nutrient levels ("*g.L^-1*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[2]] <- ggplot(T_EC)+theme_bw()+
  geom_line(aes(x=Carr, y=B_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=B_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .65), breaks = seq(0, .6, .3))+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

r[[2]] <- ggplot(T_EC)+theme_bw()+
  geom_line(aes(x=Carr, y=C_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=C_Max, colour=Gamma), size=2)+
  labs(x="", y=Cax)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .85), breaks = seq(0, .8, .4))+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

#TC ====
Yax <- expression(R[RES])
Bax <- expression(C[RES])
Cax <- expression(P[INV])

p[[3]] <- ggplot(T_TC)+ylim(0,8)+theme_bw()+
  geom_line(aes(x=Carr, y=A_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=A_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, 4.5))+
  labs(x=expression("Nutrient levels ("*g.L^-1*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[3]] <- ggplot(T_TC)+theme_bw()+
  geom_line(aes(x=Carr, y=B_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=B_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .65), breaks = seq(0, .6, .3))+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

r[[3]] <- ggplot(T_TC)+theme_bw()+
  geom_line(aes(x=Carr, y=C_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=C_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 2.5))+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

#IGPB ====
Yax <- expression(R[RES])
Bax <- expression(C[INV])
Cax <- expression(P[RES])

p[[4]] <- ggplot(T_IGPB)+theme_bw()+
  geom_line(aes(x=Carr, y=A_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=A_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), labels = c(0, "", 0.5, "", 1))+
  labs(x=expression("Nutrient levels ("*g.L^-1*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[4]] <- ggplot(T_IGPB)+theme_bw()+
  geom_line(aes(x=Carr, y=B_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=B_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .65), breaks = seq(0, .6, .3))+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

r[[4]] <- ggplot(T_IGPB)+theme_bw()+
  geom_line(aes(x=Carr, y=C_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=C_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .85), breaks = seq(0, .8, .4))+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

#IGPC ====
Yax <- expression(R[RES])
Bax <- expression(C[RES])
Cax <- expression(P[INV])

p[[5]] <- ggplot(T_IGPC)+theme_bw()+
  geom_line(aes(x=Carr, y=A_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=A_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), labels = c(0, "", 0.5, "", 1))+
  labs(x=expression("Nutrient levels ("*g.L^-1*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[5]] <- ggplot(T_IGPC)+theme_bw()+
  geom_line(aes(x=Carr, y=B_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=B_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .65), breaks = seq(0, .6, .3))+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

r[[5]] <- ggplot(T_IGPC)+theme_bw()+
  geom_line(aes(x=Carr, y=C_Min, colour=Gamma), size=2)+
  geom_line(aes(x=Carr, y=C_Max, colour=Gamma), size=2)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  scale_y_continuous(limits = c(0, .85), breaks = seq(0, .8, .4))+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x=element_blank(),
        legend.position = "none");

#Arrange & Save ====
pdf(file="Fig_S7.pdf", width=25, height=15);
grid.arrange(grobs=c(r,q,p), ncol=5)
dev.off()
######### Fig. S8: % Invasion outcomes within environmental regions along gradients of body mass ratio across trophic modules ====
# data ====
IS_AC <- read.table("./data/BMR_AC_IS.txt", h=T,  sep = "\t");
IS_Reg_AC <- read.table("./data/BMR_AC_IS_Regions.txt", h=T,  sep = "\t");
IS_AC <- IS_AC[-which(IS_AC$Alpha==1),];
IS_Reg_AC <- IS_Reg_AC[-which(IS_Reg_AC$Alpha==1),];
colnames(IS_AC)[10] <- colnames(IS_Reg_AC)[13] <- "BMR";

IS_EC <- read.table("./data/BMR_EC_IS.txt", h=T,  sep = "\t");
IS_Reg_EC <- read.table("./data/BMR_EC_IS_Regions.txt", h=T,  sep = "\t");
IS_EC <- IS_EC[-which(IS_EC$Beta==1),];
IS_Reg_EC <- IS_Reg_EC[-which(IS_Reg_EC$Beta==1),];

colnames(IS_Reg_EC)[13] <- colnames(IS_EC)[10] <- "BMR";
IS_Reg_EC[,13] <- (1/IS_Reg_EC$BMR); # Correct bmr as 1/Beta
IS_EC[,10] <- (1/IS_EC$BMR);

IS_IGP <- read.table("./data/BMR_IGP_IS_Beta.txt", h=T,  sep = "\t");
IS_Reg_IGP <- read.table("./data/BMR_IGP_IS_Beta_Regions.txt", h=T,  sep = "\t");
colnames(IS_Reg_IGP)[13] <- colnames(IS_IGP)[10] <- "BMR";

IS_TC <- read.table("./data/BMR_TC_IS_Delta.txt", h=T,  sep = "\t");
IS_Reg_TC <- read.table("./data/BMR_TC_IS_Delta_Regions.txt", h=T,  sep = "\t");
colnames(IS_Reg_TC)[11] <- colnames(IS_TC)[8] <- "BMR";

IS_IGPC <- read.table("./data/BMR_IGPC_IS_Delta.txt", h=T,  sep = "\t");
IS_Reg_IGPC <- read.table("./data/BMR_IGPC_IS_Delta_Regions.txt", h=T,  sep = "\t");
colnames(IS_Reg_IGPC)[12] <- colnames(IS_IGPC)[9] <- "BMR";

# Fig ====
Fw <- c("AC", "EC", "IGP", "TC", "IGPC");

pdf(file="Fig_S8.pdf", width=12.5, height=15);
par(mfrow = c(5, length(Fw)));
for(r in 1:5){
  for(i in 1:length(Fw) ){
    if(r==1){dataname <- get(paste0('IS_',Fw[i]))}else{dataname <- get(paste0('IS_Reg_',Fw[i]))[which(get(paste0('IS_Reg_',Fw[i]))$Regions==LETTERS[r-1]),] }
    remcol <- ifelse(r==1, 4,7); # vector to discard n columns in data when checking only for invasion outcomes;
    leftmar <- ifelse(i==1,4,0.5); rightmar <- ifelse(i==length(Fw),4,0.5);
    downmar <- ifelse(r==5,4,0.5); upmar <- ifelse(r==1,4,0.5);
    
    par(mfg=c(r,i), mar = c(downmar, leftmar, upmar, rightmar));
    plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(-1.1,1.1), ylim=c(-5,105));
    if(r==5){axis(1, at=log10(c(0.1,1,10)), labels=c(expression(10^{-1}), 1, expression(10^{1})), cex.axis=2) }else{axis(1, at=log10(c(0.1,1,10)), labels=F)}  ;
    if(r==1){axis(3, at=log10(1), labels=Fw[i], cex.axis=2.5) }else{ axis(3, at=log10(1), labels=F)}  ;
    if(i==1){axis(2, at=c(0,25,50,75,100), label=c(0,"",50,"",100), las=2, cex.axis=2) }else{ axis(2, at=c(0,25,50,75,100), label=F, las=2)};
    abline(v=log10(1), lty=3);
    
    for(j in 1:(length(colnames(dataname))-remcol) ){
      vecname <- colnames(dataname)[j];
      if(i==5){
        lines(log10(sort(dataname$BMR,decreasing=F)), dataname[[vecname]], type="l", lwd=2,
              col=tab_state[which(tab_state$State==vecname),3], lty=tab_state[which(tab_state$State==vecname),5]);
        points(log10(sort(dataname$BMR,decreasing=F)), dataname[[vecname]], lwd=2, col=tab_state[which(tab_state$State==vecname),3], pch=tab_state[which(tab_state$State==vecname),4], cex=1.5);
      }else{
        lines(log10(dataname$BMR), dataname[[vecname]], type="l", lwd=2,
              col=tab_state[which(tab_state$State==vecname),3], lty=tab_state[which(tab_state$State==vecname),5]);
        points(log10(dataname$BMR), dataname[[vecname]], lwd=2, col=tab_state[which(tab_state$State==vecname),3], pch=tab_state[which(tab_state$State==vecname),4], cex=1.5);
      }
    }
  }
}
dev.off()

######### Fig. S9: % Changes in systeme stability within environmental regions along gradients of body mass ratio across trophic modules ====
# data ====
RS_AC <- read.table("./data/BMR_AC_RS-ResInv.txt", h=T,  sep = "\t");
RS_Reg_AC <- read.table("./data/BMR_AC_RS-ResInv_Regions.txt", h=T,  sep = "\t");
RS_AC <- RS_AC[-which(RS_AC$Alpha==1),]; RS_Reg_AC <- RS_Reg_AC[-which(RS_Reg_AC$Alpha==1),];
colnames(RS_AC)[10] <- colnames(RS_Reg_AC)[11] <- "BMR";

RS_EC <- read.table("./data/BMR_EC_RS-ResInv.txt", h=T,  sep = "\t");
RS_Reg_EC <- read.table("./data/BMR_EC_RS-ResInv_Regions.txt", h=T,  sep = "\t");
RS_EC <- RS_EC[-which(RS_EC$Beta==1),]; RS_Reg_EC <- RS_Reg_EC[-which(RS_Reg_EC$Beta==1),]; 
colnames(RS_EC)[13] <- colnames(RS_Reg_EC)[14] <- "BMR";
RS_EC[,13] <- (1/RS_EC$BMR); # Correct bmr as 1/Beta
RS_Reg_EC[,14] <- (1/RS_Reg_EC$BMR);

RS_IGP <- read.table("./data/BMR_IGP_RS-ResInv_Beta.txt", h=T,  sep = "\t");
RS_Reg_IGP <- read.table("./data/BMR_IGP_RS-ResInv_Beta_Regions.txt", h=T,  sep = "\t");
colnames(RS_IGP)[12] <- colnames(RS_Reg_IGP)[12] <- "BMR";

RS_TC <- read.table("./data/BMR_TC_RS-ResInv_Delta.txt", h=T,  sep = "\t");
RS_Reg_TC <- read.table("./data/BMR_TC_RS-ResInv_Delta_Regions.txt", h=T,  sep = "\t");
colnames(RS_TC)[13] <- colnames(RS_Reg_TC)[14] <- "BMR";

RS_IGPC <- read.table("./data/BMR_IGPC_RS-ResInv_Delta.txt", h=T,  sep = "\t");
RS_Reg_IGPC <- read.table("./data/BMR_IGPC_RS-ResInv_Delta_Regions.txt", h=T,  sep = "\t");
colnames(RS_IGPC)[11] <- colnames(RS_Reg_IGPC)[12] <- "BMR";

# Fig ====
Fw <- c("AC", "EC", "IGP", "TC", "IGPC");

pdf(file="Fig_S9.pdf", width=12.5, height=15);
par(mfrow = c(5, length(Fw)));
for(r in 1:5){
  for(i in 1:length(Fw) ){
    if(r==1){dataname <- get(paste0('RS_',Fw[i]))}else{dataname <- get(paste0('RS_Reg_',Fw[i]))[which(get(paste0('RS_Reg_',Fw[i]))$Regions==LETTERS[r-1]),] }
    remcol <- ifelse(r==1, 4,5); # vector to discard n columns in data when checking only for invasion outcomes;
    leftmar <- ifelse(i==1,4,0.5); rightmar <- ifelse(i==length(Fw),4,0.5);
    downmar <- ifelse(r==5,4,0.5); upmar <- ifelse(r==1,4,0.5);
    
    par(mfg=c(r,i), mar = c(downmar, leftmar, upmar, rightmar));
    plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(-1.1,1.1), ylim=c(-5,105));
    if(r==5){axis(1, at=log10(c(0.1,1,10)), labels=c(expression(10^{-1}), 1, expression(10^{1})), cex.axis=2) }else{axis(1, at=log10(c(0.1,1,10)), labels=F)}  ;
    if(r==1){axis(3, at=log10(1), labels=Fw[i], cex.axis=2.5) }else{ axis(3, at=log10(1), labels=F)}  ;
    if(i==1){axis(2, at=c(0,25,50,75,100), label=c(0,"",50,"",100), las=2, cex.axis=2) }else{ axis(2, at=c(0,25,50,75,100), label=F, las=2)};
    abline(v=log10(1), lty=3);
    
    for(j in 1:(length(colnames(dataname))-remcol) ){
      vecname <- colnames(dataname)[j];
      if(i==5){
        lines(log10(sort(dataname$BMR,decreasing=F)), dataname[[vecname]], type="l", lwd=2,
              col=tab_stab_ResInv[which(tab_stab_ResInv$Stab==vecname),4], lty=tab_stab_ResInv[which(tab_stab_ResInv$Stab==vecname),7]);
        points(log10(sort(dataname$BMR,decreasing=F)), dataname[[vecname]], lwd=2, col=tab_stab_ResInv[which(tab_stab_ResInv$Stab==vecname),4], pch=tab_stab_ResInv[which(tab_stab_ResInv$Stab==vecname),5], cex=1.5);
      }else{
        lines(log10(dataname$BMR), dataname[[vecname]], type="l", lwd=2,
              col=tab_stab_ResInv[which(tab_stab_ResInv$Stab==vecname),4], lty=tab_stab_ResInv[which(tab_stab_ResInv$Stab==vecname),7]);
        points(log10(dataname$BMR), dataname[[vecname]], lwd=2, col=tab_stab_ResInv[which(tab_stab_ResInv$Stab==vecname),4], pch=tab_stab_ResInv[which(tab_stab_ResInv$Stab==vecname),5], cex=1.5);
      }
    }
  }
}
dev.off()

######### Fig. S10: % Observed Regime States per invasion outcomes ====
# Please go back to Fig. 4 to load the required table 
# Barplot ====
p <- ggplot(data=Tab_obs_State_Regime, aes(x=Mechanisms, y=P_obs, color=Reg))+
  geom_bar(stat="identity", color="black", fill=Tab_obs_State_Regime$Colours)+
  labs(x="Invasion mechanisms", y= expression("Regime states"~S['RES']*.*S['INV']~"(%)"))+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15), legend.position="none")

# Save
pdf(file="Fig_S10.pdf", width=8, height=7);
p;
dev.off()

######### Fig. S11: Barplots of Invasion mechanisms & Regime states 1) Overall observation; 2) across food web modules ====
# with (S11a) or without (S11b) Community's resistance
# data ====
Inv.AC <- read.table("./data/AC_Tot_data.txt", h=T,  sep = "\t");  # AC Complete Data
Inv.EC <- read.table("./data/EC_Tot_data.txt", h=T,  sep = "\t");  # EC Complete Data
Inv.TC <- read.table("./data/TC_Tot_data.txt", h=T,  sep = "\t");  # TC Complete Data
Inv.IGPB <- read.table("./data/IGPB_Tot_data.txt", h=T,  sep = "\t");  # IGPB Complete Data
Inv.IGPC <- read.table("./data/IGPC_Tot_data.txt", h=T,  sep = "\t");  # IGPC Complete Data
# 2 version of selection that either (S2) include Resistance (invasion failure) or (S3) exclude it from the observations
# Version S11a: Including Resistance in the observations ====
AC2 <-Inv.AC[which(Inv.AC$Alpha!=1),];
EC2 <- Inv.EC[which(Inv.EC$Beta!=1),];
AC <- AC2
EC <- EC2
TC <- Inv.TC
IGPB <- Inv.IGPB
IGPC <- Inv.IGPC
# Amount of observations across all topologies and in total
N_Obs_AC <- dim(AC2)[1];
N_Obs_EC <- dim(EC2)[1];
N_Obs_TC <- dim(Inv.TC)[1];
N_Obs_IGPB <- dim(Inv.IGPB)[1];
N_Obs_IGPC <- dim(Inv.IGPC)[1];
N_Obs <- N_Obs_AC+N_Obs_EC+N_Obs_TC+N_Obs_IGPC+N_Obs_IGPB;
State <- tab_state$State;
Inv.Col <- tab_state$Colours
# Version S11b: Excluding Resistance from the observations ====
AC1 <- Inv.AC[which(Inv.AC$Alpha!=1 & Inv.AC$State != "Resistance"),]; # Removing Alpha=1 & all resistance cases to invasion
EC1 <- Inv.EC[which(Inv.EC$Beta!=1 & Inv.EC$State != "Resistance"),]; # Removing Beta=1  & all resistance cases to invasion
AC <- AC1
EC <- EC1
TC <- Inv.TC[which(Inv.TC$State != "Resistance"),]; # Removing all resistance cases to invasion
IGPB <- Inv.IGPB[which(Inv.IGPB$State != "Resistance"),]; # Removing all resistance cases to invasion
IGPC <- Inv.IGPC[which(Inv.IGPC$State != "Resistance"),]; # Removing all resistance cases to invasion
# Amount of observations across all topologies and in total
N_Obs_AC <- dim(AC)[1];
N_Obs_EC <- dim(EC)[1];
N_Obs_TC <- dim(TC)[1];
N_Obs_IGPB <- dim(IGPB)[1];
N_Obs_IGPC <- dim(IGPC)[1];
N_Obs <- N_Obs_AC+N_Obs_EC+N_Obs_TC+N_Obs_IGPC+N_Obs_IGPB;
State <- tab_state$State[c(1,3:6)];
Inv.Col <- tab_state$Colours[c(1,3:6)];

### Invasion outcomes ====
# Calculus of number & proportions of observed mechanisms in each food web
# AC
Tab_obs_AC <- cbind.data.frame(State, rep('AC', length(State)), rep(0, length(State)), rep(0, length(State)), Inv.Col);
colnames(Tab_obs_AC) <- c("State", "Fw", "N_obs", "P_obs", "Colours");
for(i in 1:length(State)){
  Tab_obs_AC[i,3] <- dim(AC[AC$State==State[i],] )[1];
  Tab_obs_AC[i,4] <- Proportion(Tab_obs_AC[i,3], N_Obs_AC)
}

# EC
Tab_obs_EC <- cbind.data.frame(State, rep('EC', length(State)), rep(0, length(State)), rep(0, length(State)), Inv.Col);
colnames(Tab_obs_EC) <- c("State", "Fw", "N_obs", "P_obs", "Colours");
for(i in 1:length(State)){
  Tab_obs_EC[i,3] <- dim(EC[EC$State==State[i],] )[1];
  Tab_obs_EC[i,4] <- Proportion(Tab_obs_EC[i,3], N_Obs_EC)
}

# TC
Tab_obs_TC <- cbind.data.frame(State, rep('TC', length(State)), rep(0, length(State)), rep(0, length(State)), Inv.Col);
colnames(Tab_obs_TC) <- c("State", "Fw", "N_obs", "P_obs", "Colours");
for(i in 1:length(State)){
  Tab_obs_TC[i,3] <- dim(TC[TC$State==State[i],] )[1];
  Tab_obs_TC[i,4] <- Proportion(Tab_obs_TC[i,3], N_Obs_TC)
}

# IGPB
Tab_obs_IGPB <- cbind.data.frame(State, rep('IGPB', length(State)), rep(0, length(State)), rep(0, length(State)), Inv.Col);
colnames(Tab_obs_IGPB) <- c("State", "Fw", "N_obs", "P_obs", "Colours");
for(i in 1:length(State)){
  Tab_obs_IGPB[i,3] <- dim(IGPB[IGPB$State==State[i],] )[1];
  Tab_obs_IGPB[i,4] <- Proportion(Tab_obs_IGPB[i,3], N_Obs_IGPB)
}

# IGPC
Tab_obs_IGPC <- cbind.data.frame(State, rep('IGPC', length(State)), rep(0, length(State)), rep(0, length(State)), Inv.Col);
colnames(Tab_obs_IGPC) <- c("State", "Fw", "N_obs", "P_obs", "Colours");
for(i in 1:length(State)){
  Tab_obs_IGPC[i,3] <- dim(IGPC[IGPC$State==State[i],] )[1];
  Tab_obs_IGPC[i,4] <- Proportion(Tab_obs_IGPC[i,3], N_Obs_IGPC)
}

# Merging the different tables
Tab_obs_Fw <- rbind.data.frame(Tab_obs_AC, Tab_obs_EC, Tab_obs_TC, Tab_obs_IGPB, Tab_obs_IGPC);

# Calculus of all observations across all fw
Tab_obs <- cbind.data.frame(State, rep(0, length(State)), rep(0, length(State)), Inv.Col);
colnames(Tab_obs) <- c("State", "N_obs", "P_obs", "Colours");

for(i in 1:length(State)){
  Tab_obs[i,2] <- sum(Tab_obs_Fw$N_obs[Tab_obs_Fw$State==State[i]]) ;
  Tab_obs[i,3] <- Proportion(Tab_obs[i,2], N_Obs)
}

### Regime states ====
# Calculus of number & proportions of observed Resgime states in each food web
# AC
Tab_obs_RS_AC <- cbind.data.frame(tab_stab_ResInv$Stab, tab_stab_ResInv$State, rep('AC', 9), rep(0, 9), rep(0, 9), tab_stab_ResInv$Colours);
colnames(Tab_obs_RS_AC) <- c("Stab", "State", "Fw", "N_obs", "P_obs", "Colours");
for(i in 1:dim(Tab_obs_RS_AC)[1]){
  Tab_obs_RS_AC[i,4] <- dim(AC[AC$EQ_transit_ResInv==Tab_obs_RS_AC$Stab[i],] )[1];
  Tab_obs_RS_AC[i,5] <- Proportion(Tab_obs_RS_AC[i,4], N_Obs_AC)
}

# EC
Tab_obs_RS_EC <- cbind.data.frame(tab_stab_ResInv$Stab, tab_stab_ResInv$State, rep('EC', 9), rep(0, 9), rep(0, 9), tab_stab_ResInv$Colours);
colnames(Tab_obs_RS_EC) <- c("Stab", "State", "Fw", "N_obs", "P_obs", "Colours");
for(i in 1:dim(Tab_obs_RS_EC)[1]){
  Tab_obs_RS_EC[i,4] <- dim(EC[EC$EQ_transit_ResInv==Tab_obs_RS_EC$Stab[i],] )[1];
  Tab_obs_RS_EC[i,5] <- Proportion(Tab_obs_RS_EC[i,4], N_Obs_EC)
}

# TC
Tab_obs_RS_TC <- cbind.data.frame(tab_stab_ResInv$Stab, tab_stab_ResInv$State, rep('TC', 9), rep(0, 9), rep(0, 9), tab_stab_ResInv$Colours);
colnames(Tab_obs_RS_TC) <- c("Stab", "State", "Fw", "N_obs", "P_obs", "Colours");
for(i in 1:dim(Tab_obs_RS_TC)[1]){
  Tab_obs_RS_TC[i,4] <- dim(TC[TC$EQ_transit_ResInv==Tab_obs_RS_TC$Stab[i],] )[1];
  Tab_obs_RS_TC[i,5] <- Proportion(Tab_obs_RS_TC[i,4], N_Obs_TC)
}

# IGPB
Tab_obs_RS_IGPB <- cbind.data.frame(tab_stab_ResInv$Stab, tab_stab_ResInv$State, rep('IGPB', 9), rep(0, 9), rep(0, 9), tab_stab_ResInv$Colours);
colnames(Tab_obs_RS_IGPB) <- c("Stab", "State", "Fw", "N_obs", "P_obs", "Colours");
for(i in 1:dim(Tab_obs_RS_IGPB)[1]){
  Tab_obs_RS_IGPB[i,4] <- dim(IGPB[IGPB$EQ_transit_ResInv==Tab_obs_RS_IGPB$Stab[i],] )[1];
  Tab_obs_RS_IGPB[i,5] <- Proportion(Tab_obs_RS_IGPB[i,4], N_Obs_IGPB)
}

# IGPC
Tab_obs_RS_IGPC <- cbind.data.frame(tab_stab_ResInv$Stab, tab_stab_ResInv$State, rep('IGPC', 9), rep(0, 9), rep(0, 9), tab_stab_ResInv$Colours);
colnames(Tab_obs_RS_IGPC) <- c("Stab", "State", "Fw", "N_obs", "P_obs", "Colours");
for(i in 1:dim(Tab_obs_RS_IGPC)[1]){
  Tab_obs_RS_IGPC[i,4] <- dim(IGPC[IGPC$EQ_transit_ResInv==Tab_obs_RS_IGPC$Stab[i],] )[1];
  Tab_obs_RS_IGPC[i,5] <- Proportion(Tab_obs_RS_IGPC[i,4], N_Obs_IGPC)
}
# Merging the different tables
Tab_obs_RS_Fw <- rbind.data.frame(Tab_obs_RS_AC, Tab_obs_RS_EC, Tab_obs_RS_TC, Tab_obs_RS_IGPB, Tab_obs_RS_IGPC);

# Calculus of all observations across all fw
Tab_obs_RS <- cbind.data.frame(tab_stab_ResInv$Stab, tab_stab_ResInv$State, rep(0, 9), rep(0, 9), tab_stab_ResInv$Colours);
colnames(Tab_obs_RS) <- c("Stab", "State", "N_obs", "P_obs", "Colours");

for(i in 1:9){
  Tab_obs_RS[i,3] <- sum(Tab_obs_RS_Fw$N_obs[Tab_obs_RS_Fw$Stab==tab_stab_ResInv$Stab[i]]) ;
  Tab_obs_RS[i,4] <- Proportion(Tab_obs_RS[i,3], N_Obs)
}

### Barplots ====
p <- ggplot(data=Tab_obs, aes(x=State, y=P_obs, color=State))+
  geom_bar(stat="identity", color="black", fill=Tab_obs$Colours)+
  labs(x="Invasion mechanisms", y= "Percentage of observations (%)")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15), legend.position="none")

q <- ggplot(data=Tab_obs_Fw, aes(x=Fw, y=P_obs, color=State))+
  geom_bar(stat="identity", color="black", fill=Tab_obs_Fw$Colours)+
  labs(x="Modules", y= "Percentage of observations (%)")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15), legend.position="none")

r <- ggplot(data=Tab_obs_RS, aes(x=State, y=P_obs, color=State))+
  geom_bar(stat="identity", color="black", fill=Tab_obs_RS$Colours)+
  labs(x="Regime states", y= "Percentage of observations (%)")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15), legend.position="none")

s <- ggplot(data=Tab_obs_RS_Fw, aes(x=Fw, y=P_obs, color=State))+
  geom_bar(stat="identity", color="black", fill=Tab_obs_RS_Fw$Colours)+
  labs(x="Modules", y= "Percentage of observations (%)")+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15), legend.position="none")

### Arrange & Save ====
pdf(file="Fig_S11.pdf", width=16, height=12);
#pdf(file="Fig_S11b.pdf", width=16, height=12);
grid.arrange(p,q,r,s, nrow=2, ncol=2, layout_matrix=rbind(c(1,2),c(3,4)))
dev.off()