########################################################################
### R Script                     Figures.R
###
###   Codes used to obtain the selected figures
########################################################################
#Required packages
library(deSolve)
library(rootSolve)
library(ggplot2)
library(lattice)
library(gridExtra)

#library(devtools)
#devtools::install_github("ropenscilabs/ochRe")
library(RColorBrewer)
library(ochRe)

source("./code/Functions.R");
Fw <- c("AC", "EC", "TC", "IGPB", "IGPC");
########################################################################
####################### Color selection used across the Figures =====
##### Palette for Biodiversity change
col.bio <- c(brewer.pal(n=11, name='RdYlBu')[c(10,8)],  # -2, -1
             brewer.pal(n=9, name='Greys')[4],          # 0
             brewer.pal(n=11, name='RdYlBu')[c(5,3,1)]);# +1, +2, +3
tab_bio <-  data.frame(matrix(0, ncol=2, nrow=length(col.bio)));
colnames(tab_bio) <- c("Bio", "Colours");
tab_bio$Bio <- c(-2, -1, 0, 1, 2, 3);
tab_bio$Colours <- col.bio;
atpar_bio = seq(-2.5, 3.5, by=1); seq_bio = seq(-2, 3, by=1);

##### Palette for Community composition
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

##### Palette for Invasion states
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

##### Palette for Changes in stability regime
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

##### Palette To illustrate the outcomes of R* and P* rules ====
col.Rrule <- c(brewer.pal(n=11, name='PuOr')[3],
               brewer.pal(n=9, name='Greys')[6],
               brewer.pal(n=11, name='PuOr')[9]);

col.Prule <- c(brewer.pal(n=11, name='PuOr')[9],
               brewer.pal(n=9, name='Greys')[6],
               brewer.pal(n=11, name='PuOr')[3]);
##### Palette for the Populations Biomass extremes (Figs S10-11) ====
pal <- as.vector(unlist(ochre_palettes[8]));
col.bmr <- c(pal[2], # bmr=1
             pal[4], # bmr=4
             pal[1], # bmr=25
             pal[3]);# bmr=100
####################### FIGURE 1: Levelplots of Invasion mechanisms and Regimes states (S_res.S_inv) for Gamma==100 (Alpha= Beta= 10) ====
# AC ====
name <- paste(paste("AC","Tot", "data", sep='_'),"txt", sep='.')
F1 <- read.table( ./data/name, h=T,  sep = "\t");
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
F1 <- read.table( paste(paste("EC","Tot", "data", sep='_'),"txt", sep='.'), h=T,  sep = "\t");
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
F1 <- read.table( paste(paste("TC","Tot", "data", sep='_'),"txt", sep='.'), h=T,  sep = "\t");
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
F1 <- read.table( paste(paste("IGPB","Tot", "data", sep='_'),"txt", sep='.'), h=T,  sep = "\t");
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
F1 <- read.table( paste(paste("IGPC","Tot", "data", sep='_'),"txt", sep='.'), h=T,  sep = "\t");
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
pdf(file="./plots/Fig_1.pdf");
grid.arrange(p1, p3, p5, p7, p9, p2, p4, p6, p8, p10, ncol=5, nrow=2, newpage = T)
dev.off()
####################### FIGURE 2: GGplots Invasion mechanisms along BMR gradients ====
# AC ====
Res.bmr.is <- read.table("./data/BMR_AC_IS.txt", h=T,  sep = "\t");
Res.bmr.is <- Res.bmr.is[,-1];

p1 <- ggplot(Res.bmr.is)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Alpha, y=Resistance), colour=tab_state$Colours[2], size=1.2)+
  geom_line(aes(x=Alpha, y=Substitution), colour=tab_state$Colours[3], size=1.2)+
  geom_line(aes(x=Alpha, y=Occupancy), colour=tab_state$Colours[5], size=1.2)+
  geom_point(aes(x=Alpha, y=Resistance), colour=tab_state$Colours[2], shape=tab_state$Symbols[2], size=2.5)+
  geom_point(aes(x=Alpha, y=Substitution), colour=tab_state$Colours[3], shape=tab_state$Symbols[3], size=2.5)+
  geom_point(aes(x=Alpha, y=Occupancy), colour=tab_state$Colours[5], shape=tab_state$Symbols[5], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3)

# EC ====
Res.bmr.is <- read.table("./data/BMR_EC_IS.txt", h=T,  sep = "\t");
Res.bmr.is <- Res.bmr.is[,-1]; Res.bmr.is$Beta <- (1/Res.bmr.is$Beta);

p2 <- ggplot(Res.bmr.is)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Beta, y=Vulnerability), colour=tab_state$Colours[1], size=1.2)+
  geom_line(aes(x=Beta, y=Resistance), colour=tab_state$Colours[2], size=1.2)+
  geom_line(aes(x=Beta, y=Substitution), colour=tab_state$Colours[3], size=1.2)+
  geom_line(aes(x=Beta, y=Occupancy), colour=tab_state$Colours[5], size=1.2)+
  geom_point(aes(x=Beta, y=Vulnerability), colour=tab_state$Colours[1], shape=tab_state$Symbols[1], size=2.5)+
  geom_point(aes(x=Beta, y=Resistance), colour=tab_state$Colours[2], shape=tab_state$Symbols[2], size=2.5)+
  geom_point(aes(x=Beta, y=Substitution), colour=tab_state$Colours[3], shape=tab_state$Symbols[3], size=2.5)+
  geom_point(aes(x=Beta, y=Occupancy), colour=tab_state$Colours[5], shape=tab_state$Symbols[5], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3)

# IGP ====
Res.bmr.is <- read.table("./data/BMR_IGP_IS_Beta.txt", h=T,  sep = "\t");

p3 <- ggplot(Res.bmr.is)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Beta, y=Resistance), colour=tab_state$Colours[2], size=1.2)+
  geom_line(aes(x=Beta, y=Substitution), colour=tab_state$Colours[3], size=1.2)+
  geom_line(aes(x=Beta, y=Occupancy), colour=tab_state$Colours[5], size=1.2)+
  geom_line(aes(x=Beta, y=Integration), colour=tab_state$Colours[6], size=1.2)+
  geom_point(aes(x=Beta, y=Resistance), colour=tab_state$Colours[2], shape=tab_state$Symbols[2], size=2.5)+
  geom_point(aes(x=Beta, y=Substitution), colour=tab_state$Colours[3], shape=tab_state$Symbols[3], size=2.5)+
  geom_point(aes(x=Beta, y=Occupancy), colour=tab_state$Colours[5], shape=tab_state$Symbols[5], size=2.5)+
  geom_point(aes(x=Beta, y=Integration), colour=tab_state$Colours[6], shape=tab_state$Symbols[6], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3)

# TC ==== 
Res.bmr.is.gam <- read.table("./data/BMR_TC_IS_Gamma.txt", h=T,  sep = "\t");
Res.bmr.is.del <- read.table("./data/BMR_TC_IS_Delta.txt", h=T,  sep = "\t");

p4 <- ggplot(Res.bmr.is.gam)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=Resistance), colour=tab_state$Colours[2], size=1.2)+
  geom_line(aes(x=Gamma, y=Rescue), colour=tab_state$Colours[4], size=1.2)+
  geom_line(aes(x=Gamma, y=Integration), colour=tab_state$Colours[6], size=1.2)+
  geom_point(aes(x=Gamma, y=Resistance), colour=tab_state$Colours[2], shape=tab_state$Symbols[2], size=2.5)+
  geom_point(aes(x=Gamma, y=Rescue), colour=tab_state$Colours[4], shape=tab_state$Symbols[4], size=2.5)+
  geom_point(aes(x=Gamma, y=Integration), colour=tab_state$Colours[6], shape=tab_state$Symbols[6], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(1,2.5,5,10,25,50,100), labels=c(expression(10^0),"","",1,"","",expression(10^2)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))

p5 <- ggplot(Res.bmr.is.del)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Ddelta, y=Resistance), colour=tab_state$Colours[2], size=1.2)+
  geom_line(aes(x=Ddelta, y=Rescue), colour=tab_state$Colours[4], size=1.2)+
  geom_line(aes(x=Ddelta, y=Integration), colour=tab_state$Colours[6], size=1.2)+
  geom_point(aes(x=Ddelta, y=Resistance), colour=tab_state$Colours[2], shape=tab_state$Symbols[2], size=2.5)+
  geom_point(aes(x=Ddelta, y=Rescue), colour=tab_state$Colours[4], shape=tab_state$Symbols[4], size=2.5)+
  geom_point(aes(x=Ddelta, y=Integration), colour=tab_state$Colours[6], shape=tab_state$Symbols[6], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3)

# IGPC ====
Res.bmr.is.gam <- read.table("./data/BMR_IGPC_IS_Gamma.txt", h=T,  sep = "\t");
Res.bmr.is.del <- read.table("./data/BMR_IGPC_IS_Delta.txt", h=T,  sep = "\t");

p6 <- ggplot(Res.bmr.is.gam)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=Resistance), colour=tab_state$Colours[2], size=1.2)+
  geom_line(aes(x=Gamma, y=Substitution), colour=tab_state$Colours[3], size=1.2)+
  geom_line(aes(x=Gamma, y=Occupancy), colour=tab_state$Colours[5], size=1.2)+
  geom_line(aes(x=Gamma, y=Integration), colour=tab_state$Colours[6], size=1.2)+
  geom_point(aes(x=Gamma, y=Resistance), colour=tab_state$Colours[2], shape=tab_state$Symbols[2], size=2.5)+
  geom_point(aes(x=Gamma, y=Substitution), colour=tab_state$Colours[3], shape=tab_state$Symbols[3], size=2.5)+
  geom_point(aes(x=Gamma, y=Occupancy), colour=tab_state$Colours[5], shape=tab_state$Symbols[5], size=2.5)+
  geom_point(aes(x=Gamma, y=Integration), colour=tab_state$Colours[6], shape=tab_state$Symbols[6], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(1,2.5,5,10,25,50,100), labels=c(expression(10^0),"","",expression(10^1),"","",expression(10^2)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))

p7 <- ggplot(Res.bmr.is.del)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Ddelta, y=Resistance), colour=tab_state$Colours[2], size=1.2)+
  geom_line(aes(x=Ddelta, y=Substitution), colour=tab_state$Colours[3], size=1.2)+
  geom_line(aes(x=Ddelta, y=Occupancy), colour=tab_state$Colours[5], size=1.2)+
  geom_line(aes(x=Ddelta, y=Integration), colour=tab_state$Colours[6], size=1.2)+
  geom_point(aes(x=Ddelta, y=Resistance), colour=tab_state$Colours[2], shape=tab_state$Symbols[2], size=2.5)+
  geom_point(aes(x=Ddelta, y=Substitution), colour=tab_state$Colours[3], shape=tab_state$Symbols[3], size=2.5)+
  geom_point(aes(x=Ddelta, y=Occupancy), colour=tab_state$Colours[5], shape=tab_state$Symbols[5], size=2.5)+
  geom_point(aes(x=Ddelta, y=Integration), colour=tab_state$Colours[6], shape=tab_state$Symbols[6], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3)

# Arrange & Save ====
pdf(file="./plots/Fig_2.pdf");
grid.arrange(p1,p2,p3,p5,p7 ,nrow=2, ncol=3,
             layout_matrix=rbind(c(1,2,3),c(4,5,NA)))
dev.off()
####################### FIGURE 3: GGplots Regime State along BMR gradients ====
# AC ====
Res.bmr.ri <- read.table("./data/BMR_AC_RS-ResInv.txt", h=T,  sep = "\t");

p1 <- ggplot(Res.bmr.ri)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Alpha, y=Null.Null), colour=tab_stab_ResInv$Colours[1], linetype=tab_stab_ResInv$Linetypes[1], size=1.2)+
  geom_line(aes(x=Alpha, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], linetype=tab_stab_ResInv$Linetypes[3], size=1.2)+
  geom_line(aes(x=Alpha, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], linetype=tab_stab_ResInv$Linetypes[5], size=1.2)+
  geom_line(aes(x=Alpha, y=Cycles.Eq), colour=tab_stab_ResInv$Colours[8], linetype=tab_stab_ResInv$Linetypes[8], size=1.2)+
  geom_line(aes(x=Alpha, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], linetype=tab_stab_ResInv$Linetypes[9], size=1.2)+
  geom_point(aes(x=Alpha, y=Null.Null), colour=tab_stab_ResInv$Colours[1], shape=tab_stab_ResInv$Symbols[1], size=2.5)+
  geom_point(aes(x=Alpha, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], shape=tab_stab_ResInv$Symbols[3], size=2.5)+
  geom_point(aes(x=Alpha, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], shape=tab_stab_ResInv$Symbols[5], size=2.5)+
  geom_point(aes(x=Alpha, y=Cycles.Eq), colour=tab_stab_ResInv$Colours[8], shape=tab_stab_ResInv$Symbols[8], size=2.5)+
  geom_point(aes(x=Alpha, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], shape=tab_stab_ResInv$Symbols[9], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3)

# EC ====
Res.bmr.ri <- read.table("./data/BMR_EC_RS-ResInv.txt", h=T,  sep = "\t");
Res.bmr.ri$Beta <- (1/Res.bmr.ri$Beta);

p2 <- ggplot(Res.bmr.ri)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Beta, y=Null.Null), colour=tab_stab_ResInv$Colours[1], linetype=tab_stab_ResInv$Linetypes[1], size=1.2)+
  geom_line(aes(x=Beta, y=Eq.Null), colour=tab_stab_ResInv$Colours[4], linetype=tab_stab_ResInv$Linetypes[4], size=1.2)+
  geom_line(aes(x=Beta, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], linetype=tab_stab_ResInv$Linetypes[5], size=1.2)+
  geom_line(aes(x=Beta, y=Eq.Cycles), colour=tab_stab_ResInv$Colours[6], linetype=tab_stab_ResInv$Linetypes[6], size=1.2)+
  geom_line(aes(x=Beta, y=Cycles.Null), colour=tab_stab_ResInv$Colours[7], linetype=tab_stab_ResInv$Linetypes[7], size=1.2)+
  geom_line(aes(x=Beta, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], linetype=tab_stab_ResInv$Linetypes[9], size=1.2)+
  geom_point(aes(x=Beta, y=Null.Null), colour=tab_stab_ResInv$Colours[1], shape=tab_stab_ResInv$Symbols[1], size=2.5)+
  geom_point(aes(x=Beta, y=Eq.Null), colour=tab_stab_ResInv$Colours[4], shape=tab_stab_ResInv$Symbols[4], size=2.5)+
  geom_point(aes(x=Beta, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], shape=tab_stab_ResInv$Symbols[5], size=2.5)+
  geom_point(aes(x=Beta, y=Eq.Cycles), colour=tab_stab_ResInv$Colours[6], shape=tab_stab_ResInv$Symbols[6], size=2.5)+
  geom_point(aes(x=Beta, y=Cycles.Null), colour=tab_stab_ResInv$Colours[7], shape=tab_stab_ResInv$Symbols[7], size=2.5)+
  geom_point(aes(x=Beta, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], shape=tab_stab_ResInv$Symbols[9], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3)

# IGP ====
Res.bmr.ri <- read.table("./data/BMR_IGP_RS-ResInv_Beta.txt", h=T,  sep = "\t");

p3 <- ggplot(Res.bmr.ri)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Beta, y=Null.Null), colour=tab_stab_ResInv$Colours[1], linetype=tab_stab_ResInv$Linetypes[1], size=1.2)+
  geom_line(aes(x=Beta, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], linetype=tab_stab_ResInv$Linetypes[2], size=1.2)+
  geom_line(aes(x=Beta, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], linetype=tab_stab_ResInv$Linetypes[3], size=1.2)+
  geom_line(aes(x=Beta, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], linetype=tab_stab_ResInv$Linetypes[5], size=1.2)+
  geom_line(aes(x=Beta, y=Cycles.Eq), colour=tab_stab_ResInv$Colours[8], linetype=tab_stab_ResInv$Linetypes[8], size=1.2)+
  geom_line(aes(x=Beta, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], linetype=tab_stab_ResInv$Linetypes[9], size=1.2)+
  geom_point(aes(x=Beta, y=Null.Null), colour=tab_stab_ResInv$Colours[1], shape=tab_stab_ResInv$Symbols[1], size=2.5)+
  geom_point(aes(x=Beta, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], shape=tab_stab_ResInv$Symbols[2], size=2.5)+
  geom_point(aes(x=Beta, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], shape=tab_stab_ResInv$Symbols[3], size=2.5)+
  geom_point(aes(x=Beta, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], shape=tab_stab_ResInv$Symbols[5], size=2.5)+
  geom_point(aes(x=Beta, y=Cycles.Eq), colour=tab_stab_ResInv$Colours[8], shape=tab_stab_ResInv$Symbols[8], size=2.5)+
  geom_point(aes(x=Beta, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], shape=tab_stab_ResInv$Symbols[9], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3)

# TC ====
Res.bmr.ri.gam <- read.table("./data/BMR_TC_RS-ResInv_Gamma.txt", h=T,  sep = "\t");
Res.bmr.ri.del <- read.table("./data/BMR_TC_RS-ResInv_Delta.txt", h=T,  sep = "\t");

# GGplot Gamma gradient
p4 <- ggplot(Res.bmr.ri.gam)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=Null.Null), colour=tab_stab_ResInv$Colours[1], linetype=tab_stab_ResInv$Linetypes[1], size=1.2)+
  geom_line(aes(x=Gamma, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], linetype=tab_stab_ResInv$Linetypes[2], size=1.2)+
  geom_line(aes(x=Gamma, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], linetype=tab_stab_ResInv$Linetypes[3], size=1.2)+
  geom_line(aes(x=Gamma, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], linetype=tab_stab_ResInv$Linetypes[5], size=1.2)+
  geom_line(aes(x=Gamma, y=Eq.Cycles), colour=tab_stab_ResInv$Colours[6], linetype=tab_stab_ResInv$Linetypes[6], size=1.2)+
  geom_line(aes(x=Gamma, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], linetype=tab_stab_ResInv$Linetypes[9], size=1.2)+
  geom_point(aes(x=Gamma, y=Null.Null), colour=tab_stab_ResInv$Colours[1], shape=tab_stab_ResInv$Symbols[1], size=2.5)+
  geom_point(aes(x=Gamma, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], shape=tab_stab_ResInv$Symbols[2], size=2.5)+
  geom_point(aes(x=Gamma, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], shape=tab_stab_ResInv$Symbols[3], size=2.5)+
  geom_point(aes(x=Gamma, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], shape=tab_stab_ResInv$Symbols[5], size=2.5)+
  geom_point(aes(x=Gamma, y=Eq.Cycles), colour=tab_stab_ResInv$Colours[6], shape=tab_stab_ResInv$Symbols[6], size=2.5)+
  geom_point(aes(x=Gamma, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], shape=tab_stab_ResInv$Symbols[9], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(1,2.5,5,10,25,50,100), labels=c(expression(10^0),"","",expression(10^1),"","",expression(10^2)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))

# GGplot Delta gradient
p5 <- ggplot(Res.bmr.ri.del)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Ddelta, y=Null.Null), colour=tab_stab_ResInv$Colours[1], linetype=tab_stab_ResInv$Linetypes[1], size=1.2)+
  geom_line(aes(x=Ddelta, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], linetype=tab_stab_ResInv$Linetypes[2], size=1.2)+
  geom_line(aes(x=Ddelta, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], linetype=tab_stab_ResInv$Linetypes[3], size=1.2)+
  geom_line(aes(x=Ddelta, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], linetype=tab_stab_ResInv$Linetypes[5], size=1.2)+
  geom_line(aes(x=Ddelta, y=Eq.Cycles), colour=tab_stab_ResInv$Colours[6], linetype=tab_stab_ResInv$Linetypes[6], size=1.2)+
  geom_line(aes(x=Ddelta, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], linetype=tab_stab_ResInv$Linetypes[9], size=1.2)+
  
  geom_point(aes(x=Ddelta, y=Null.Null), colour=tab_stab_ResInv$Colours[1], shape=tab_stab_ResInv$Symbols[1], size=2.5)+
  geom_point(aes(x=Ddelta, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], shape=tab_stab_ResInv$Symbols[2], size=2.5)+
  geom_point(aes(x=Ddelta, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], shape=tab_stab_ResInv$Symbols[3], size=2.5)+
  geom_point(aes(x=Ddelta, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], shape=tab_stab_ResInv$Symbols[5], size=2.5)+
  geom_point(aes(x=Ddelta, y=Eq.Cycles), colour=tab_stab_ResInv$Colours[6], shape=tab_stab_ResInv$Symbols[6], size=2.5)+
  geom_point(aes(x=Ddelta, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], shape=tab_stab_ResInv$Symbols[9], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.2, linetype=3)

# IGPC ====
Res.bmr.ri.gam <- read.table("./data/BMR_IGPC_RS-ResInv_Gamma.txt", h=T,  sep = "\t");
Res.bmr.ri.del <- read.table("./data/BMR_IGPC_RS-ResInv_Delta.txt", h=T,  sep = "\t");

# GGplot Gamma gradient
p6 <- ggplot(Res.bmr.ri.gam)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=Null.Null), colour=tab_stab_ResInv$Colours[1], linetype=tab_stab_ResInv$Linetypes[1], size=1.2)+
  geom_line(aes(x=Gamma, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], linetype=tab_stab_ResInv$Linetypes[2], size=1.2)+
  geom_line(aes(x=Gamma, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], linetype=tab_stab_ResInv$Linetypes[3], size=1.2)+
  geom_line(aes(x=Gamma, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], linetype=tab_stab_ResInv$Linetypes[5], size=1.2)+
  geom_line(aes(x=Gamma, y=Cycles.Eq), colour=tab_stab_ResInv$Colours[8], linetype=tab_stab_ResInv$Linetypes[8], size=1.2)+
  geom_line(aes(x=Gamma, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], linetype=tab_stab_ResInv$Linetypes[9], size=1.2)+
  geom_point(aes(x=Gamma, y=Null.Null), colour=tab_stab_ResInv$Colours[1], shape=tab_stab_ResInv$Symbols[1], size=2.5)+
  geom_point(aes(x=Gamma, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], shape=tab_stab_ResInv$Symbols[2], size=2.5)+
  geom_point(aes(x=Gamma, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], shape=tab_stab_ResInv$Symbols[3], size=2.5)+
  geom_point(aes(x=Gamma, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], shape=tab_stab_ResInv$Symbols[5], size=2.5)+
  geom_point(aes(x=Gamma, y=Cycles.Eq), colour=tab_stab_ResInv$Colours[8], shape=tab_stab_ResInv$Symbols[8], size=2.5)+
  geom_point(aes(x=Gamma, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], shape=tab_stab_ResInv$Symbols[9], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(1,2.5,5,10,25,50,100), labels=c(expression(10^0),"","",expression(10^1),"","",expression(10^2)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))

# GGplot Delta gradient
p7 <- ggplot(Res.bmr.ri.del)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Ddelta, y=Null.Null), colour=tab_stab_ResInv$Colours[1], linetype=tab_stab_ResInv$Linetypes[1], size=1.2)+
  geom_line(aes(x=Ddelta, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], linetype=tab_stab_ResInv$Linetypes[2], size=1.2)+
  geom_line(aes(x=Ddelta, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], linetype=tab_stab_ResInv$Linetypes[3], size=1.2)+
  geom_line(aes(x=Ddelta, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], linetype=tab_stab_ResInv$Linetypes[5], size=1.2)+
  geom_line(aes(x=Ddelta, y=Cycles.Eq), colour=tab_stab_ResInv$Colours[8], linetype=tab_stab_ResInv$Linetypes[8], size=1.2)+
  geom_line(aes(x=Ddelta, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], linetype=tab_stab_ResInv$Linetypes[9], size=1.2)+
  geom_point(aes(x=Ddelta, y=Null.Null), colour=tab_stab_ResInv$Colours[1], shape=tab_stab_ResInv$Symbols[1], size=2.5)+
  geom_point(aes(x=Ddelta, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], shape=tab_stab_ResInv$Symbols[2], size=2.5)+
  geom_point(aes(x=Ddelta, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], shape=tab_stab_ResInv$Symbols[3], size=2.5)+
  geom_point(aes(x=Ddelta, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], shape=tab_stab_ResInv$Symbols[5], size=2.5)+
  geom_point(aes(x=Ddelta, y=Cycles.Eq), colour=tab_stab_ResInv$Colours[8], shape=tab_stab_ResInv$Symbols[8], size=2.5)+
  geom_point(aes(x=Ddelta, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], shape=tab_stab_ResInv$Symbols[9], size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.2, linetype=3)

# Arrange & Save ====
pdf(file="./plots/Fig_3.pdf");
grid.arrange(p1,p2,p3,p5,p7 ,nrow=2, ncol=3,
             layout_matrix=rbind(c(1,2,3),c(4,5,NA)))
dev.off()
####################### FIGURE 4: Barplots Observed Regime States per Mechanisms  ====
### Loading the datasets
Inv.AC <- read.table("./data/AC_Tot_data.txt", h=T,  sep = "\t");  # AC Complete Data
Inv.EC <- read.table("./data/EC_Tot_data.txt", h=T,  sep = "\t");  # EC Complete Data
Inv.TC <- read.table("./data/TC_Tot_data.txt", h=T,  sep = "\t");  # TC Complete Data
Inv.IGPB <- read.table("./data/IGPB_Tot_data.txt", h=T,  sep = "\t");  # IGPB Complete Data
Inv.IGPC <- read.table("./data/IGPC_Tot_data.txt", h=T,  sep = "\t");  # IGPC Complete Data
AC2 <-Inv.AC[which(Inv.AC$Alpha!=1),];
EC2 <- Inv.EC[which(Inv.EC$Beta!=1),];
AC <- AC2
EC <- EC2
TC <- Inv.TC
IGPB <- Inv.IGPB
IGPC <- Inv.IGPC
State <- tab_state$State;
### Proportions of regime states for each invasion mechanisms ====
data <- rbind.data.frame(cbind(AC2$State, AC2$EQ_transit_ResInv),
                         cbind(EC2$State, EC2$EQ_transit_ResInv), 
                         cbind(TC$State, TC$EQ_transit_ResInv),
                         cbind(IGPB$State, IGPB$EQ_transit_ResInv),
                         cbind(IGPC$State, IGPC$EQ_transit_ResInv))
colnames(data) <- c("State", "Regimes")
N_Obs <- dim(data)[1];
State <- tab_state$State;

Tab_obs_State_Regime <- as.data.frame(matrix(0, ncol=4))
colnames(Tab_obs_State_Regime) <- c("Mec", "Lev", "Obs", "Prob")
for(i in 1:length(State)){
  Lev <- levels(factor(data$Regimes[data$State==State[i]]));
  Obs <- as.numeric(summary(as.factor(data$Regimes[data$State==State[i]])));
  Mec <- rep(State[i], length(Lev));
  Prob <- Proportion(Obs,sum(Obs))
  sub <- cbind.data.frame(Mec, Lev, Obs, Prob)
  Tab_obs_State_Regime <- rbind(Tab_obs_State_Regime, sub )
}
Tab_obs_State_Regime <- Tab_obs_State_Regime[-1,]

Tab_obs_State_Regime <- cbind(Tab_obs_State_Regime, rep(NA, dim(Tab_obs_State_Regime)[1]), rep(NA, dim(Tab_obs_State_Regime)[1]), rep(NA, dim(Tab_obs_State_Regime)[1]), rep(NA, dim(Tab_obs_State_Regime)[1]))
colnames(Tab_obs_State_Regime) <- c("Mechanisms", "Regimes", "N_obs", "P_obs", "Colours", "Reg", "Reg.Lev", "Div.Lev")

for(i in 1:dim(tab_stab_ResInv)[1]){
  Tab_obs_State_Regime$Colours[Tab_obs_State_Regime$Regimes==tab_stab_ResInv$Stab[i]] <- tab_stab_ResInv$Colours[i]
  Tab_obs_State_Regime$Reg[Tab_obs_State_Regime$Regimes==tab_stab_ResInv$Stab[i]] <- tab_stab_ResInv$State[i]
}

Destabilizing.S <- c("O.N", "E.O", "E.N");
Neutral.S <- c("N.N", "O.O", "E.E");
Stabilizing.S <- c("N.E", "N.O", "O.E");
Div.loss <- c("Vulnerability");
Div.neut <- c("Resistance", "Substitution");
Div.gain <- c("Integration", "Rescue", "Occupancy");

Tab_obs_State_Regime$Div.Lev[Tab_obs_State_Regime$Mechanisms=="Vulnerability"] <- LETTERS[1];
#Div.neut <- c("Substitution");
#Tab_obs_State_Regime$Reg.Lev[Tab_obs_State_Regime$Reg==Neutral.S[1]] <- LETTERS[2];

for(i in 1:3){
  if(i<3){ Tab_obs_State_Regime$Div.Lev[Tab_obs_State_Regime$Mechanisms==Div.neut[i]] <- LETTERS[2];}
  Tab_obs_State_Regime$Reg.Lev[Tab_obs_State_Regime$Reg==Destabilizing.S[i]] <- LETTERS[1];
  Tab_obs_State_Regime$Reg.Lev[Tab_obs_State_Regime$Reg==Neutral.S[i]] <- LETTERS[2];
  Tab_obs_State_Regime$Reg.Lev[Tab_obs_State_Regime$Reg==Stabilizing.S[i]] <- LETTERS[3];
  Tab_obs_State_Regime$Div.Lev[Tab_obs_State_Regime$Mechanisms==Div.gain[i]] <- LETTERS[3];
}

#### Barplot ====
p <- ggplot(data=Tab_obs_State_Regime, aes(x=Mechanisms, y=P_obs, color=Reg))+
  geom_bar(stat="identity", color="black", fill=Tab_obs_State_Regime$Colours)+
  labs(x="Invasion mechanisms", y= expression("Regime states"~S['RES']*.*S['INV']~"(%)"))+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15), legend.position="none")

# Save ====
pdf(file="./plots/Fig_4-v0.pdf", width=8, height=7);
p;
dev.off()
####################### FIGURE S1: Structure of resident communities along gradients of abiotic condition and size structure ====
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

# Arrange & Save ====
pdf(file="./plots/Fig_S1.pdf");
grid.arrange(p1, p2, ncol=2, nrow=1, newpage = T)
dev.off()
####################### FIGURE S2 & S3: Barplots of Invasion mechanisms & Regime states 1) Overall observation; 2) across food web modules ====
# with (S2) or without (S3) Community's resistance
### Loading the datasets
Inv.AC <- read.table("./data/AC_Tot_data.txt", h=T,  sep = "\t");  # AC Complete Data
Inv.EC <- read.table("./data/EC_Tot_data.txt", h=T,  sep = "\t");  # EC Complete Data
Inv.TC <- read.table("./data/TC_Tot_data.txt", h=T,  sep = "\t");  # TC Complete Data
Inv.IGPB <- read.table("./data/IGPB_Tot_data.txt", h=T,  sep = "\t");  # IGPB Complete Data
Inv.IGPC <- read.table("./data/IGPC_Tot_data.txt", h=T,  sep = "\t");  # IGPC Complete Data
# 2 version of selection that either (S2) include Resistance (invasion failure) or (S3) exclude it from the observations
# Version S2: Including Resistance in the observations ====
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
# Version S3: Excluding Resistance from the observations ====
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

### Invasion mechanisms ====
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
pdf(file="./plots/Fig_S2.pdf", width=16, height=12);
pdf(file="./plots/Fig_S3.pdf", width=16, height=12);
grid.arrange(p,q,r,s, nrow=2, ncol=2, layout_matrix=rbind(c(1,2),c(3,4)))
dev.off()

####################### FIGURE S4: Invasion mechanisms & Regime states in AC & EC along competing species BMR ====
ratio <- c(0.1, 0.2, 0.5, 2, 5, 10);
# AC ====
Inv <- read.table("./data/AC_Tot_data.txt", h=T,  sep = "\t");  # AC Complete Data
AC <- Inv[which(Inv$Alpha!=1 & Inv$Beta==10),];
State_code <- rep(0, dim(AC)[1]);
AC <- cbind.data.frame(AC, State_code);

for(j in 1:dim(tab_state)[1]){
  if(length(AC$State_code[which(AC$State == tab_state$State[j])]) != 0){
    AC$State_code[which(AC$State == tab_state$State[j]) ] <- tab_state$Code[j];
  }else{next}
}

p <- list(); q <- list();
for(i in 1:length(ratio)){
  sub <- AC[which(AC$Alpha==ratio[i]),];
  
  ylabel <- ifelse(i==1, expression("Temperature ("*degree*C*')'), "")
  p[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(State_code)))+theme_bw()+
    geom_tile()+
    labs(x="", y= ylabel)+
    scale_fill_manual(values = tab_state$Colour, breaks=tab_state$Code,
                      labels=tab_state$State, name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
  
  q[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(Zone_EQtransit_ResInv)))+theme_bw()+
    geom_tile()+
    labs(x="", y= ylabel)+
    scale_fill_manual(values = tab_stab_ResInv$Colour, breaks=tab_stab_ResInv$Code,
                      labels=tab_stab_ResInv$Stab, name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
}
# EC ====
Inv <- read.table("./data/EC_Tot_data.txt", h=T,  sep = "\t");  # EC Complete Data
EC <- Inv[which(Inv$Beta!=1 & Inv$Alpha==10),];
State_code <- rep(0, dim(EC)[1]);
EC$Beta <- 1/EC$Beta;
EC <- cbind.data.frame(EC, State_code);

for(j in 1:dim(tab_state)[1]){
  if(length(EC$State_code[which(EC$State == tab_state$State[j])]) != 0){
    EC$State_code[which(EC$State == tab_state$State[j]) ] <- tab_state$Code[j];
  }else{next}
}

r <- list(); s <- list();
for(i in 1:length(ratio)){
  sub <- EC[which(EC$Beta==ratio[i]),];
  
  ylabel <- ifelse(i==1, expression("Temperature ("*degree*C*')'), "")
  r[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(State_code)))+theme_bw()+
    geom_tile()+
    labs(x="", y= ylabel)+
    scale_fill_manual(values = tab_state$Colour, breaks=tab_state$Code,
                      labels=tab_state$State, name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
  
  s[[i]] <- ggplot(sub, aes(x=Car, y=Temp, fill=factor(Zone_EQtransit_ResInv)))+theme_bw()+
    geom_tile()+
    labs(x=expression("Nutrient levels ("*g.L^-1*')', ""), y= ylabel)+
    scale_fill_manual(values = tab_stab_ResInv$Colour, breaks=tab_stab_ResInv$Code,
                      labels=tab_stab_ResInv$Stab, name=NULL, guide="none")+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=15))
}
# Arrange & Save ====
pdf(file="./plots/Fig_S4.pdf", width=20, height=11);
grid.arrange(grobs=c(p,q, r, s), ncol=6, nrow=4)
dev.off()
####################### FIGURE S5: P* & R* rules in AC, EC ====
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
  
  ylabel <- ifelse(i==1, expression("Temperature ("*degree*C*')'), "")
  
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
pdf(file="./plots/Fig_S5.pdf", width=20, height=11);
grid.arrange(grobs=c(p,q, r, s), ncol=6, nrow=4)
dev.off()

####################### FIGURE S6: R* rule & Initial biomass growth of IG-prey in IGP ====
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
pdf(file="./plots/Fig_S6.pdf", width=20, height=11);
grid.arrange(grobs=c(s,v,t,w,u,x), ncol=6, nrow=3)
dev.off()

####################### FIGURE S7: Delta D (Biodiversity change) & Delta S (change in stability regime) along BMR gradients ====
# AC ====
Res.bmr.is <- read.table("./data/BMR_AC_IS.txt", h=T,  sep = "\t");
Res.bmr.is <- Res.bmr.is[-5,]; Res.bmr.is <- Res.bmr.is[,-1];
Res.bmr.ri <- read.table("./data/BMR_AC_RS-ResInv.txt", h=T,  sep = "\t");
Res.bmr.ri <- Res.bmr.ri[-5,];

p1 <- ggplot(Res.bmr.is)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Alpha, y=Bio_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Alpha, y=Bio_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Alpha, y=Bio_loss), size=1.2)+
  geom_point(aes(x=Alpha, y=Bio_gain), shape=15, size=2.5)+
  geom_point(aes(x=Alpha, y=Bio_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Alpha, y=Bio_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3)

p2 <- ggplot(Res.bmr.ri)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Alpha, y=DeltaS_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Alpha, y=DeltaS_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Alpha, y=DeltaS_loss), size=1.2)+
  geom_point(aes(x=Alpha, y=DeltaS_gain), shape=15, size=2.5)+
  geom_point(aes(x=Alpha, y=DeltaS_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Alpha, y=DeltaS_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3)

# EC ====
Res.bmr.is <- read.table("./data/BMR_EC_IS.txt", h=T,  sep = "\t");
Res.bmr.is <- Res.bmr.is[-5,]; Res.bmr.is <- Res.bmr.is[,-1]; Res.bmr.is$Beta <- (1/Res.bmr.is$Beta);
Res.bmr.ri <- read.table("./data/BMR_EC_RS-ResInv.txt", h=T,  sep = "\t");
Res.bmr.ri <- Res.bmr.ri[-5,]; Res.bmr.ri$Beta <- (1/Res.bmr.ri$Beta);

p3 <- ggplot(Res.bmr.is)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Beta, y=Bio_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Beta, y=Bio_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Beta, y=Bio_loss), size=1.2)+
  geom_point(aes(x=Beta, y=Bio_gain), shape=15, size=2.5)+
  geom_point(aes(x=Beta, y=Bio_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Beta, y=Bio_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3)

p4 <- ggplot(Res.bmr.ri)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Beta, y=DeltaS_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Beta, y=DeltaS_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Beta, y=DeltaS_loss), size=1.2)+
  geom_point(aes(x=Beta, y=DeltaS_gain), shape=15, size=2.5)+
  geom_point(aes(x=Beta, y=DeltaS_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Beta, y=DeltaS_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.2, linetype=3);

# IGP ====
Res.bmr.is <- read.table("./data/BMR_IGP_IS_Beta.txt", h=T,  sep = "\t");
Res.bmr.ri <- read.table("./data/BMR_IGP_RS-ResInv_Beta.txt", h=T,  sep = "\t");
Res.bmr.ri <- Res.bmr.ri[-c(8:14),]

p5 <- ggplot(Res.bmr.is)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Beta, y=Bio_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Beta, y=Bio_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Beta, y=Bio_loss), size=1.2)+
  geom_point(aes(x=Beta, y=Bio_gain), shape=15, size=2.5)+
  geom_point(aes(x=Beta, y=Bio_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Beta, y=Bio_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3)

p6 <- ggplot(Res.bmr.ri)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Beta, y=DeltaS_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Beta, y=DeltaS_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Beta, y=DeltaS_loss), size=1.2)+
  geom_point(aes(x=Beta, y=DeltaS_gain), shape=15, size=2.5)+
  geom_point(aes(x=Beta, y=DeltaS_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Beta, y=DeltaS_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.2, linetype=3)

# TC ==== 
Res.bmr.is.gam <- read.table("./data/BMR_TC_IS_Gamma.txt", h=T,  sep = "\t");
Res.bmr.is.del <- read.table("./data/BMR_TC_IS_Delta.txt", h=T,  sep = "\t");
Res.bmr.ri.gam <- read.table("./data/BMR_TC_RS-ResInv_Gamma.txt", h=T,  sep = "\t");
Res.bmr.ri.del <- read.table("./data/BMR_TC_RS-ResInv_Delta.txt", h=T,  sep = "\t");
# Diversity change
p7 <- ggplot(Res.bmr.is.gam)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=Bio_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Gamma, y=Bio_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Gamma, y=Bio_loss), size=1.2)+
  geom_point(aes(x=Gamma, y=Bio_gain), shape=15, size=2.5)+
  geom_point(aes(x=Gamma, y=Bio_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Gamma, y=Bio_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(1,2.5,5,10,25,50,100), labels=c(expression(10^0),"","",expression(10^1),"","",expression(10^2)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20));

p8 <- ggplot(Res.bmr.is.del)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Ddelta, y=Bio_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Ddelta, y=Bio_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Ddelta, y=Bio_loss), size=1.2)+
  geom_point(aes(x=Ddelta, y=Bio_gain), shape=15, size=2.5)+
  geom_point(aes(x=Ddelta, y=Bio_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Ddelta, y=Bio_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.0, linetype=3);
# Regime state change
p9 <- ggplot(Res.bmr.ri.gam)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=DeltaS_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Gamma, y=DeltaS_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Gamma, y=DeltaS_loss), size=1.2)+
  geom_point(aes(x=Gamma, y=DeltaS_gain), shape=15, size=2.5)+
  geom_point(aes(x=Gamma, y=DeltaS_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Gamma, y=DeltaS_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(1,2.5,5,10,25,50,100), labels=c(expression(10^0),"","",expression(10^1),"","",expression(10^2)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20));

p10 <- ggplot(Res.bmr.ri.del)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Ddelta, y=DeltaS_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Ddelta, y=DeltaS_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Ddelta, y=DeltaS_loss), size=1.2)+
  geom_point(aes(x=Ddelta, y=DeltaS_gain), shape=15, size=2.5)+
  geom_point(aes(x=Ddelta, y=DeltaS_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Ddelta, y=DeltaS_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.2, linetype=3);

# IGPC ====
Res.bmr.is.gam <- read.table("./data/BMR_IGPC_IS_Gamma.txt", h=T,  sep = "\t");
Res.bmr.is.del <- read.table("./data/BMR_IGPC_IS_Delta.txt", h=T,  sep = "\t");
Res.bmr.ri.gam <- read.table("./data/BMR_IGPC_RS-ResInv_Gamma.txt", h=T,  sep = "\t");
Res.bmr.ri.del <- read.table("./data/BMR_IGPC_RS-ResInv_Delta.txt", h=T,  sep = "\t");

# Diversity change
p11 <- ggplot(Res.bmr.is.gam)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=Bio_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Gamma, y=Bio_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Gamma, y=Bio_loss), size=1.2)+
  geom_point(aes(x=Gamma, y=Bio_gain), shape=15, size=2.5)+
  geom_point(aes(x=Gamma, y=Bio_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Gamma, y=Bio_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(1,2.5,5,10,25,50,100), labels=c(expression(10^0),"","",expression(10^1),"","",expression(10^2)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20));

p12 <- ggplot(Res.bmr.is.del)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Ddelta, y=Bio_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Ddelta, y=Bio_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Ddelta, y=Bio_loss), size=1.2)+
  geom_point(aes(x=Ddelta, y=Bio_gain), shape=15, size=2.5)+
  geom_point(aes(x=Ddelta, y=Bio_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Ddelta, y=Bio_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.2, linetype=3);

# Regime state change
p13 <- ggplot(Res.bmr.ri.gam)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=DeltaS_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Gamma, y=DeltaS_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Gamma, y=DeltaS_loss), size=1.2)+
  geom_point(aes(x=Gamma, y=DeltaS_gain), shape=15, size=2.5)+
  geom_point(aes(x=Gamma, y=DeltaS_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Gamma, y=DeltaS_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(1,2.5,5,10,25,50,100), labels=c(expression(10^0),"","",expression(10^1),"","",expression(10^2)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20));

p14 <- ggplot(Res.bmr.ri.del)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Ddelta, y=DeltaS_gain), linetype=2, size=1.2)+
  geom_line(aes(x=Ddelta, y=DeltaS_neutral), linetype=3, size=1.2)+
  geom_line(aes(x=Ddelta, y=DeltaS_loss), size=1.2)+
  geom_point(aes(x=Ddelta, y=DeltaS_gain), shape=15, size=2.5)+
  geom_point(aes(x=Ddelta, y=DeltaS_neutral), shape=16, size=2.5)+
  geom_point(aes(x=Ddelta, y=DeltaS_loss), shape=17, size=2.5)+
  labs(x="", y="")+
  scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10), labels=c(expression(10^-1),"","",1,"","",expression(10^1)))+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  geom_vline(xintercept = 1, size=1.2, linetype=3);

# Arrange & Save ====
pdf(file="./plots/Fig_S7.pdf", width=20, height=10);
grid.arrange(p1,p3,p5,p8,p12,
             p2,p4,p6,p10,p14 ,nrow=2, ncol=5)
dev.off()
####################### FIGURE S8: Invasion mechanisms & Regime states in TC & IGPP along BMR gradients ====
# TC ====
setwd(dir1);
Tc.is.bet <- read.table("./data/BMR_TC_IS_Beta.txt", h=T,  sep = "\t");
Tc.ri.bet <- read.table("./data/BMR_TC_RS-ResInv_Beta.txt", h=T,  sep = "\t");
Tc.is.gam <- read.table("./data/BMR_TC_IS_Gamma.txt", h=T,  sep = "\t");
Tc.ri.gam <- read.table("./data/BMR_TC_RS-ResInv_Gamma.txt", h=T,  sep = "\t");

# Beta gradient
p.bio <- ggplot(Tc.is.bet)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Beta, y=Integration), colour=tab_state$Colours[6], size=1.2)+
  geom_line(aes(x=Beta, y=Rescue), colour=tab_state$Colours[4], size=1.2)+
  geom_line(aes(x=Beta, y=Resistance), colour=tab_state$Colours[2], size=1.2)+
  geom_point(aes(x=Beta, y=Integration), colour=tab_state$Colours[6], shape=tab_state$Symbols[6], size=2.5)+
  geom_point(aes(x=Beta, y=Rescue), colour=tab_state$Colours[4], shape=tab_state$Symbols[4], size=2.5)+
  geom_point(aes(x=Beta, y=Resistance), colour=tab_state$Colours[2], shape=tab_state$Symbols[2], size=2.5)+
  labs(x="", y="")+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  scale_x_log10(breaks = c(1,2.5,5,10), labels=c(expression(10^0),"","",expression(10^1)))

p.ri <- ggplot(Tc.ri.bet)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Beta, y=Null.Null), colour=tab_stab_ResInv$Colours[1], linetype=3, size=1.2)+
  geom_line(aes(x=Beta, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], size=1.2)+
  geom_line(aes(x=Beta, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], size=1.2)+
  geom_line(aes(x=Beta, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], linetype=3, size=1.2)+
  geom_line(aes(x=Beta, y=Eq.Cycles), colour=tab_stab_ResInv$Colours[6], size=1.2)+
  geom_line(aes(x=Beta, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], linetype=3, size=1.2)+
  geom_point(aes(x=Beta, y=Null.Null), colour=tab_stab_ResInv$Colours[1], shape=tab_stab_ResInv$Symbols[1], size=2.5)+
  geom_point(aes(x=Beta, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], shape=tab_stab_ResInv$Symbols[2], size=2.5)+
  geom_point(aes(x=Beta, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], shape=tab_stab_ResInv$Symbols[3], size=2.5)+
  geom_point(aes(x=Beta, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], shape=tab_stab_ResInv$Symbols[5], size=2.5)+
  geom_point(aes(x=Beta, y=Eq.Cycles), colour=tab_stab_ResInv$Colours[6], shape=tab_stab_ResInv$Symbols[6], size=2.5)+
  geom_point(aes(x=Beta, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], shape=tab_stab_ResInv$Symbols[9], size=2.5)+
  labs(x="", y="")+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  scale_x_log10(breaks = c(1,2.5,5,10), labels=c(expression(10^0),"","",expression(10^1)))

# Gamma gradient
q.bio <- ggplot(Tc.is.gam)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=Resistance), colour=tab_state$Colours[2], size=1.2)+
  geom_line(aes(x=Gamma, y=Rescue), colour=tab_state$Colours[4], size=1.2)+
  geom_line(aes(x=Gamma, y=Integration), colour=tab_state$Colours[6], size=1.2)+
  geom_point(aes(x=Gamma, y=Resistance), colour=tab_state$Colours[2], shape=tab_state$Symbols[2], size=2.5)+
  geom_point(aes(x=Gamma, y=Rescue), colour=tab_state$Colours[4], shape=tab_state$Symbols[4], size=2.5)+
  geom_point(aes(x=Gamma, y=Integration), colour=tab_state$Colours[6], shape=tab_state$Symbols[6], size=2.5)+
  labs(x="", y="")+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  scale_x_log10(breaks = c(1,2.5,5,10,25,50,100), labels=c(expression(10^0),"","",expression(10^1),"","",expression(10^2)))

q.ri <- ggplot(Tc.ri.gam)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=Null.Null), colour=tab_stab_ResInv$Colours[1], linetype=3, size=1.2)+
  geom_line(aes(x=Gamma, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], size=1.2)+
  geom_line(aes(x=Gamma, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], size=1.2)+
  geom_line(aes(x=Gamma, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], linetype=3, size=1.2)+
  geom_line(aes(x=Gamma, y=Eq.Cycles), colour=tab_stab_ResInv$Colours[6], size=1.2)+
  geom_line(aes(x=Gamma, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], linetype=3, size=1.2)+
  geom_point(aes(x=Gamma, y=Null.Null), colour=tab_stab_ResInv$Colours[1], shape=tab_stab_ResInv$Symbols[1], size=2.5)+
  geom_point(aes(x=Gamma, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], shape=tab_stab_ResInv$Symbols[2], size=2.5)+
  geom_point(aes(x=Gamma, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], shape=tab_stab_ResInv$Symbols[3], size=2.5)+
  geom_point(aes(x=Gamma, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], shape=tab_stab_ResInv$Symbols[5], size=2.5)+
  geom_point(aes(x=Gamma, y=Eq.Cycles), colour=tab_stab_ResInv$Colours[6], shape=tab_stab_ResInv$Symbols[6], size=2.5)+
  geom_point(aes(x=Gamma, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], shape=tab_stab_ResInv$Symbols[9], size=2.5)+
  labs(x="", y="")+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  scale_x_log10(breaks = c(1,2.5,5,10,25,50,100), labels=c(expression(10^0),"","",expression(10^1),"","",expression(10^2)))

# IGPP ====
Igp.is.bet <- read.table("./data/BMR_IGPC_IS_Beta.txt", h=T,  sep = "\t");
Igp.ri.bet <- read.table("./data/BMR_IGPC_RS-ResInv_Beta.txt", h=T,  sep = "\t");
Igp.is.gam <- read.table("./data/BMR_IGPC_IS_Gamma.txt", h=T,  sep = "\t");
Igp.ri.gam <- read.table("./data/BMR_IGPC_RS-ResInv_Gamma.txt", h=T,  sep = "\t");

# Beta gradient
r.bio <- ggplot(Igp.is.bet)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Beta, y=Resistance), colour=tab_state$Colours[2], size=1.2)+
  geom_line(aes(x=Beta, y=Substitution), colour=tab_state$Colours[3],  size=1.2)+
  geom_line(aes(x=Beta, y=Occupancy), colour=tab_state$Colours[5], size=1.2)+
  geom_line(aes(x=Beta, y=Integration), colour=tab_state$Colours[6], size=1.2)+
  geom_point(aes(x=Beta, y=Resistance), colour=tab_state$Colours[2], shape=tab_state$Symbols[2], size=2.5)+
  geom_point(aes(x=Beta, y=Substitution), colour=tab_state$Colours[3], shape=tab_state$Symbols[3], size=2.5)+
  geom_point(aes(x=Beta, y=Occupancy), colour=tab_state$Colours[5], shape=tab_state$Symbols[5], size=2.5)+
  geom_point(aes(x=Beta, y=Integration), colour=tab_state$Colours[6], shape=tab_state$Symbols[6], size=2.5)+
  labs(x="", y="")+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  scale_x_log10(breaks = c(1,2.5,5,10), labels=c(expression(10^0),"","",expression(10^1)))

r.ri <- ggplot(Igp.ri.bet)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Beta, y=Null.Null), colour=tab_stab_ResInv$Colours[1], linetype=3, size=1.2)+
  geom_line(aes(x=Beta, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], size=1.2)+
  geom_line(aes(x=Beta, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], size=1.2)+
  geom_line(aes(x=Beta, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], linetype=3, size=1.2)+
  geom_line(aes(x=Beta, y=Cycles.Eq), colour=tab_stab_ResInv$Colours[8], size=1.2)+
  geom_line(aes(x=Beta, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], linetype=3, size=1.2)+
  geom_point(aes(x=Beta, y=Null.Null), colour=tab_stab_ResInv$Colours[1], shape=tab_stab_ResInv$Symbols[1], size=2.5)+
  geom_point(aes(x=Beta, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], shape=tab_stab_ResInv$Symbols[2], size=2.5)+
  geom_point(aes(x=Beta, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], shape=tab_stab_ResInv$Symbols[3], size=2.5)+
  geom_point(aes(x=Beta, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], shape=tab_stab_ResInv$Symbols[5], size=2.5)+
  geom_point(aes(x=Beta, y=Cycles.Eq), colour=tab_stab_ResInv$Colours[8], shape=tab_stab_ResInv$Symbols[8], size=2.5)+
  geom_point(aes(x=Beta, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], shape=tab_stab_ResInv$Symbols[9], size=2.5)+
  labs(x="", y="")+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  scale_x_log10(breaks = c(1,2.5,5,10), labels=c(expression(10^0),"","",expression(10^1)))

# Gamma gradient
s.bio <- ggplot(Igp.is.gam)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=Resistance), colour=tab_state$Colours[2], size=1.2)+
  geom_line(aes(x=Gamma, y=Substitution), colour=tab_state$Colours[3],  size=1.2)+
  geom_line(aes(x=Gamma, y=Occupancy), colour=tab_state$Colours[5], size=1.2)+
  geom_line(aes(x=Gamma, y=Integration), colour=tab_state$Colours[6], size=1.2)+
  geom_point(aes(x=Gamma, y=Resistance), colour=tab_state$Colours[2], shape=tab_state$Symbols[2], size=2.5)+
  geom_point(aes(x=Gamma, y=Substitution), colour=tab_state$Colours[3], shape=tab_state$Symbols[3], size=2.5)+
  geom_point(aes(x=Gamma, y=Occupancy), colour=tab_state$Colours[5], shape=tab_state$Symbols[5], size=2.5)+
  geom_point(aes(x=Gamma, y=Integration), colour=tab_state$Colours[6], shape=tab_state$Symbols[6], size=2.5)+
  labs(x="", y="")+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  scale_x_log10(breaks = c(1,2.5,5,10,25,50,100), labels=c(expression(10^0),"","",expression(10^1),"","",expression(10^2)))

s.ri <- ggplot(Igp.ri.gam)+ylim(0,100)+theme_bw()+
  geom_line(aes(x=Gamma, y=Null.Null), colour=tab_stab_ResInv$Colours[1], linetype=3, size=1.2)+
  geom_line(aes(x=Gamma, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], size=1.2)+
  geom_line(aes(x=Gamma, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], size=1.2)+
  geom_line(aes(x=Gamma, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], linetype=3, size=1.2)+
  geom_line(aes(x=Gamma, y=Cycles.Eq), colour=tab_stab_ResInv$Colours[8], size=1.2)+
  geom_line(aes(x=Gamma, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], linetype=3, size=1.2)+
  geom_point(aes(x=Gamma, y=Null.Null), colour=tab_stab_ResInv$Colours[1], shape=tab_stab_ResInv$Symbols[1], size=2.5)+
  geom_point(aes(x=Gamma, y=Null.Eq), colour=tab_stab_ResInv$Colours[2], shape=tab_stab_ResInv$Symbols[2], size=2.5)+
  geom_point(aes(x=Gamma, y=Null.Cycles), colour=tab_stab_ResInv$Colours[3], shape=tab_stab_ResInv$Symbols[3], size=2.5)+
  geom_point(aes(x=Gamma, y=Eq.Eq), colour=tab_stab_ResInv$Colours[5], shape=tab_stab_ResInv$Symbols[5], size=2.5)+
  geom_point(aes(x=Gamma, y=Cycles.Eq), colour=tab_stab_ResInv$Colours[8], shape=tab_stab_ResInv$Symbols[8], size=2.5)+
  geom_point(aes(x=Gamma, y=Cycles.Cycles), colour=tab_stab_ResInv$Colours[9], shape=tab_stab_ResInv$Symbols[9], size=2.5)+
  labs(x="", y="")+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  scale_x_log10(breaks = c(1,2.5,5,10,25,50,100), labels=c(expression(10^0),"","",expression(10^1),"","",expression(10^2)))

# Arrange & Save ====
pdf(file="./plots/Fig_S8.pdf", width=18, height=10);
grid.arrange(p.bio, q.bio, r.bio, s.bio,
             p.ri, q.ri, r.ri, s.ri, ncol=4, nrow=2)
dev.off()

####################### FIGURE S9: Barplots of Regime states ~ Invading mechanisms ====
# Require to load the table obtained for Fig 4
# and removing Resistance from the outcome:
Tab_obs_State_Regime <- Tab_obs_State_Regime[which(Tab_obs_State_Regime$Mechanisms!="Resistance"),]

## Delta S regarding to Invasion Mechanisms
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

### Barplots ====
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
pdf(file="./plots/Fig_S9.pdf", width=15, height=9);

grid.arrange(q, r, ncol=2, nrow=1)

dev.off()

####################### FIGURE S10-S11: Population biomass extremes along abiotic gradients ====
# Creating the subsets from the original dataframes: selection for gradients of N=5 and T=10
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

######### FIGURE S10 Population biomass extremes along temperature gradient for a fixed N=5 (selection for Year=10) ====
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
### Version initial but not selected: using geom_points ====
#AC
Yax <- expression(R[RES]);
Bax <- expression(R[INV]);
Cax <- expression(C[RES]);

p[[1]] <- ggplot(N_AC)+ylim(0,1)+theme_bw()+
  geom_point(aes(x=Temperature, y=A_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=A_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x=expression("Temperature ("*degree*C*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[1]] <- ggplot(N_AC)+ylim(0,.75)+theme_bw()+
  geom_point(aes(x=Temperature, y=B_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=B_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

r[[1]] <- ggplot(N_AC)+ylim(0,.75)+theme_bw()+
  geom_point(aes(x=Temperature, y=C_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=C_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

#EC
Yax <- expression(R[RES])
Bax <- expression(C[INV])
Cax <- expression(C[RES])

p[[2]] <- ggplot(N_EC)+ylim(0,1)+theme_bw()+
  geom_point(aes(x=Temperature, y=A_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=A_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x=expression("Temperature ("*degree*C*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[2]] <- ggplot(N_EC)+ylim(0,.75)+theme_bw()+
  geom_point(aes(x=Temperature, y=B_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=B_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

r[[2]] <- ggplot(N_EC)+ylim(0,.75)+theme_bw()+
  geom_point(aes(x=Temperature, y=C_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=C_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

#TC
Yax <- expression(R[RES])
Bax <- expression(C[RES])
Cax <- expression(P[INV])

p[[3]] <- ggplot(N_TC)+ylim(0,5)+theme_bw()+
  geom_point(aes(x=Temperature, y=A_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=A_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x=expression("Temperature ("*degree*C*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[3]] <- ggplot(N_TC)+ylim(0,.75)+theme_bw()+
  geom_point(aes(x=Temperature, y=B_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=B_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

r[[3]] <- ggplot(N_TC)+ylim(0,3)+theme_bw()+
  geom_point(aes(x=Temperature, y=C_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=C_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

#IGPB
Yax <- expression(R[RES])
Bax <- expression(C[INV])
Cax <- expression(P[RES])

p[[4]] <- ggplot(N_IGPB)+ylim(0,1)+theme_bw()+
  geom_point(aes(x=Temperature, y=A_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=A_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x=expression("Temperature ("*degree*C*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[4]] <- ggplot(N_IGPB)+ylim(0,.75)+theme_bw()+
  geom_point(aes(x=Temperature, y=B_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=B_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

r[[4]] <- ggplot(N_IGPB)+ylim(0,.75)+theme_bw()+
  geom_point(aes(x=Temperature, y=C_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=C_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

#IGPC
Yax <- expression(R[RES])
Bax <- expression(C[RES])
Cax <- expression(P[INV])

p[[5]] <- ggplot(N_IGPC)+ylim(0,1)+theme_bw()+
  geom_point(aes(x=Temperature, y=A_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=A_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x=expression("Temperature ("*degree*C*')'), y=Yax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

q[[5]] <- ggplot(N_IGPC)+ylim(0,.75)+theme_bw()+
  geom_point(aes(x=Temperature, y=B_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=B_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x="", y=Bax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

r[[5]] <- ggplot(N_IGPC)+ylim(0,.75)+theme_bw()+
  geom_point(aes(x=Temperature, y=C_Min, colour=Gamma), size=2.5)+
  geom_point(aes(x=Temperature, y=C_Max, colour=Gamma), size=2.5)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])+
  labs(x="", y=Cax)+
  theme(axis.text=element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none");

theme(axis.text=element_text(size=20),
      axis.title = element_text(size=20),
      legend.position = "right",
      legend.justification = "right",
      legend.background = element_rect(fill = "white", color = "black"),
      legend.text = element_text(size=15),
      legend.title = element_text(size=15));
a<- ggplot(N_IGPC)+ylim(0,.75)+theme_bw()+
  aes(x=Temperature, y=C_Min, colour=Gamma)+
  aes(x=Temperature, y=C_Max, colour=Gamma)+
  geom_point(size=1.5)+
  labs(x="", y=Cax)+
  scale_color_manual(breaks=c("100", "25", "4", "1"), values=col.bmr[4:1])

a <- a+theme(legend.position = "right",
             legend.justification = "right",
             legend.background = element_rect(fill = "white", color = "black"))
#Arrange & Save
pdf(file="./Plots/Fig_S10.pdf");
grid.arrange(grobs=c(r,q,p), ncol=5);
dev.off()

### Version Selected: Replacing the points by lines ====
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
pdf(file="./Plots/Fig_S10.pdf", width=25, height=15);
grid.arrange(grobs=c(r,q,p), ncol=5);
dev.off()
######### FIGURE S11 Population biomass extremes along nutrient gradient for a fixed T=10 (selection for Year=10) ====
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
pdf(file="./Plots/Fig_S11.pdf", width=25, height=15);
grid.arrange(grobs=c(r,q,p), ncol=5)
dev.off()