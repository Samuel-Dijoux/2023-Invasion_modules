########################################################################
### R Script                      BMR_MeanProp.R
###   To calculate the averaged percentages of
###       _a. changes in community stability in each module after 5000 years (resident and invaded systems)
###       _b. changes in biodiversity, community composition & community stability transition between resident and invaded systems
########################################################################
library(lattice)

source("./Functions.R");

# Species Body mass ratios across modules
Alpha <- c(1, 2, 5, 10);
Alpha2<- c(0.1, 0.2, 0.5, 0.75, 1, 1.5, 2, 5, 10); 
Beta  <- c(1, 2, 5, 10);
Beta2 <- c(0.1, 0.2, 0.75, 1, 1.5, 2, 5, 10);

Gamma <- c(1, 1.5, 2, 2.5, 3, 3.75, 4, 5, 7.5, 10, 15, 20, 25, 50, 100);
Delta <- c(0.1, 0.2, 0.4, 0.5, 1, 2, 2.5, 5, 10);
Delta2 <- c(0.01, 0.02, 0.04, 0.05, 0.1, 0.2, 0.25, 0.4, 0.5, 1, 2, 2.5, 5, 10);
Delta3 <- c(0.1, 0.2, 0.25, 0.4, 0.5, 1, 2, 2.5, 4, 5, 10, 20, 25, 50, 100);
# 
Fw <- c('AC', 'EC', 'IGPB', 'IGPC', 'TC');
BioD <- c(3, 2, 1, 0, -1, -2);
Zone.i <- c(1, 2, 3); Zone.f <- c(0, 1, 2, 3);
State <- c("Integration", "Occupancy", "Rescue", "Substitution", "Resistance", "Vulnerability");

############# 2spAB; C-R system without invasion ====
Res <- read.table("./data/CR_Tot_Alpha.txt", h=T,  sep = "\t");
Alpha <- c(1, 2, 5, 10);# BMR Cr:Rr

#### % Stability (initial vs after 5000 years) ====
Stab <- as.vector(names(summary(as.factor(Res$EQ_transit_InitFin)) ));
prop <- as.data.frame(matrix(0, nrow= length(Alpha), ncol=length(Stab)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Alpha)){
  sub <- Res[which(Res$Alpha == Alpha[i]),];
  Stat_size <- c();
  for(j in 1:length(Stab)){
    Stat_size <- append(Stat_size, dim(sub[which(sub$EQ_transit_InitFin == Stab[j]),])[1]);
  }
  prop_vec <- Proportion(Stat_size,tot);
  prop[i,] <- prop_vec;
}
colnames(prop) <- Stab;
prop <- cbind.data.frame(prop, Alpha);

############# 2spAC; C-R system without invasion ====
Res <- read.table("./data/CR_Tot_Gamma.txt", h=T,  sep = "\t");
Gamma <- c(1, 1.5, 2, 2.5, 3, 3.75, 4, 5, 7.5, 10, 15, 20, 25, 50, 100);
#### % Stability (at 5000 years) ====
Stab <- as.vector(names(summary(as.factor(Res$Dyn_state_f)) ));
prop <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(Stab)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Gamma)){
  sub <- Res[which(Res$Gamma == Gamma[i]),];
  Stat_size <- c();
  for(j in 1:length(Stab)){
    Stat_size <- append(Stat_size, dim(sub[which(sub$Dyn_state_f == Stab[j]),])[1]);
  }
  prop_vec <- Proportion(Stat_size,tot);
  prop[i,] <- prop_vec;
}

colnames(prop) <- Stab;
Cycle <- prop$Cr.Rr.cycles;
Eq <- (prop$Cr.Rr.eq+prop$R.eq)
prop <- cbind.data.frame(prop[,c(1,2,4,3)], Eq, Cycle, Gamma);

write.table(prop, file="./data/BMR_2spAC_RS-F.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
#### % Stability (initial vs after 5000 years) ====
Stab <- as.vector(names(summary(as.factor(Res$EQ_transit_InitFin)) ));
prop <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(Stab)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Gamma)){
  sub <- Res[which(Res$Gamma == Gamma[i]),];
  Stat_size <- c();
  for(j in 1:length(Stab)){
    Stat_size <- append(Stat_size, dim(sub[which(sub$EQ_transit_InitFin == Stab[j]),])[1]);
  }
  prop_vec <- Proportion(Stat_size,tot);
  prop[i,] <- prop_vec;
}
colnames(prop) <- Stab;
prop <- cbind.data.frame(prop, Gamma);

setwd(dir1);
write.table(prop, file="./data/BMR_2spAC_RS-IF.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
############# AC; Apparent competition ====
AC <- read.table("./data/AC_Tot_data.txt", h=T,  sep = "\t");
Alpha <- c(0.1, 0.2, 0.5, 0.75, 1, 1.5, 2, 5, 10); # BMR Ri:Rr
Beta  <- c(1, 2, 5, 10);  # BMR Cr:Ri

#### % Biodiversity change ====
Bio <- as.numeric(names(summary(as.factor(AC$BioD)) ));
prop <- as.data.frame(matrix(0, nrow= length(Alpha), ncol=length(Bio)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Alpha)){
  sub <- AC[which(AC$Alpha == Alpha[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Beta)) ))
  mean_vec <- c();
  for(j in 1:length(Bio)){
    bio_propvec <- c();
    for(k in 1:length(bmr_vec)){
      bio_prop <- Proportion( dim(sub[which(sub$Beta == bmr_vec[k] & sub$BioD == Bio[j]),])[1], tot);
      bio_propvec <- append(bio_propvec, bio_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(bio_propvec)); # Average proportion over for each BioD and Alpha
  }
  prop[i,] <- mean_vec;
}
colnames(prop) <- Bio;
prop <- cbind.data.frame(prop, Alpha);

setwd(dir1);
write.table(prop, file="./data/BMR_AC_Bio.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Community state (EQf) ==== 
Compo <- as.character(names(summary(as.factor(AC$Dyn_state_f)) ));
prop <- as.data.frame(matrix(0, nrow= length(Alpha), ncol=length(Compo)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Alpha)){
  sub <- AC[which(AC$Alpha == Alpha[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Beta)) ))
  mean_vec <- c();
  for(j in 1:length(Compo)){
    compo_propvec <- c();
    for(k in 1:length(bmr_vec)){
      compo_prop <- Proportion( dim(sub[which(sub$Beta == bmr_vec[k] & sub$Dyn_state_f == Compo[j]),])[1], tot);
      compo_propvec <- append(compo_propvec, compo_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(compo_propvec)); # Average proportion over for each BioD and Alpha
  }
  prop[i,] <- mean_vec;
}
colnames(prop) <- Compo;
prop <- cbind.data.frame(prop, Alpha);

write.table(prop, file="./data/BMR_AC_EQf.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
#### % Invasion states ====
Mec <- as.character(names(summary(as.factor(AC$State)) ));
prop <- as.data.frame(matrix(0, nrow= length(Alpha), ncol=length(Mec)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Alpha)){
  sub <- AC[which(AC$Alpha == Alpha[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Beta)) ))
  mean_vec <- c();
  for(j in 1:length(Mec)){
    mec_propvec <- c();
    for(k in 1:length(bmr_vec)){
      mec_prop <- Proportion( dim(sub[which(sub$Beta == bmr_vec[k] & sub$State == Mec[j]),])[1], tot);
      mec_propvec <- append(mec_propvec, mec_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(mec_propvec)); # Average proportion over for each BioD and Alpha
  }
  prop[i,] <- mean_vec;
}
colnames(prop) <- Mec;

Bio_gain <- apply(cbind.data.frame(prop$Integration, prop$Occupancy, prop$Rescue), 1, sum);
Bio_neutral <- apply(cbind.data.frame(prop$Resistance, prop$Substitution), 1, sum);
Bio_loss <- prop$Vulnerability

prop <- cbind.data.frame(prop, Bio_gain, Bio_neutral, Bio_loss, Alpha);

write.table(prop, file="./data/BMR_AC_IS.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
#### % Stability (initial vs after 5000 years) ====
EQIF <- as.character(names(summary(as.factor(AC$EQ_transit_InitFin)) ));
prop <- as.data.frame(matrix(0, nrow= length(Alpha), ncol=length(EQIF)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Alpha)){
  sub <- AC[which(AC$Alpha == Alpha[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Beta)) ))
  mean_vec <- c();
  for(j in 1:length(EQIF)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Beta == bmr_vec[k] & sub$EQ_transit_InitFin == EQIF[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each BioD and Alpha
  }
  prop[i,] <- mean_vec;
}
colnames(prop) <- EQIF;
prop <- cbind.data.frame(prop, Alpha);

write.table(prop, file="./data/BMR_AC_RS-IF.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
#### % Stability (resident vs invaded community after 5000 years) ====
EQRI <- as.character(names(summary(as.factor(AC$EQ_transit_ResInv)) ));
prop <- as.data.frame(matrix(0, nrow= length(Alpha), ncol=length(EQRI)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Alpha)){
  sub <- AC[which(AC$Alpha == Alpha[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Beta)) ))
  mean_vec <- c();
  for(j in 1:length(EQRI)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Beta == bmr_vec[k] & sub$EQ_transit_ResInv == EQRI[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each BioD and Alpha
  }
  prop[i,] <- mean_vec;
}
colnames(prop) <- EQRI;

DeltaS_gain <- apply(cbind.data.frame(prop$Cycles.Eq, prop$Null.Cycles), 1, sum);
DeltaS_neutral <- apply(cbind.data.frame(prop$Null.Null, prop$Cycles.Cycles, prop$Eq.Eq), 1, sum);
DeltaS_loss <- prop$Cycles.Null;

prop <- cbind.data.frame(prop, DeltaS_gain, DeltaS_neutral, DeltaS_loss, Alpha);

write.table(prop, file="./data/BMR_AC_RS-ResInv.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);
############# EC; Ecploitative competition ====
EC <- read.table("./data/EC_Tot_data.txt", h=T,  sep = "\t");
Alpha <- c(1, 2, 5, 10); # BMR Ci:Rr
Beta <- c(0.1, 0.2, 0.5, 0.75, 1, 1.5, 2, 5, 10); # BMR Cr:Ci

#### % Biodiversity change ====
Bio <- as.numeric(names(summary(as.factor(EC$BioD)) ));
prop <- as.data.frame(matrix(0, nrow= length(Alpha), ncol=length(Bio)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Beta)){
  sub <- EC[which(EC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Bio)){
    bio_propvec <- c();
    for(k in 1:length(bmr_vec)){
      bio_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$BioD == Bio[j]),])[1], tot);
      bio_propvec <- append(bio_propvec, bio_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(bio_propvec)); # Average proportion over for each BioD and Alpha
  }
  prop[i,] <- mean_vec;
}
colnames(prop) <- Bio;
prop <- cbind.data.frame(prop, Beta);

write.table(prop, file="./data/BMR_EC_Bio.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Community change (EQf) ====
Compo <- as.character(names(summary(as.factor(EC$Dyn_state_f)) ));
prop <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(Compo)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Beta)){
  sub <- EC[which(EC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Compo)){
    compo_propvec <- c();
    for(k in 1:length(bmr_vec)){
      compo_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$Dyn_state_f == Compo[j]),])[1], tot);
      compo_propvec <- append(compo_propvec, compo_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(compo_propvec)); # Average proportion over for each BioD and Alpha
  }
  prop[i,] <- mean_vec;
}
colnames(prop) <- Compo;
prop <- cbind.data.frame(prop, Beta);

write.table(prop, file="./data/BMR_EC_EQf.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Invasion states ====
Mec <- as.character(names(summary(as.factor(EC$State)) ));
prop <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(Mec)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Beta)){
  sub <- EC[which(EC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Mec)){
    mec_propvec <- c();
    for(k in 1:length(bmr_vec)){
      mec_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$State == Mec[j]),])[1], tot);
      mec_propvec <- append(mec_propvec, mec_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(mec_propvec)); # Average proportion over for each BioD and Alpha
  }
  prop[i,] <- mean_vec;
}
colnames(prop) <- Mec;

Bio_gain <- apply(cbind.data.frame(prop$Integration, prop$Occupancy, prop$Rescue), 1, sum);
Bio_neutral <- apply(cbind.data.frame(prop$Resistance, prop$Substitution), 1, sum);
Bio_loss <- prop$Vulnerability
prop <- cbind.data.frame(prop, Bio_gain, Bio_neutral, Bio_loss, Beta);

write.table(prop, file="./data/BMR_EC_IS.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Stability (initial vs after 5000 years) ====
EQIF <- as.character(names(summary(as.factor(EC$EQ_transit_InitFin)) ));
prop <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(EQIF)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Beta)){
  sub <- EC[which(EC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(EQIF)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$EQ_transit_InitFin == EQIF[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each BioD and Alpha
  }
  prop[i,] <- mean_vec;
}
colnames(prop) <- EQIF;
prop <- cbind.data.frame(prop, Beta);

write.table(prop, file="./data/BMR_EC_RS-IF.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Stability (resident vs invaded community after 5000 years) ====
EQRI <- as.character(names(summary(as.factor(EC$EQ_transit_ResInv)) ));
prop <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(EQRI)));
tot <- 80200; # Total of data point (T*Car) for a given bmr

for(i in 1:length(Beta)){
  sub <- EC[which(EC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(EQRI)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$EQ_transit_ResInv == EQRI[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each BioD and Alpha
  }
  prop[i,] <- mean_vec;
}
colnames(prop) <- EQRI;

DeltaS_gain <- apply(cbind.data.frame(prop$Cycles.Eq, prop$Null.Cycles, prop$Null.Eq), 1, sum);
DeltaS_neutral <- apply(cbind.data.frame(prop$Null.Null, prop$Cycles.Cycles, prop$Eq.Eq), 1, sum);
DeltaS_loss <- apply(cbind.data.frame(prop$Cycles.Null, prop$Eq.Cycles, prop$Eq.Null), 1, sum);
prop <- cbind.data.frame(prop, DeltaS_gain, DeltaS_neutral, DeltaS_loss, Beta);

write.table(prop, file="./data/BMR_EC_RS-ResInv.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

############# TC; Trophic chain ====
TC <- read.table("./data/TC_Tot_data.txt", h=T,  sep = "\t");
TC <- TC[which(TC$BioD != "NA"),] #removing bug lines
Delta <- rep(0, length=dim(TC)[1] )
TC <- cbind.data.frame(TC, Delta)
TC$Delta <- (TC$Beta/TC$Alpha);

Alpha <- c(1, 2, 5, 10); # BMR Cr:Rr
Beta <- c(1, 2, 5, 10);  # BMR Pi:Cr
Gamma <- c(1, 2, 4, 5, 10, 20, 25, 50, 100)
Ddelta <- c(0.1, 0.2, 0.4, 0.5, 1, 2, 2.5, 5, 10);
tot <- 80200; # Total of data point (T*Car) for a given bmr
#### % Biodiversity change ====
Bio <- as.numeric(names(summary(as.factor(TC$BioD)) ));

# Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(Bio)));
for(i in 1:length(Gamma)){
  sub <- TC[which(TC$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(Bio)){
    bio_propvec <- c();
    for(k in 1:length(bmr_vec)){
      bio_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$BioD == Bio[j]),])[1], tot);
      bio_propvec <- append(bio_propvec, bio_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(bio_propvec)); # Average proportion over for each BioD and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- Bio;
prop1 <- cbind.data.frame(prop1, Gamma);

write.table(prop1, file="./data/BMR_TC_Bio_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(Bio)));
for(i in 1:length(Ddelta)){
  sub <- TC[which(TC$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(Bio)){
    bio_propvec <- c();
    for(k in 1:length(bmr_vec)){
      bio_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$BioD == Bio[j]),])[1], tot);
      bio_propvec <- append(bio_propvec, bio_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(bio_propvec)); # Average proportion over for each BioD and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- Bio;
prop2 <- cbind.data.frame(prop2, Ddelta);

setwd(dir1);
write.table(prop2, file="./data/BMR_TC_Bio_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(Bio)));

for(i in 1:length(Beta)){
  sub <- TC[which(TC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Bio)){
    bio_propvec <- c();
    for(k in 1:length(bmr_vec)){
      bio_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$BioD == Bio[j]),])[1], tot);
      bio_propvec <- append(bio_propvec, bio_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(bio_propvec)); # Average proportion over for each BioD and Alpha
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- Bio;
prop3 <- cbind.data.frame(prop3, Beta);

setwd(dir1);
write.table(prop3, file="./data/BMR_TC_Bio_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Community change (EQf) ====
Compo <- as.character(names(summary(as.factor(TC$Dyn_state_f)) ));

# Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(Compo)));
for(i in 1:length(Gamma)){
  sub <- TC[which(TC$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(Compo)){
    comp_propvec <- c();
    for(k in 1:length(bmr_vec)){
      comp_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$Dyn_state_f == Compo[j]),])[1], tot);
      comp_propvec <- append(comp_propvec, comp_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(comp_propvec)); # Average proportion over for each Dyn_state_f and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- Compo;
prop1 <- cbind.data.frame(prop1, Gamma);

write.table(prop1, file="./data/BMR_TC_EQf_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(Compo)));
for(i in 1:length(Ddelta)){
  sub <- TC[which(TC$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(Compo)){
    comp_propvec <- c();
    for(k in 1:length(bmr_vec)){
      comp_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$Dyn_state_f == Compo[j]),])[1], tot);
      comp_propvec <- append(comp_propvec, comp_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(comp_propvec)); # Average proportion over for each Dyn_state_f and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- Compo;
prop2 <- cbind.data.frame(prop2, Ddelta);

write.table(prop2, file="./data/BMR_TC_EQf_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(Compo)));
for(i in 1:length(Beta)){
  sub <- TC[which(TC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Compo)){
    comp_propvec <- c();
    for(k in 1:length(bmr_vec)){
      comp_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$Dyn_state_f == Compo[j]),])[1], tot);
      comp_propvec <- append(comp_propvec, comp_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(comp_propvec)); # Average proportion over for each Dyn_state_f and Alphha
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- Compo;
prop3 <- cbind.data.frame(prop3, Beta);

write.table(prop3, file="./data/BMR_TC_EQf_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Invasion states ====
Mec <- as.character(names(summary(as.factor(TC$State)) ));

# Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(Mec)));
for(i in 1:length(Gamma)){
  sub <- TC[which(TC$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(Mec)){
    mec_propvec <- c();
    for(k in 1:length(bmr_vec)){
      mec_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$State == Mec[j]),])[1], tot);
      mec_propvec <- append(mec_propvec, mec_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(mec_propvec)); # Average proportion over for each Mec and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- Mec;

Bio_gain <- apply(cbind.data.frame(prop1$Integration, prop1$Rescue), 1, sum);
Bio_neutral <- prop1$Resistance;
Bio_loss <- prop1$Vulnerability
prop1 <- cbind.data.frame(prop1, Bio_gain, Bio_neutral, Bio_loss, Gamma);

write.table(prop1, file="./data/BMR_TC_IS_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(Mec)));
for(i in 1:length(Ddelta)){
  sub <- TC[which(TC$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(Mec)){
    mec_propvec <- c();
    for(k in 1:length(bmr_vec)){
      mec_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$State == Mec[j]),])[1], tot);
      mec_propvec <- append(mec_propvec, mec_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(mec_propvec)); # Average proportion over for each Mec and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- Mec;

Bio_gain <- apply(cbind.data.frame(prop2$Integration, prop2$Rescue), 1, sum);
Bio_neutral <- prop2$Resistance;
Bio_loss <- prop2$Vulnerability

prop2 <- cbind.data.frame(prop2, Bio_gain, Bio_neutral, Bio_loss, Ddelta);

write.table(prop2, file="./data/BMR_TC_IS_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(Mec)));
for(i in 1:length(Beta)){
  sub <- TC[which(TC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Mec)){
    mec_propvec <- c();
    for(k in 1:length(bmr_vec)){
      mec_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$State == Mec[j]),])[1], tot);
      mec_propvec <- append(mec_propvec, mec_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(mec_propvec)); # Average proportion over for each Mec and Alpha
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- Mec;

Bio_gain <- apply(cbind.data.frame(prop3$Integration, prop3$Rescue), 1, sum);
Bio_neutral <- prop3$Resistance;
Bio_loss <- prop3$Vulnerability

prop3 <- cbind.data.frame(prop3, Bio_gain, Bio_neutral, Bio_loss, Beta);

write.table(prop3, file="./data/BMR_TC_IS_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Stability (initial vs after 5000 years) ====
EQIF <- as.character(names(summary(as.factor(TC$EQ_transit_InitFin)) ));
# over Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(EQIF)));
for(i in 1:length(Gamma)){
  sub <- TC[which(TC$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(EQIF)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$EQ_transit_InitFin == EQIF[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQIF and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- EQIF;
prop1 <- cbind.data.frame(prop1, Gamma);

write.table(prop1, file="./data/BMR_TC_RS-IF_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# over Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(EQIF)));

for(i in 1:length(Ddelta)){
  sub <- TC[which(TC$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(EQIF)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$EQ_transit_InitFin == EQIF[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQIF and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- EQIF;
prop2 <- cbind.data.frame(prop2, Ddelta);

write.table(prop2, file="./data/BMR_TC_RS-IF_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# over Gamma gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(EQIF)));

for(i in 1:length(Beta)){
  sub <- TC[which(TC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(EQIF)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$EQ_transit_InitFin == EQIF[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQIF and Alpha
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- EQIF;
prop3 <- cbind.data.frame(prop3, Beta);

write.table(prop3, file="./data/BMR_TC_RS-IF_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Stability (resident vs invaded community after 5000 years) ====
EQRI <- as.character(names(summary(as.factor(TC$EQ_transit_ResInv)) ));
# over Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(EQRI)));
for(i in 1:length(Gamma)){
  sub <- TC[which(TC$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(EQRI)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$EQ_transit_ResInv == EQRI[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQRI and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- EQRI;

DeltaS_gain <- apply(cbind.data.frame(prop1$Cycles.Eq, prop1$Null.Cycles, prop1$Null.Eq), 1, sum);
DeltaS_neutral <- apply(cbind.data.frame(prop1$Null.Null, prop1$Cycles.Cycles, prop1$Eq.Eq), 1, sum);
DeltaS_loss <- apply(cbind.data.frame(prop1$Cycles.Null, prop1$Eq.Cycles, prop1$Eq.Null), 1, sum);
prop1 <- cbind.data.frame(prop1, DeltaS_gain, DeltaS_neutral, DeltaS_loss, Gamma);

write.table(prop1, file="./data/BMR_TC_RS-ResInv_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# over Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(EQRI)));
for(i in 1:length(Ddelta)){
  sub <- TC[which(TC$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(EQRI)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$EQ_transit_ResInv == EQRI[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQRI and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- EQRI;
DeltaS_gain <- apply(cbind.data.frame(prop2$Cycles.Eq, prop2$Null.Cycles, prop2$Null.Eq), 1, sum);
DeltaS_neutral <- apply(cbind.data.frame(prop2$Null.Null, prop2$Cycles.Cycles, prop2$Eq.Eq), 1, sum);
DeltaS_loss <- apply(cbind.data.frame(prop2$Cycles.Null, prop2$Eq.Cycles, prop2$Eq.Null), 1, sum);
prop2 <- cbind.data.frame(prop2, DeltaS_gain, DeltaS_neutral, DeltaS_loss, Ddelta);

write.table(prop2, file="./data/BMR_TC_RS-ResInv_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# over Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(EQRI)));
for(i in 1:length(Beta)){
  sub <- TC[which(TC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(EQRI)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$EQ_transit_ResInv == EQRI[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQRI and Alpha
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- EQRI;
DeltaS_gain <- apply(cbind.data.frame(prop3$Cycles.Eq, prop3$Null.Cycles, prop3$Null.Eq), 1, sum);
DeltaS_neutral <- apply(cbind.data.frame(prop3$Null.Null, prop3$Cycles.Cycles, prop3$Eq.Eq), 1, sum);
DeltaS_loss <- apply(cbind.data.frame(prop3$Cycles.Null, prop3$Eq.Cycles, prop3$Eq.Null), 1, sum);
prop3 <- cbind.data.frame(prop3, DeltaS_gain, DeltaS_neutral, DeltaS_loss, Beta);

write.table(prop3, file="./data/BMR_TC_RS-ResInv_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

############# IGP; Intraguild-Predation ====
IGPB <- read.table("./data/IGPB_Tot_data.txt", h=T,  sep = "\t");
Delta <- rep(0, length=dim(IGPB)[1] )
IGPB <- cbind.data.frame(IGPB, Delta)
IGPB$Beta <- (1/IGPB$Beta);
IGPB$Delta <- (IGPB$Beta/IGPB$Alpha);

IGPC <- read.table("./data/IGPC_Tot_data.txt", h=T,  sep = "\t");
Delta <- rep(0, length=dim(IGPC)[1] )
IGPC <- cbind.data.frame(IGPC, Delta)
IGPC$Delta <- (IGPC$Beta/IGPC$Alpha);

Beta <- c(0.1,0.2,0.5, 1, 2, 5, 10);
tot <- 80200; # Total of data point (T*Car) for a given bmr

# Rename columns with predator => consumer to create IGP
names(IGPB)
colnames(IGPB)[c(11,14,17,20,24,28)] <- c("M_Cr", "Cr_i", "Cr_f", "Cr_fin", "Cr_extime", "Cr_eq");
colnames(IGPC)[c(11,14,17,20,24,28)] <- c("M_Ci", "Ci_i", "Ci_f", "Ci_fin", "Ci_extime", "Ci_eq");
#igpc.names <- names(IGPC)[c(1:9,11,10,12,14,13,15,17,16,18,20,19,21:22,24,23,25:26,28,27,29:42)];
IGPC <- cbind.data.frame(IGPC[,c(1:9,11,10,12,14,13,15,17,16,18,20,19,21:22,24,23,25:26,28,27,29:42)])

IGP <- rbind.data.frame(IGPB, IGPC)

#### % Invasion states ====
Mec <- as.character(names(summary(as.factor(IGP$State)) ));

prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(Mec)));
for(i in 1:length(Beta)){
  sub <- IGP[which(IGP$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Mec)){
    mec_propvec <- c();
    tot1 <- ifelse(Beta[i] == 1, tot*2, tot);
    for(k in 1:length(bmr_vec)){
      mec_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$State == Mec[j]),])[1], tot1);
      mec_propvec <- append(mec_propvec, mec_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(mec_propvec)); # Average proportion over for each Mec and Beta
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- Mec;

Bio_gain <- apply(cbind.data.frame(prop3$Integration, prop3$Occupancy, prop3$Rescue), 1, sum);
Bio_neutral <- apply(cbind.data.frame(prop3$Resistance, prop3$Substitution), 1, sum);
Bio_loss <- prop3$Vulnerability
prop3 <- cbind.data.frame(prop3, Bio_gain, Bio_neutral, Bio_loss, Beta);

write.table(prop3, file="./data/BMR_IGP_IS_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Stability (resident vs invaded community after 5000 years) ====
EQRI <- as.character(names(summary(as.factor(IGP$EQ_transit_ResInv)) ));

# Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(EQRI)));
for(i in 1:length(Beta)){
  sub <- IGP[which(IGP$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  tot1 <- ifelse(Beta[i] == 1, tot*2, tot);
  for(j in 1:length(EQRI)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$EQ_transit_ResInv == EQRI[j]),])[1], tot1);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQRI and Beta
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- EQRI;

DeltaS_gain <- apply(cbind.data.frame(prop3$Cycles.Eq, prop3$Null.Cycles, prop3$Cycles.Eq), 1, sum);
DeltaS_neutral <- apply(cbind.data.frame(prop3$Null.Null, prop3$Cycles.Cycles, prop3$Eq.Eq), 1, sum);
DeltaS_loss <- apply(cbind.data.frame(prop3$Cycles.Null), 1, sum);
prop3 <- cbind.data.frame(prop3, DeltaS_gain, DeltaS_neutral, DeltaS_loss, Beta);

write.table(prop3, file="./data/BMR_IGP_RS-ResInv_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

############# IGPB; Intraguild-prey ====
IGPB <- read.table("./data/IGPB_Tot_data.txt", h=T,  sep = "\t");
Delta <- rep(0, length=dim(IGPB)[1] )
IGPB <- cbind.data.frame(IGPB, Delta)
IGPB$Delta <- (IGPB$Beta/IGPB$Alpha);

Alpha <- c(1, 2, 5, 10); # BMR Cr:Rr
Beta <- c(1, 2, 5, 10);  # BMR Pi:Cr
Gamma <- c(1, 2, 4, 5, 10, 20, 25, 50, 100)
Ddelta <- c(0.1, 0.2, 0.4, 0.5, 1, 2, 2.5, 5, 10);
tot <- 80200; # Total of data point (T*Car) for a given bmr

#### % Biodiversity change ====
Bio <- as.numeric(names(summary(as.factor(IGPB$BioD)) ));
# Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(Bio)));
for(i in 1:length(Gamma)){
  sub <- IGPB[which(IGPB$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(Bio)){
    bio_propvec <- c();
    for(k in 1:length(bmr_vec)){
      bio_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$BioD == Bio[j]),])[1], tot);
      bio_propvec <- append(bio_propvec, bio_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(bio_propvec)); # Average proportion over for each BioD and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- Bio;
prop1 <- cbind.data.frame(prop1, Gamma);

write.table(prop1, file="./data/BMR_IGPB_Bio_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(Bio)));
for(i in 1:length(Ddelta)){
  sub <- IGPB[which(IGPB$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(Bio)){
    bio_propvec <- c();
    for(k in 1:length(bmr_vec)){
      bio_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$BioD == Bio[j]),])[1], tot);
      bio_propvec <- append(bio_propvec, bio_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(bio_propvec)); # Average proportion over for each BioD and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- Bio;
prop2 <- cbind.data.frame(prop2, Ddelta);

write.table(prop2, file="./data/BMR_IGPB_Bio_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# 1/Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(Bio)));
for(i in 1:length(Beta)){
  sub <- IGPB[which(IGPB$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Bio)){
    bio_propvec <- c();
    for(k in 1:length(bmr_vec)){
      bio_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$BioD == Bio[j]),])[1], tot);
      bio_propvec <- append(bio_propvec, bio_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(bio_propvec)); # Average proportion over for each BioD and Gamma
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- Bio;
prop3 <- cbind.data.frame(prop3, (1/Beta));
colnames(prop3) <- c(names(prop3)[1:4], "Beta")

write.table(prop3, file="./data/BMR_IGPB_Bio_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Community change (EQf) ====
Compo <- as.character(names(summary(as.factor(IGPB$Dyn_state_f)) ));
# Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(Compo)));
for(i in 1:length(Gamma)){
  sub <- IGPB[which(IGPB$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(Compo)){
    comp_propvec <- c();
    for(k in 1:length(bmr_vec)){
      comp_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$Dyn_state_f == Compo[j]),])[1], tot);
      comp_propvec <- append(comp_propvec, comp_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(comp_propvec)); # Average proportion over for each Dyn_state_f and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- Compo;
prop1 <- cbind.data.frame(prop1, Gamma);

write.table(prop1, file="./data/BMR_IGPB_EQf_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(Compo)));
for(i in 1:length(Ddelta)){
  sub <- IGPB[which(IGPB$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(Compo)){
    comp_propvec <- c();
    for(k in 1:length(bmr_vec)){
      comp_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$Dyn_state_f == Compo[j]),])[1], tot);
      comp_propvec <- append(comp_propvec, comp_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(comp_propvec)); # Average proportion over for each Dyn_state_f and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- Compo;
prop2 <- cbind.data.frame(prop2, Ddelta);

write.table(prop2, file="./data/BMR_IGPB_EQf_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# 1/Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(Compo)));
for(i in 1:length(Beta)){
  sub <- IGPB[which(IGPB$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Compo)){
    comp_propvec <- c();
    for(k in 1:length(bmr_vec)){
      comp_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$Dyn_state_f == Compo[j]),])[1], tot);
      comp_propvec <- append(comp_propvec, comp_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(comp_propvec)); # Average proportion over for each Dyn_state_f and Alpha
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- Compo;
prop3 <- cbind.data.frame(prop3, (1/Beta));
colnames(prop3) <- c(names(prop3)[1:6], "Beta")

write.table(prop3, file="./data/BMR_IGPB_EQf_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Invasion states ====
Mec <- as.character(names(summary(as.factor(IGPB$State)) ));

# Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(Mec)));
for(i in 1:length(Gamma)){
  sub <- IGPB[which(IGPB$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(Mec)){
    mec_propvec <- c();
    for(k in 1:length(bmr_vec)){
      mec_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$State == Mec[j]),])[1], tot);
      mec_propvec <- append(mec_propvec, mec_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(mec_propvec)); # Average proportion over for each Mec and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- Mec;

Bio_gain <- apply(cbind.data.frame(prop1$Integration, prop1$Occupancy, prop1$Rescue), 1, sum);
Bio_neutral <- apply(cbind.data.frame(prop1$Resistance, prop1$Substitution), 1, sum);
Bio_loss <- prop1$Vulnerability;
prop1 <- cbind.data.frame(prop1, Bio_gain, Bio_neutral, Bio_loss, Gamma);

write.table(prop1, file="./data/BMR_IGPB_IS_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(Mec)));
for(i in 1:length(Ddelta)){
  sub <- IGPB[which(IGPB$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(Mec)){
    mec_propvec <- c();
    for(k in 1:length(bmr_vec)){
      mec_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$State == Mec[j]),])[1], tot);
      mec_propvec <- append(mec_propvec, mec_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(mec_propvec)); # Average proportion over for each Mec and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- Mec;

Bio_gain <- apply(cbind.data.frame(prop2$Integration, prop2$Occupancy, prop2$Rescue), 1, sum);
Bio_neutral <- apply(cbind.data.frame(prop2$Resistance, prop2$Substitution), 1, sum);
Bio_loss <- prop2$Vulnerability;
prop2 <- cbind.data.frame(prop2, Bio_gain, Bio_neutral, Bio_loss, Ddelta);

write.table(prop2, file="./data/BMR_IGPB_IS_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# 1/Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(Mec)));
for(i in 1:length(Beta)){
  sub <- IGPB[which(IGPB$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Mec)){
    mec_propvec <- c();
    for(k in 1:length(bmr_vec)){
      mec_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$State == Mec[j]),])[1], tot);
      mec_propvec <- append(mec_propvec, mec_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(mec_propvec)); # Average proportion over for each Mec and Beta
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- Mec;

Bio_gain <- apply(cbind.data.frame(prop3$Integration, prop3$Occupancy, prop3$Rescue), 1, sum);
Bio_neutral <- apply(cbind.data.frame(prop3$Resistance, prop3$Substitution), 1, sum);
Bio_loss <- prop3$Vulnerability;
prop3 <- cbind.data.frame(prop3, Bio_gain, Bio_neutral, Bio_loss, (1/Beta));
colnames(prop3) <- c(names(prop3)[1:9], "Beta")

write.table(prop3, file="./data/BMR_IGPB_IS_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Stability (initial vs after 5000 years) ====
EQIF <- as.character(names(summary(as.factor(IGPB$EQ_transit_InitFin)) ));
# Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(EQIF)));

for(i in 1:length(Gamma)){
  sub <- IGPB[which(IGPB$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(EQIF)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$EQ_transit_InitFin == EQIF[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQIF and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- EQIF;
prop1 <- cbind.data.frame(prop1, Gamma);

write.table(prop1, file="./data/BMR_IGPB_RS-IF_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(EQIF)));
for(i in 1:length(Ddelta)){
  sub <- IGPB[which(IGPB$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(EQIF)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$EQ_transit_InitFin == EQIF[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQIF and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- EQIF;
prop2 <- cbind.data.frame(prop2, Ddelta);

write.table(prop2, file="./data/BMR_IGPB_RS-IF_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# 1/Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(EQIF)));
for(i in 1:length(Beta)){
  sub <- IGPB[which(IGPB$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(EQIF)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$EQ_transit_InitFin == EQIF[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQIF and Beta
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- EQIF;
prop3 <- cbind.data.frame(prop3, (1/Beta));
colnames(prop3) <- c(names(prop3)[1:length(EQIF)], "Beta")

write.table(prop3, file="./data/BMR_IGPB_RS-IF_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Stability (resident vs invaded community after 5000 years) ====
EQRI <- as.character(names(summary(as.factor(IGPB$EQ_transit_ResInv)) ));
# Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(EQRI)));
for(i in 1:length(Gamma)){
  sub <- IGPB[which(IGPB$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(EQRI)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$EQ_transit_ResInv == EQRI[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQRI and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- EQRI;

DeltaS_gain <- prop1$Null.Cycles
DeltaS_neutral <- apply(cbind.data.frame(prop1$Null.Null, prop1$Cycles.Cycles, prop1$Eq.Eq), 1, sum);
DeltaS_loss <- prop1$Cycles.Null
prop1 <- cbind.data.frame(prop1, DeltaS_gain, DeltaS_neutral, DeltaS_loss, Gamma);

write.table(prop1, file="./data/BMR_IGPB_RS-ResInv_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(EQRI)));
for(i in 1:length(Ddelta)){
  sub <- IGPB[which(IGPB$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(EQRI)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$EQ_transit_ResInv == EQRI[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQRI and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- EQRI;

DeltaS_gain <- prop2$Null.Cycles
DeltaS_neutral <- apply(cbind.data.frame(prop2$Null.Null, prop2$Cycles.Cycles, prop2$Eq.Eq), 1, sum);
DeltaS_loss <- prop2$Cycles.Null
prop2 <- cbind.data.frame(prop2, DeltaS_gain, DeltaS_neutral, DeltaS_loss, Ddelta);

write.table(prop2, file="./data/BMR_IGPB_RS-ResInv_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# 1/Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(EQRI)));
for(i in 1:length(Beta)){
  sub <- IGPB[which(IGPB$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(EQRI)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$EQ_transit_ResInv == EQRI[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQRI and Beta
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- EQRI;

DeltaS_gain <- prop3$Null.Cycles;
DeltaS_neutral <- apply(cbind.data.frame(prop3$Null.Null, prop3$Cycles.Cycles, prop3$Eq.Eq), 1, sum);
DeltaS_loss <- prop3$Cycles.Null;
prop3 <- cbind.data.frame(prop3, DeltaS_gain, DeltaS_neutral, DeltaS_loss, (1/Beta));
colnames(prop3) <- c(names(prop3)[1:9], "Beta")

write.table(prop3, file="./data/BMR_IGPB_RS-ResInv_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

############# IGPC; Intraguild-predator ====
IGPC <- read.table("./data/IGPC_Tot_data.txt", h=T,  sep = "\t");
Delta <- rep(0, length=dim(IGPC)[1] )
IGPC <- cbind.data.frame(IGPC, Delta)
IGPC$Delta <- (IGPC$Beta/IGPC$Alpha);

Alpha <- c(1, 2, 5, 10); # BMR Cr:Rr
Beta <- c(1, 2, 5, 10);  # BMR Pi:Cr
Gamma <- c(1, 2, 4, 5, 10, 20, 25, 50, 100)
Ddelta <- c(0.1, 0.2, 0.4, 0.5, 1, 2, 2.5, 5, 10);
tot <- 80200; # Total of data point (T*Car) for a given bmr

#### % Biodiversity change ====
Bio <- as.numeric(names(summary(as.factor(IGPC$BioD)) ));
# Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(Bio)));

for(i in 1:length(Gamma)){
  sub <- IGPC[which(IGPC$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(Bio)){
    bio_propvec <- c();
    for(k in 1:length(bmr_vec)){
      bio_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$BioD == Bio[j]),])[1], tot);
      bio_propvec <- append(bio_propvec, bio_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(bio_propvec)); # Average proportion over for each BioD and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- Bio;
prop1 <- cbind.data.frame(prop1, Gamma);

write.table(prop1, file="./data/BMR_IGPC_Bio_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(Bio)));
for(i in 1:length(Ddelta)){
  sub <- IGPC[which(IGPC$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(Bio)){
    bio_propvec <- c();
    for(k in 1:length(bmr_vec)){
      bio_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$BioD == Bio[j]),])[1], tot);
      bio_propvec <- append(bio_propvec, bio_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(bio_propvec)); # Average proportion over for each BioD and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- Bio;
prop2 <- cbind.data.frame(prop2, Ddelta);

write.table(prop2, file="./data/BMR_IGPC_Bio_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(Bio)));
for(i in 1:length(Beta)){
  sub <- IGPC[which(IGPC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Bio)){
    bio_propvec <- c();
    for(k in 1:length(bmr_vec)){
      bio_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$BioD == Bio[j]),])[1], tot);
      bio_propvec <- append(bio_propvec, bio_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(bio_propvec)); # Average proportion over for each BioD and Gamma
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- Bio;
prop3 <- cbind.data.frame(prop3, Beta);

write.table(prop3, file="./data/BMR_IGPC_Bio_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Community change (EQf) ====
Compo <- as.character(names(summary(as.factor(IGPC$Dyn_state_f)) ));
# Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(Compo)));
for(i in 1:length(Gamma)){
  sub <- IGPC[which(IGPC$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(Compo)){
    comp_propvec <- c();
    for(k in 1:length(bmr_vec)){
      comp_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$Dyn_state_f == Compo[j]),])[1], tot);
      comp_propvec <- append(comp_propvec, comp_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(comp_propvec)); # Average proportion over for each Compo and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- Compo;
prop1 <- cbind.data.frame(prop1, Gamma);

write.table(prop1, file="./data/BMR_IGPC_EQf_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(Compo)));
for(i in 1:length(Ddelta)){
  sub <- IGPC[which(IGPC$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(Compo)){
    comp_propvec <- c();
    for(k in 1:length(bmr_vec)){
      comp_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$Dyn_state_f == Compo[j]),])[1], tot);
      comp_propvec <- append(comp_propvec, comp_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(comp_propvec)); # Average proportion over for each Compo and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- Compo;
prop2 <- cbind.data.frame(prop2, Ddelta);

write.table(prop2, file="./data/BMR_IGPC_EQf_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(Compo)));
for(i in 1:length(Beta)){
  sub <- IGPC[which(IGPC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Compo)){
    comp_propvec <- c();
    for(k in 1:length(bmr_vec)){
      comp_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$Dyn_state_f == Compo[j]),])[1], tot);
      comp_propvec <- append(comp_propvec, comp_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(comp_propvec)); # Average proportion over for each Compo and Gamma
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- Compo;
prop3 <- cbind.data.frame(prop3, Beta);

write.table(prop3, file="./data/BMR_IGPC_EQf_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Invasion states ====
Mec <- as.character(names(summary(as.factor(IGPC$State)) ));
# Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(Mec)));
for(i in 1:length(Gamma)){
  sub <- IGPC[which(IGPC$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(Mec)){
    mec_propvec <- c();
    for(k in 1:length(bmr_vec)){
      mec_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$State == Mec[j]),])[1], tot);
      mec_propvec <- append(mec_propvec, mec_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(mec_propvec)); # Average proportion over for each Mec and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- Mec;

Bio_gain <- apply(cbind.data.frame(prop1$Integration, prop1$Occupancy), 1, sum);
Bio_neutral <- apply(cbind.data.frame(prop1$Resistance, prop1$Substitution), 1, sum);
Bio_loss <- prop1$Vulnerability
prop1 <- cbind.data.frame(prop1, Bio_gain, Bio_neutral, Bio_loss, Gamma);

write.table(prop1, file="./data/BMR_IGPC_IS_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(Mec)));
for(i in 1:length(Ddelta)){
  sub <- IGPC[which(IGPC$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(Mec)){
    mec_propvec <- c();
    for(k in 1:length(bmr_vec)){
      mec_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$State == Mec[j]),])[1], tot);
      mec_propvec <- append(mec_propvec, mec_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(mec_propvec)); # Average proportion over for each Mec and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- Mec;

Bio_gain <- apply(cbind.data.frame(prop2$Integration, prop2$Occupancy), 1, sum);
Bio_neutral <- apply(cbind.data.frame(prop2$Resistance, prop2$Substitution), 1, sum);
Bio_loss <- prop2$Vulnerability
prop2 <- cbind.data.frame(prop2, Bio_gain, Bio_neutral, Bio_loss, Ddelta);

write.table(prop2, file="./data/BMR_IGPC_IS_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(Mec)));
for(i in 1:length(Beta)){
  sub <- IGPC[which(IGPC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(Mec)){
    mec_propvec <- c();
    for(k in 1:length(bmr_vec)){
      mec_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$State == Mec[j]),])[1], tot);
      mec_propvec <- append(mec_propvec, mec_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(mec_propvec)); # Average proportion over for each Mec and Beta
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- Mec;

Bio_gain <- apply(cbind.data.frame(prop3$Integration, prop3$Occupancy), 1, sum);
Bio_neutral <- apply(cbind.data.frame(prop3$Resistance, prop3$Substitution), 1, sum);
Bio_loss <- prop3$Vulnerability
prop3 <- cbind.data.frame(prop3, Bio_gain, Bio_neutral, Bio_loss, Beta);

write.table(prop3, file="./data/BMR_IGPC_IS_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Stability (initial vs after 5000 years) ====
EQIF <- as.character(names(summary(as.factor(IGPC$EQ_transit_InitFin)) ));
# over Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(EQIF)));
for(i in 1:length(Gamma)){
  sub <- IGPC[which(IGPC$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(EQIF)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$EQ_transit_InitFin == EQIF[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQIF and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- EQIF;
prop1 <- cbind.data.frame(prop1, Gamma);

write.table(prop1, file="./data/BMR_IGPC_RS-IF_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(EQIF)));
for(i in 1:length(Ddelta)){
  sub <- IGPC[which(IGPC$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(EQIF)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$EQ_transit_InitFin == EQIF[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQIF and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- EQIF;
prop2 <- cbind.data.frame(prop2, Ddelta);

write.table(prop2, file="./data/BMR_IGPC_RS-IF_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(EQIF)));
for(i in 1:length(Beta)){
  sub <- IGPC[which(IGPC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(EQIF)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$EQ_transit_InitFin == EQIF[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQIF and Beta
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- EQIF;
prop3 <- cbind.data.frame(prop3, Beta);

write.table(prop3, file="./data/BMR_IGPC_RS-IF_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

#### % Stability (resident vs invaded community after 5000 years) ====
EQRI <- as.character(names(summary(as.factor(IGPC$EQ_transit_ResInv)) ));
# Gamma gradient
prop1 <- as.data.frame(matrix(0, nrow= length(Gamma), ncol=length(EQRI)));
for(i in 1:length(Gamma)){
  sub <- IGPC[which(IGPC$Gamma == Gamma[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Delta)) ))
  mean_vec <- c();
  for(j in 1:length(EQRI)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Delta == bmr_vec[k] & sub$EQ_transit_ResInv == EQRI[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQRI and Delta
  }
  prop1[i,] <- mean_vec;
}
colnames(prop1) <- EQRI;

DeltaS_gain <- apply(cbind.data.frame(prop1$Cycles.Eq, prop1$Null.Cycles, prop1$Null.Eq), 1, sum);
DeltaS_neutral <- apply(cbind.data.frame(prop1$Null.Null, prop1$Cycles.Cycles, prop1$Eq.Eq), 1, sum);
DeltaS_loss <- prop1$Cycles.Null;
prop1 <- cbind.data.frame(prop1, DeltaS_gain, DeltaS_neutral, DeltaS_loss, Gamma);

write.table(prop1, file="./data/BMR_IGPC_RS-ResInv_Gamma.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Delta gradient
prop2 <- as.data.frame(matrix(0, nrow= length(Ddelta), ncol=length(EQRI)));
for(i in 1:length(Ddelta)){
  sub <- IGPC[which(IGPC$Delta == Ddelta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Gamma)) ))
  mean_vec <- c();
  for(j in 1:length(EQRI)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Gamma == bmr_vec[k] & sub$EQ_transit_ResInv == EQRI[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQRI and Gamma
  }
  prop2[i,] <- mean_vec;
}
colnames(prop2) <- EQRI;

DeltaS_gain <- apply(cbind.data.frame(prop2$Cycles.Eq, prop2$Null.Cycles, prop2$Null.Eq), 1, sum);
DeltaS_neutral <- apply(cbind.data.frame(prop2$Null.Null, prop2$Cycles.Cycles, prop2$Eq.Eq), 1, sum);
DeltaS_loss <- prop2$Cycles.Null;
prop2 <- cbind.data.frame(prop2, DeltaS_gain, DeltaS_neutral, DeltaS_loss, Ddelta);

write.table(prop2, file="./data/BMR_IGPC_RS-ResInv_Delta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);

# Beta gradient
prop3 <- as.data.frame(matrix(0, nrow= length(Beta), ncol=length(EQRI)));
for(i in 1:length(Beta)){
  sub <- IGPC[which(IGPC$Beta == Beta[i]),];
  bmr_vec <- as.numeric(names(summary(as.factor(sub$Alpha)) ))
  mean_vec <- c();
  for(j in 1:length(EQRI)){
    eq_propvec <- c();
    for(k in 1:length(bmr_vec)){
      eq_prop <- Proportion( dim(sub[which(sub$Alpha == bmr_vec[k] & sub$EQ_transit_ResInv == EQRI[j]),])[1], tot);
      eq_propvec <- append(eq_propvec, eq_prop); # Proportions for all varied bmr
    }
    mean_vec <- append(mean_vec, mean(eq_propvec)); # Average proportion over for each EQRI and Beta
  }
  prop3[i,] <- mean_vec;
}
colnames(prop3) <- EQRI;

DeltaS_gain <- apply(cbind.data.frame(prop3$Cycles.Eq, prop3$Null.Cycles, prop3$Null.Eq), 1, sum);
DeltaS_neutral <- apply(cbind.data.frame(prop3$Null.Null, prop3$Cycles.Cycles, prop3$Eq.Eq), 1, sum);
DeltaS_loss <- prop3$Cycles.Null;
prop3 <- cbind.data.frame(prop3, DeltaS_gain, DeltaS_neutral, DeltaS_loss, Beta);

write.table(prop3, file="./data/BMR_IGPC_RS-ResInv_Beta.txt", dec= ".", sep = "\t", row.names=FALSE, col.names=TRUE);