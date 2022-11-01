####################### Functions.R
# script gathering the functions required for the simulations and those used for other analyses.

######################## Parameter settings according to temperature- and mass- dependencies ====
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

######################## Function to calculate Proportions ====
Proportion <- function(xx, yy){
  prop <- (xx/yy)*100;
  return(prop)
}

######################## Functions for Temperature- and mass dependent rates ====
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

########################  Functions to calculate population biomasses at equilibrium ====
#### R* rule function for EC module
R.Rule.Expcomp <- function(Temp, Car, MRr, MCi, MCr){
  #### temperature- and mass-dependent growth rate of basal prey #
  r_r= exp(Ir)*MRr^Sr*exp(Ear*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  #### temperature- and mass-dependent carrying capacity of basal prey #
  K_r = Car*MRr^Sk*exp(Eak*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  
  #### temperature- and mass-dependent metabolic rate of consumer#
  x_r= exp(Ix)*MCr^Sx*exp(Eax*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  x_i= exp(Ix)*MCi^Sx*exp(Eax*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  #### temperature- and mass-dependent handling time of consumer on resident or invading prey
  Th_b_r = exp(Iy) * MCr^Sy_pred * MRr^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  Th_b_i = exp(Iy) * MCi^Sy_pred * MRr^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent component of handling time (quadratic correction, eq 2.13 in Binzer et al. 2012)
  Th_m_r = exp(IThm + S1_Thm * log(MCr/MRr) + S2_Thm * (log(MCr/MRr))^2);
  Th_m_i = exp(IThm + S1_Thm * log(MCi/MRr) + S2_Thm * (log(MCi/MRr))^2);
  # temperature-dependent component of handling time (quadratic correction, eq 2.14 in Binzer et al. 2012)
  Th_T = exp(IThT + S1_ThT * (Temp) + S2_ThT * (Temp)^2);
  # combination of the different functions
  th_r = Th_b_r*Th_m_r*Th_T;
  th_i = Th_b_i*Th_m_i*Th_T;
  
  #### temperature- and mass-dependent attack rate of consumer on resident or invading prey
  # mass-dependent half-saturation constant (baseline allometry, eq 2.10 in Binzer et al. 2012)
  a_b_r = exp(I_B0) * MCr^(SB0_pred) * MRr^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  a_b_i = exp(I_B0) * MCi^(SB0_pred) * MRr^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent attack rate (quadratic correction, eq 2.12 in Binzer et al. 2012)
  a_m_r = exp(Iam + S1_am * log(MCr/MRr) + S2_am * (log(MCr/MRr))^2);
  a_m_i = exp(Iam + S1_am * log(MCi/MRr) + S2_am * (log(MCi/MRr))^2);
  # combination of the different functions
  a_r = a_b_r*a_m_r;
  a_i = a_b_i*a_m_i;
  
  # Resource equilirbria
  R_star_Cr = x_r/(a_r*(e - x_r*th_r));
  R_star_Ci = x_i/(a_i*(e - x_i*th_i));
  
  # Consumer equilibria
  Cr_star = (r_r/(K_r*a_r))*(1+ (a_r*th_r*R_star_Cr))*(K_r - R_star_Cr);
  Ci_star = (r_r/(K_r*a_i))*(1+ (a_i*th_i*R_star_Ci))*(K_r - R_star_Ci);
  
  res<-cbind(Temp, Car, MRr, MCi, MCr, R_star_Cr, R_star_Ci, Cr_star, Ci_star);
  #
  return(res)
}

#### P* rule function for AC module
P.Rule.Appcomp <- function(Temp, Car, MRr, MRi, MCr){
  #### temperature- and mass-dependent growth rate of basal prey #
  r_r= exp(Ir)*MRr^Sr*exp(Ear*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  r_i= exp(Ir)*MRi^Sr*exp(Ear*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  #### temperature- and mass-dependent carrying capacity of basal prey #
  K_r = Car*MRr^Sk*exp(Eak*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  K_i = Car*MRi^Sk*exp(Eak*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  
  #### temperature- and mass-dependent metabolic rate of consumer#
  x2= exp(Ix)*MCr^Sx*exp(Eax*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  
  #### temperature- and mass-dependent handling time of consumer on resident or invading prey
  Th_b_r = exp(Iy) * MCr^Sy_pred * MRr^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  Th_b_i = exp(Iy) * MCr^Sy_pred * MRi^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent component of handling time (quadratic correction, eq 2.13 in Binzer et al. 2012)
  Th_m_r = exp(IThm + S1_Thm * log(MCr/MRr) + S2_Thm * (log(MCr/MRr))^2);
  Th_m_i = exp(IThm + S1_Thm * log(MCr/MRi) + S2_Thm * (log(MCr/MRi))^2);
  # temperature-dependent component of handling time (quadratic correction, eq 2.14 in Binzer et al. 2012)
  Th_T = exp(IThT + S1_ThT * (Temp) + S2_ThT * (Temp)^2);
  # combination of the different functions
  th_r = Th_b_r*Th_m_r*Th_T;
  th_i = Th_b_i*Th_m_i*Th_T;
  
  #### temperature- and mass-dependent attack rate of consumer on resident or invading prey
  # mass-dependent half-saturation constant (baseline allometry, eq 2.10 in Binzer et al. 2012)
  a_b_r = exp(I_B0) * MCr^(SB0_pred) * MRr^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  a_b_i = exp(I_B0) * MCr^(SB0_pred) * MRi^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent attack rate (quadratic correction, eq 2.12 in Binzer et al. 2012)
  a_m_r = exp(Iam + S1_am * log(MCr/MRr) + S2_am * (log(MCr/MRr))^2);
  a_m_i = exp(Iam + S1_am * log(MCr/MRi) + S2_am * (log(MCr/MRi))^2);
  # combination of the different functions
  a_r = a_b_r*a_m_r;
  a_i = a_b_i*a_m_i;
  
  # Resource equilirbria
  Rr_star = x2/(a_r*(e - x2*th_r));
  Ri_star = x2/(a_i*(e - x2*th_i));
  
  # Consumer equilibria
  C_star_Rr = (r_r/(K_r*a_r))*(1+ (a_r*th_r*Rr_star))*(K_r - Rr_star);
  C_star_Ri = (r_i/(K_i*a_i))*(1+ (a_i*th_i*Ri_star))*(K_i - Ri_star);
  
  res<-cbind(Temp, Car, MRr, MRi, MCr, Rr_star, Ri_star, C_star_Rr, C_star_Ri);
  #
  return(res)
}

#### P* rule function
PR.Rules.IGP <- function(Temp, Car, MRr, MCr, MCi){
  #### temperature- and mass-dependent growth rate of basal prey #
  r_r= exp(Ir)*MRr^Sr*exp(Ear*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  #### temperature- and mass-dependent carrying capacity of basal prey #
  K_r = Car*MRr^Sk*exp(Eak*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  
  #### temperature- and mass-dependent metabolic rate of consumer#
  x_r= exp(Ix)*MCr^Sx*exp(Eax*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  x_i= exp(Ix)*MCi^Sx*exp(Eax*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  
  #### temperature- and mass-dependent handling time of consumer feeding on resident resource
  Th_b_Cr = exp(Iy) * MCr^Sy_pred * MRr^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  Th_b_Ci = exp(Iy) * MCi^Sy_pred * MRr^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent component of handling time (quadratic correction, eq 2.13 in Binzer et al. 2012)
  Th_m_Cr = exp(IThm + S1_Thm * log(MCr/MRr) + S2_Thm * (log(MCr/MRr))^2);
  Th_m_Ci = exp(IThm + S1_Thm * log(MCi/MRr) + S2_Thm * (log(MCi/MRr))^2);
  # temperature-dependent component of handling time (quadratic correction, eq 2.14 in Binzer et al. 2012)
  Th_T = exp(IThT + S1_ThT * (Temp) + S2_ThT * (Temp)^2);
  # combination of the different functions
  th_Cr_R = Th_b_Cr*Th_m_Cr*Th_T;
  th_Ci_R = Th_b_Ci*Th_m_Ci*Th_T;
  
  #### temperature- and mass-dependent handling time of IG-predator feeding on IG-prey
  Th_b_CrCi = exp(Iy) * MCr^Sy_pred * MCi^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  Th_b_CiCr = exp(Iy) * MCi^Sy_pred * MCr^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent component of handling time (quadratic correction, eq 2.13 in Binzer et al. 2012)
  Th_m_CrCi = exp(IThm + S1_Thm * log(MCr/MCi) + S2_Thm * (log(MCr/MCi))^2);
  Th_m_CiCr = exp(IThm + S1_Thm * log(MCi/MCr) + S2_Thm * (log(MCi/MCr))^2);
  # combination of the different functions
  th_Cr_Ci = Th_b_CrCi*Th_m_CrCi*Th_T;
  th_Ci_Cr = Th_b_CiCr*Th_m_CiCr*Th_T;
  
  #### temperature- and mass-dependent attack rate of consumer on resident resource
  # mass-dependent half-saturation constant (baseline allometry, eq 2.10 in Binzer et al. 2012)
  a_b_Cr = exp(I_B0) * MCr^(SB0_pred) * MRr^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  a_b_Ci = exp(I_B0) * MCi^(SB0_pred) * MRr^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent attack rate (quadratic correction, eq 2.12 in Binzer et al. 2012)
  a_m_Cr = exp(Iam + S1_am * log(MCr/MRr) + S2_am * (log(MCr/MRr))^2);
  a_m_Ci = exp(Iam + S1_am * log(MCi/MRr) + S2_am * (log(MCi/MRr))^2);
  # combination of the different functions
  a_Cr = a_b_Cr*a_m_Cr;
  a_Ci = a_b_Ci*a_m_Ci;
  
  #### temperature- and mass-dependent attack rate of IG-predator feeding on IG-prey
  # mass-dependent half-saturation constant (baseline allometry, eq 2.10 in Binzer et al. 2012)
  a_b_CrCi = exp(I_B0) * MCr^(SB0_pred) * MCi^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  a_b_CiCr = exp(I_B0) * MCi^(SB0_pred) * MCr^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent attack rate (quadratic correction, eq 2.12 in Binzer et al. 2012)
  a_m_CrCi = exp(Iam + S1_am * log(MCr/MCi) + S2_am * (log(MCr/MCi))^2);
  a_m_CiCr = exp(Iam + S1_am * log(MCi/MCr) + S2_am * (log(MCi/MCr))^2);
  # combination of the different functions
  a_CrCi = a_b_CrCi*a_m_CrCi;
  a_CiCr = a_b_CiCr*a_m_CiCr;
  
  # Resource equilirbrium in CR interaction
  R_star_Cr = x_r/(a_Cr*(e - x_r*th_Cr_R));
  R_star_Ci = x_i/(a_Ci*(e - x_i*th_Ci_R));
  
  # Consumer equilibrium in CR interaction
  Cr_star = (r_r/(K_r*a_Cr))*(1+ (a_Cr*th_Cr_R*R_star_Cr))*(K_r - R_star_Cr);
  Ci_star = (r_r/(K_r*a_Ci))*(1+ (a_Ci*th_Ci_R*R_star_Ci))*(K_r - R_star_Ci);
  
  # IG-predator equilibrium in PC interaction
  Pi_Ir = (1/a_CiCr)*( ((e*a_Cr*R_star_Cr)/(1+a_Cr*th_Cr_R*R_star_Cr))- x_r)*(1+a_CiCr*th_Ci_Cr*Cr_star);
  
  # Initial Per-capita growth of IG-prey 
  dIi <- (e*a_Ci*(R_star_Cr*1.02)/(1+a_Ci*th_Ci_R*(R_star_Cr*1.02))) - x_i*1e-6 - Cr_star*(e*a_CrCi*1e-6/(1+a_CrCi*th_Cr_Ci*1e-6))
  dIr <- (e*a_Cr*(R_star_Cr*1.02)/(1+a_Cr*th_Cr_R*(R_star_Cr*1.02))) - x_r*Cr_star*1.02 - 1e-6*(e*a_CiCr*Cr_star*1.02/(1+a_CiCr*th_Ci_Cr*Cr_star*1.02))
  
  # Initial Per-capita growth of IG-predator
  dPi <- (e*(((a_Ci*R_star_Cr*1.02)/(1+a_Ci*th_Ci_R*R_star_Cr*1.02))+((a_CiCr*Cr_star*1.02)/(1+a_CiCr*th_Ci_Cr*Cr_star*1.02))) - x_i )
  
  f_PiRr <- (a_Ci*R_star_Cr*1.02)/(1+a_Ci*th_Ci_R*R_star_Cr*1.02)
  f_PiCr <- (a_CiCr*Cr_star*1.02)/(1+a_CiCr*th_Ci_Cr*Cr_star*1.02)
  
  res<-cbind(Temp, Car, MRr, MCi, MCr, R_star_Cr, R_star_Ci, Cr_star, Ci_star, dIi, dIr, dPi, f_PiRr, f_PiCr);
  return(res);
}
######################## Functions for species equilibra and biomass at equilibrium ====
################ Equilibrium ====
########## Equilibrium for bitrophic system
##### Interaction between species A & B
EQUI_2spAB <- function(Temp,Car,Alpha){
  #### Species masses #
  Mx = 0.001; (My = Mx * Alpha);
  #### temperature- and mass-dependent growth rate of basal prey #
  r1= exp(Ir)*Mx^Sr*exp(Ear*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  #### temperature- and mass-dependent carrying capacity of basal prey #
  K1 = Car*Mx^Sk*exp(Eak*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  #### temperature- and mass-dependent metabolic rate of intermediate predator #
  x2= exp(Ix)*My^Sx*exp(Eax*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  
  #### temperature- and mass-dependent handling time  of intermediate predator that prey on basal species #
  Th_basalY = exp(Iy) * My^Sy_pred * Mx^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent component of handling time (quadratic correction, eq 2.13 in Binzer et al. 2012)
  Th_massY = exp(IThm + S1_Thm * log(My/Mx) + S2_Thm * (log(My/Mx))^2);
  # temperature-dependent component of handling time (quadratic correction, eq 2.14 in Binzer et al. 2012)
  Th_TempY = exp(IThT + S1_ThT * (Temp) + S2_ThT * (Temp)^2);
  # combination of the different functions
  th21 = Th_basalY*Th_massY*Th_TempY;
  
  #### temperature- and mass-dependent attack rate of intermediate predator that prey on basal species #
  # mass-dependent half-saturation constant (baseline allometry, eq 2.10 in Binzer et al. 2012)
  a_baseY = exp(I_B0) * My^(SB0_pred) * Mx^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent attack rate (quadratic correction, eq 2.12 in Binzer et al. 2012)
  a_massY = exp(Iam + S1_am * log(My/Mx) + S2_am * (log(My/Mx))^2);
  # combination of the different functions
  a21 = a_baseY*a_massY;
  
  #### isocline of basal species #
  Bx = x2/(a21*(e - x2*th21));
  
  #### isocline of Intermediate consumer #
  By = (r1/(K1*a21))*(1+ (a21*th21*Bx))*(K1 - Bx);   
  
  res<-cbind(Temp,Car,Bx,By,Mx,My,Alpha);
  return(res)
}
##### Interaction between species A & C
EQUI_2spAC <- function(Temp,Car,Gamma){
  #### Species masses
  Mx = 0.001; (Mz = Mx * Gamma);
  #### temperature- and mass-dependent growth rate of basal prey #
  r1= exp(Ir)*Mx^Sr*exp(Ear*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  #### temperature- and mass-dependent carrying capacity of basal prey #
  K1 = Car*Mx^Sk*exp(Eak*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  #### temperature- and mass-dependent metabolic rate of top predator #
  x3= exp(Ix)*Mz^Sx*exp(Eax*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  
  #### temperature- and mass-dependent handling time  of top predator that prey on basal species #
  Th_basalY = exp(Iy) * Mz^Sy_pred * Mx^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent component of handling time (quadratic correction, eq 2.13 in Binzer et al. 2012)
  Th_massY = exp(IThm + S1_Thm * log(Mz/Mx) + S2_Thm * (log(Mz/Mx))^2);
  # temperature-dependent component of handling time (quadratic correction, eq 2.14 in Binzer et al. 2012)
  Th_TempY = exp(IThT + S1_ThT * (Temp) + S2_ThT * (Temp)^2);
  # combination of the different functions
  th31 = Th_basalY*Th_massY*Th_TempY;
  
  ####temperature- and mass-dependent attack rate of top predator that prey on basal species #
  # mass-dependent half-saturation constant (baseline allometry, eq 2.10 in Binzer et al. 2012)
  a_baseY = exp(I_B0) * Mz^(SB0_pred) * Mx^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # size-dependent attack rate (quadratic correction, eq 2.12 in Binzer et al. 2012)
  a_massY = exp(Iam + S1_am * log(Mz/Mx) + S2_am * (log(Mz/Mx))^2);
  # combination of the different functions
  a31 = a_baseY*a_massY;
  
  #### isocline of basal species #
  Bx = x3/(a31*(e - x3*th31));
  
  #### isocline of top predator #
  Bz = (r1/(K1*a31))*(1+ (a31*th31*Bx))*(K1 - Bx);   
  
  res<-cbind(Temp,Car,Bx,Bz,Mx,Mz,Gamma);
  return(res);
}
##### Interaction between species B & C
EQUI_2spBC <- function(Temp,Car,Alpha,Beta){
  #### Species masses
  My = 0.001*Alpha; (Mz = My * Beta);
  #### temperature- and mass-dependent growth rate of basal prey #
  r2= exp(Ir)*My^Sr*exp(Ear*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  #### temperature- and mass-dependent carrying capacity of basal prey #
  K2 = Car*My^Sk*exp(Eak*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  #### temperature- and mass-dependent metabolic rate of top predator #
  x3= exp(Ix)*Mz^Sx*exp(Eax*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  
  #### temperature- and mass-dependent handling time  of top predator that prey on basal species #
  Th_basalY = exp(Iy) * Mz^Sy_pred * My^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent component of handling time (quadratic correction, eq 2.13 in Binzer et al. 2012)
  Th_massY = exp(IThm + S1_Thm * log(Mz/My) + S2_Thm * (log(Mz/My))^2);
  # temperature-dependent component of handling time (quadratic correction, eq 2.14 in Binzer et al. 2012)
  Th_TempY = exp(IThT + S1_ThT * (Temp) + S2_ThT * (Temp)^2);
  # combination of the different functions
  th32 = Th_basalY*Th_massY*Th_TempY;
  
  ####temperature- and mass-dependent attack rate of top predator that prey on basal species #
  # mass-dependent half-saturation constant (baseline allometry, eq 2.10 in Binzer et al. 2012)
  a_baseY = exp(I_B0) * Mz^(SB0_pred) * My^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # size-dependent attack rate (quadratic correction, eq 2.12 in Binzer et al. 2012)
  a_massY = exp(Iam + S1_am * log(Mz/My) + S2_am * (log(Mz/My))^2);
  # combination of the different functions
  a32 = a_baseY*a_massY;
  
  #### isocline of basal species #
  By = x3/(a32*(e - x3*th32));
  
  #### isocline of top predator #
  Bz = (r2/(K2*a32))*(1+ (a32*th32*By))*(K2 - By);   
  
  res<-cbind(Temp,Car,By,Bz,My,Mz,Alpha,Beta);
  return(res);
}
########## Equilibrium for tritrophic chain)
EQUI_3sp <- function(Temp,Car,Alpha,Beta,Gamma){
  #### Species masses
  Mx = 0.001; (My = Mx*Alpha); (Mz = My*Beta);
  #### temperature- and mass-dependent growth rate of basal prey #
  r= exp(Ir)*Mx^Sr*exp(Ear*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  #### temperature- and mass-dependent carrying capacity of basal prey #
  K = Car*Mx^Sk*exp(Eak*(T0-(Temp+T0K))/(Boltz*(Temp+T0K) * T0));
  #### temperature- and mass-dependent metabolic rate of intermediate predator #
  x2= exp(Ix)*My^Sx*exp(Eax*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  #### temperature- and mass-dependent metabolic rate of top predator #
  x3= exp(Ix)*Mz^Sx*exp(Eax*(T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  
  #### temperature- and mass-dependent handling time of intermediate predator that prey on basal species #
  Th_basalY = exp(Iy) * My^Sy_pred * Mx^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent component of handling time (quadratic correction, eq 2.13 in Binzer et al. 2012)
  Th_massY = exp(IThm + S1_Thm * log(My/Mx) + S2_Thm * (log(My/Mx))^2);
  # temperature-dependent component of handling time (quadratic correction, eq 2.14 in Binzer et al. 2012)
  Th_TempY = exp(IThT + S1_ThT * (Temp) + S2_ThT * (Temp)^2);
  # combination of the different functions
  th21 = Th_basalY*Th_massY*Th_TempY;
  
  #### temperature- and mass-dependent handling time of top predator on the intermediary predator #
  Th_baseZy= exp(Iy) * Mz^Sy_pred * My^Sy_prey * exp(Eay * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  Th_massZy = exp(IThm + S1_Thm * log(Mz/My) + S2_Thm * (log(Mz/My))^2);
  Th_TempZy = exp(IThT + S1_ThT * (Temp) + S2_ThT * (Temp)^2);
  th32 = Th_baseZy*Th_massZy*Th_TempZy;
  
  #### temperature- and mass-dependent attack rate of intermediate predator that prey on basal species #
  # mass-dependent half-saturation constant (baseline allometry, eq 2.10 in Binzer et al. 2012)
  a_baseY = exp(I_B0) * My^(SB0_pred) * Mx^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  # mass-dependent attack rate (quadratic correction, eq 2.12 in Binzer et al. 2012)
  a_massY = exp(Iam + S1_am * log(My/Mx) + S2_am * (log(My/Mx))^2);
  # combination of the different functions
  a21 = a_baseY*a_massY;
  
  #### temperature- and mass-dependent attack rate of top predator that prey on intermediary predator #
  a_baseZy = exp(I_B0) * Mz^(SB0_pred) * My^(SB0_prey) * exp(EaB0 * (T0-(Temp+T0K))/(Boltz*(Temp+T0K)*T0));
  a_massZy = exp(Iam + S1_am * log(Mz/My) + S2_am * (log(Mz/My))^2);
  a32 = a_baseZy*a_massZy;
  
  #### Determinant of the polynomial form of the solution when solving basal species isocline #
  Det= (K*r*a21*th21-r)^2+4*K*r*a21*th21*(r-a21*(x3/(a32*(e-x3*th32)))); 
  
  #### isocline of basal species #
  Bx = (r*(K*a21*th21-1) + sqrt(Det))/(2*r*a21*th21);
  
  #### isocline of intermediate consumer #
  By = x3/(a32*(e-th32*x3));
  
  #### isocline of top predator #
  Num1 = (1+a32*th32*(x3/(a32*(e-th32*x3))));
  Num2 = ( (e*a21-x2*a21*th21)*((r*(K*a21*th21-1) + sqrt(Det))/(2*r*a21*th21))-x2);
  Denom = (a32*(1+a21*th21*((r*(K*a21*th21-1) + sqrt(Det))/(2*r*a21*th21) )));
  Bz = (Num1*Num2)/Denom;
  
  res<-cbind(Temp,Car,Bx,By,Bz,Mx,My,Mz,Alpha,Beta,Gamma);
  return(res);
}

################ Functions for species biomass at equilibrium ====
########## Densities for 2 species equilibrium  ====
#### Interaction between species A & B
EQ_2spAB<-function(T_min,T_max,Alpha){
  
  Result<-data.frame(Temp=numeric(0), Car=numeric(0),
                     Bx=numeric(0), By=numeric(0), Mx=numeric(0), My=numeric(0), Alpha=numeric(0));
  
  for (i in seq(T_min,T_max, by=0.1)){
    for (j in seq(0.1,20, by=0.1)){
      Tempfile <- as.data.frame(EQUI_2spAB(i, j, Alpha));
      Result <- rbind(Result, Tempfile);
    }
  }
  res<-data.frame(Result);
  return(res);
}
#### Interaction between species A & C
EQ_2spAC<-function(T_min,T_max,Gamma){
  
  Result<-data.frame(Temp=numeric(0), Car=numeric(0),
                     Bx=numeric(0), Bz=numeric(0), Mx=numeric(0), Mz=numeric(0), Gamma=numeric(0));
  
  for (i in seq(T_min,T_max, by=0.1)){
    for (j in seq(0.1,20, by=0.1)){
      Tempfile <- as.data.frame(EQUI_2spAC(i, j, Gamma));
      Result <- rbind(Result, Tempfile);
    }
  }
  res<-data.frame(Result);
  return(res);
}
#### Interaction between species B & C
EQ_2spBC<-function(T_min,T_max,Alpha,Beta){
  
  Result<-data.frame(Temp=numeric(0), Car=numeric(0),
                     By=numeric(0), Bz=numeric(0), My=numeric(0), Mz=numeric(0),
                     Alpha=numeric(0), Beta=numeric(0));
  
  for (i in seq(T_min,T_max, by=0.1)){
    for (j in seq(0.1,20, by=0.1)){
      Tempfile <- as.data.frame(EQUI_2spBC(i, j, Alpha, Beta));
      Result <- rbind(Result, Tempfile);
    }
  }
  res<-data.frame(Result);
  return(res);
}

########## Densities for 3 species equilibrium  ====
EQ_3sp<-function(T_min,T_max,Alpha,Beta,Gamma){
  
  Result<-data.frame(Temp=numeric(0), Car=numeric(0),
                     Bx=numeric(0), By=numeric(0), Bz=numeric(0),
                     Mx=numeric(0), My=numeric(0), Mz=numeric(0),
                     Alpha=numeric(0), Beta=numeric(0), Gamma=numeric(0));
  
  for (i in seq(T_min,T_max, by=0.1)){
    for (j in seq(0.1,20, by=0.1)){
      Tempfile <- as.data.frame(EQUI_3sp(i, j, Alpha,Beta, Gamma));
      Result <- rbind(Result, Tempfile);
    }
  }
  res<-data.frame(Result);
  return(res);
}

######################## Functions for each ODE system: ====
#### Control (2sp & 3sp) without species invasion: ABtc, ACtc, TrC
#### System with species invasion in bitrophic system:
#### TC (top predator), ExpComp (competitor), AppComp (apparent competitor),
#### & IGP (IGPB & IGPC, invasion of IG prey or IG predator

################ 2 species systems  ====
########## Interaction A-B system
TC_2spAB <- function(Time,Temp,Car,Bx,By,Mx,My,Alpha){
  ODEsyst <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fB <- aba*A/(1+aba*A*thba);      # Type II Functional response of B on prey A
      dA <- ra*A*(1-A/Ka)-B*fB;        # Biomass changes of basal species A
      dB <- e*B*fB-xb*B;               # Biomass changes of top species B
      list(c(dA,dB));
    })
  }
  
  ## event triggered if state variable = Threshold value
  Tr<-10^(-12)                         # Threshold value for extinction
  rootfun <- function(t, state, pars) {
    return(state - Tr);
  }
  
  ## sets state variable = 0
  eventfun <- function(t, state, pars) {
    which ((state-Tr) <= 0);
    STATE <- state;
    STATE[which ((state-Tr) <=0)] <- 0;
    return(STATE);
  }
  
  parameters <- c(
    ra  <- r(Temp, Mx),                # Growth rate of A
    Ka  <- K(Temp, Mx, Car),           # Carrying capacity of A
    xb  <- x2(Temp, My),               # Metabolism of B
    thba <- th(Temp, My, Mx),          # Handling time of B on A
    aba  <- a(Temp, My, Mx)            # Attack rate of B on A
  )
  
  Y = Time;                               # Number of year for the simulation
  Yearlength=Y*31557600;
  times <- seq(0,Yearlength, by= 86400);  # Time period in seconds converted in Y year per days 
  state <- c(A=Bx, B=By);                 # Initial biomasses
  
  out <- as.data.frame(ode(times= times, func = ODEsyst, y = state, parms = parameters,
                           events = list(func = eventfun, root = TRUE), rootfun = rootfun));
  
  A_ext_time <- ifelse(length(out[which(out[,2] == 0), 1])>0, min(out[which(out[,2] == 0), 1]), NA)
  B_ext_time <- ifelse(length(out[which(out[,3] == 0), 1])>0, min(out[which(out[,3] == 0), 1]), NA)
  Time_ext <- c(A_ext_time, B_ext_time)
  res<-append(tail(out,1), Time_ext);
  
  ## The last outcomes are used for final biomasses
  out_tail <- tail(out[,2:3],3653);             # take the last 10 years from out data
  Di <- seq(1, 3653, by=365.25);                # First day of the last years
  Df <- Di-1; Df <- Df[-1]; Df <- c(Df,3653);   # Last day of the last years
  Di <- ceiling(Di); Df <- ceiling(Df);
  
  pop <- c(); i <- 0;
  for(i in 1:length(Di)){
    j <- 0;
    out_sub <- out_tail[Di[i]:Df[i],];
    for(j in 1:2){
      min <- min(out_sub[,j]);   # Calculate min density of species i
      max <- max(out_sub[,j]);   # Calculate max density of species i
      mean <- mean(out_sub[,j]); # Calculate mean density of species i 
      sd <- sd(out_sub[,j]);     # Calculate sd of density of species i
      cv  <- sd/mean;            # Calculate the coeff. of variation of species density i
      measures <- c(min, max, mean, sd, cv);
      pop <- append(pop, measures);
      j<-j+1;
    }
    i <- i+1;
  }
  res<-append(res,pop);
}
########## Interaction A-C system
TC_2spAC <- function(Time,Temp,Car,Bx,Bz,Mx,Mz,Gamma){
  ODEsyst <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fC <- aca*A/(1+aca*A*thca);       # Type II Functional response of B on prey A
      dA <- ra*A*(1-A/Ka)-C*fC;         # Biomass changes of basal species A
      dC <- e*C*fC-xc*C;                # Biomass changes of top species C
      list(c(dA,dC));
    })
  }
  
  ## event triggered if state variable = Threshold value
  Tr<-10^(-12)                          # Threshold value for extinction
  rootfun <- function(t, state, pars) {
    return(state - Tr);
  }
  
  ## sets state variable = 0
  eventfun <- function(t, state, pars) {
    which ((state-Tr) <= 0);
    STATE <- state;
    STATE[which ((state-Tr) <=0)] <- 0;
    return(STATE);
  }
  
  parameters <- c(
    ra   <- r(Temp, Mx),                 # Growth rate of A
    Ka   <- K(Temp, Mx, Car),            # Carrying capacity of A
    xc  <- x2(Temp, Mz),                # Metabolism of C
    thca <- th(Temp, Mz, Mx),            # Handling time of C on A
    aca  <- a(Temp, Mz, Mx)              # Attack rate of C on A
  );
  
  Y = Time;                              # Number of year for the simulation
  Yearlength=Y*31557600;
  times <- seq(0,Yearlength, by= 86400); # Time period in seconds converted in Y year per days 
  state <- c(A=Bx, C=Bz);                # Initial biomasses
  
  out <- as.data.frame(ode(times= times, func = ODEsyst, y = state, parms = parameters,
                           events = list(func = eventfun, root = TRUE), rootfun = rootfun));
  
  A_ext_time <- ifelse(length(out[which(out[,2] == 0), 1])>0, min(out[which(out[,2] == 0), 1]), NA)
  C_ext_time <- ifelse(length(out[which(out[,3] == 0), 1])>0, min(out[which(out[,3] == 0), 1]), NA)
  Time_ext <- c(A_ext_time, C_ext_time)
  res<-append(tail(out,1), Time_ext);
  
  ## The last outcomes are used for final biomasses
  out_tail <- tail(out[,2:3],3653);             # take the last 10 years from out data
  Di <- seq(1, 3653, by=365.25);                # First day of the last years
  Df <- Di-1; Df <- Df[-1]; Df <- c(Df,3653);   # Last day of the last years
  Di <- ceiling(Di); Df <- ceiling(Df);
  
  pop <- c(); i <- 0;
  for(i in 1:length(Di)){
    j <- 0;
    out_sub <- out_tail[Di[i]:Df[i],];
    for(j in 1:2){
      min <- min(out_sub[,j]);   # Calculate min density of species i
      max <- max(out_sub[,j]);   # Calculate max density of species i
      mean <- mean(out_sub[,j]); # Calculate mean density of species i 
      sd <- sd(out_sub[,j]);     # Calculate sd of density of species i
      cv  <- sd/mean;            # Calculate the coeff. of variation of species density i
      measures <- c(min, max, mean, sd, cv);
      pop <- append(pop, measures);
      j<-j+1;
    }
    i <- i+1;
  }
  res<-append(res,pop);
}

################ 3 species systems  ====
########## Trophic chain system
TC <- function(Time,Temp,Car,Bx,By,Bz,Mx,My,Mz,Alpha,Beta,Gamma){
  ## ODE system
  ODEsyst <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fB <- aba*A/(1+aba*A*thba);      # Type II Functional response of B on prey A
      fC <- acb*B/(1+acb*B*thcb);      # Type II Functional response of C on prey B
      dA <- ra*A*(1-A/Ka)-B*fB;        # Biomass changes of basal species A
      dB <- e*B*fB-xb*B-C*fC;          # Biomass changes of top species B
      dC <- e*C*fC-xc*C;               # Biomass changes of top species C
      list(c(dA,dB,dC));
    })
  }
  
  ## event triggered if state variable = Threshold value
  Tr<-10^(-12)                         # Threshold value for extinction
  rootfun <- function(t, state, pars) {
    return(state - Tr);
  }
  
  ## sets state variable = 0
  eventfun <- function(t, state, pars) {
    which ((state-Tr) <= 0);
    STATE <- state;
    STATE[which ((state-Tr) <=0)] <- 0;
    return(STATE);
  }
  
  parameters <- c(
    ra   <- r(Temp, Mx),               # Growth rate of A
    Ka   <- K(Temp, Mx, Car),          # Carrying capacity of A
    xb  <- x2(Temp, My),               # Metabolism of B
    xc  <- x2(Temp, Mz),               # Metabolism of C
    thba <- th(Temp, My, Mx),          # Handling time of B on A
    thcb <- th(Temp, Mz, My),          # Handling time of C on A
    aba  <- a(Temp, My, Mx),           # Attack rate of B on A
    acb  <- a(Temp, Mz, My)            # Attack rate of C on B
  );
  
  Y = Time;                            # Number of year for the simulation
  Yearlength=Y*31557600;
  times <- seq(0,Yearlength, by= 86400);   # Time period in seconds converted in Y year per days 
  state <- c(A=Bx, B=By, C=Bz);            # Initial biomasses
  
  out <- as.data.frame(ode(times= times, func = ODEsyst, y = state, parms = parameters,
                           events = list(func = eventfun, root = TRUE), rootfun = rootfun));
  
  A_ext_time <- ifelse(length(out[which(out[,2] == 0), 1])>0, min(out[which(out[,2] == 0), 1]), NA)
  B_ext_time <- ifelse(length(out[which(out[,3] == 0), 1])>0, min(out[which(out[,3] == 0), 1]), NA)
  C_ext_time <- ifelse(length(out[which(out[,4] == 0), 1])>0, min(out[which(out[,4] == 0), 1]), NA)
  Time_ext <- c(A_ext_time, B_ext_time, C_ext_time)
  res<-append(tail(out,1), Time_ext);
  
  ## The last outcomes are used for final biomasses
  out_tail <- tail(out[,2:4],3653);             # take the last 10 years from out data
  Di <- seq(1, 3653, by=365.25);                # First day of the last years
  Df <- Di-1; Df <- Df[-1]; Df <- c(Df,3653);   # Last day of the last years
  Di <- ceiling(Di); Df <- ceiling(Df);
  
  pop <- c(); i <- 0;
  for(i in 1:length(Di)){
    j <- 0;
    out_sub <- out_tail[Di[i]:Df[i],];
    for(j in 1:3){
      min <- min(out_sub[,j]);   # Calculate min density of species i
      max <- max(out_sub[,j]);   # Calculate max density of species i
      mean <- mean(out_sub[,j]); # Calculate mean density of species i 
      sd <- sd(out_sub[,j]);     # Calculate sd of density of species i
      cv  <- sd/mean;            # Calculate the coeff. of variation of species density i
      measures <- c(min, max, mean, sd, cv);
      pop <- append(pop, measures);
      j<-j+1;
    }
    i <- i+1;
  }
  res<-append(res,pop);
}
########## Exploitative Competition
ExpComp<-function(Time,Temp,Car,Bx,By,Bz,Mx,My,Mz,Alpha,Beta,Gamma){
  ODEsyst <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fB <- aba*A/(1+aba*A*thba);      # Type II Functional response of B on prey A
      fC <- aca*A/(1+aca*A*thca);      # Type II Functional response of C on prey A
      dA <- ra*A*(1-A/Ka)-(B*fB+C*fC); # Biomass changes of basal species A
      dB <- e*B*fB-xb*B;               # Biomass changes of top species B
      dC <- e*C*fC-xc*C;               # Biomass changes of top species C
      list(c(dA,dB,dC));
    })
  }
  
  ## event triggered if state variable = Threshold value
  Tr<-10^(-12)                         # Threshold value for extinction
  rootfun <- function(t, state, pars) {
    return(state - Tr);
  }
  
  ## sets state variable = 0
  eventfun <- function(t, state, pars) {
    which ((state-Tr) <= 0);
    STATE <- state;
    STATE[which ((state-Tr) <=0)] <- 0;
    return(STATE);
  }
  
  parameters <- c(
    ra  <- r(Temp,Mx),                  # Growth rate of A
    Ka  <- K(Temp,Mx,Car),              # Carrying capacity of A
    xb <- x2(Temp,My),                  # Metabolism of B
    xc <- x2(Temp,Mz),                  # Metabolism of C
    thba <- th(Temp,My,Mx),             # Handling time of B on A
    thca <- th(Temp,Mz,My),             # Handling time of C on A
    aba  <- a(Temp,My,Mx),              # Attack rate of B on A
    aca  <- a(Temp,Mz,Mx)               # Attack rate of C on B
  );
  
  Y = Time;                             # Number of year for the simulation
  Yearlength=Y*31557600;
  times <- seq(0,Yearlength, by= 86400);  # Time period in seconds converted in Y year per days 
  state <- c(A=Bx, B=By, C=Bz);           # Starting biomasses
  
  out <- as.data.frame(ode(times= times, func = ODEsyst, y = state, parms = parameters,
                           events = list(func = eventfun, root = TRUE), rootfun = rootfun));
  
  A_ext_time <- ifelse(length(out[which(out[,2] == 0), 1])>0, min(out[which(out[,2] == 0), 1]), NA)
  B_ext_time <- ifelse(length(out[which(out[,3] == 0), 1])>0, min(out[which(out[,3] == 0), 1]), NA)
  C_ext_time <- ifelse(length(out[which(out[,4] == 0), 1])>0, min(out[which(out[,4] == 0), 1]), NA)
  Time_ext <- c(A_ext_time, B_ext_time, C_ext_time)
  res<-append(tail(out,1), Time_ext);
  
  ## The last outcomes are used for final biomasses
  out_tail <- tail(out[,2:4],3653);             # take the last 10 years from out data
  Di <- seq(1, 3653, by=365.25);                # First day of the last years
  Df <- Di-1; Df <- Df[-1]; Df <- c(Df,3653);   # Last day of the last years
  Di <- ceiling(Di); Df <- ceiling(Df);
  
  pop <- c(); i <- 0;
  for(i in 1:length(Di)){
    j <- 0;
    out_sub <- out_tail[Di[i]:Df[i],];
    for(j in 1:3){
      min <- min(out_sub[,j]);   # Calculate min density of species i
      max <- max(out_sub[,j]);   # Calculate max density of species i
      mean <- mean(out_sub[,j]); # Calculate mean density of species i 
      sd <- sd(out_sub[,j]);     # Calculate sd of density of species i
      cv  <- sd/mean;            # Calculate the coeff. of variation of species density i
      measures <- c(min, max, mean, sd, cv);
      pop <- append(pop, measures);
      j<-j+1;
    }
    i <- i+1;
  }
  res<-append(res,pop);
}
########## Apparent Competition
AppComp<-function(Time,Temp,Car,Bx,By,Bz,Mx,My,Mz,Alpha,Beta,Gamma){
  ODEsyst <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fCA <- aca*A/(1+aca*A*thca);      # Type II Functional response of C on prey A
      fCB <- acb*B/(1+acb*B*thcb);      # Type II Functional response of C on prey B
      dA <- ra*A*(1-A/Ka)-C*fCA;        # Biomass changes of basal species A
      dB <- rb*B*(1-B/Kb)-C*fCB;        # Biomass changes of top species B
      dC <- e*C*(fCA+fCB)-xc*C;         # Biomass changes of top species C
      list(c(dA,dB,dC));
    })
  }
  
  ## event triggered if state variable = Threshold value
  Tr<-10^(-12)                          # Threshold value for extinction
  rootfun <- function(t, state, pars) {
    return(state - Tr);
  }
  
  ## sets state variable = 0
  eventfun <- function(t, state, pars) {
    which ((state-Tr) <= 0);
    STATE <- state;
    STATE[which ((state-Tr) <=0)] <- 0;
    return(STATE);
  }
  
  parameters <- c(
    ra  <- r(Temp, Mx),                  # Growth rate of A
    Ka  <- K(Temp, Mx, Car),             # Carrying capacity of A
    rb  <- r(Temp, My),                  # Growth rate of B
    Kb  <- K(Temp, My, Car),             # Carrying capacity of B
    xc  <- x2(Temp, Mz),                 # Metabolism of C
    thca <- th(Temp, Mz, Mx),            # Handling time of C on A
    thcb <- th(Temp, Mz, My),            # Handling time of C on B
    aca  <- a(Temp, Mz, Mx),             # Attack rate of C on A
    acb  <- a(Temp, Mz, My)              # Attack rate of C on B
  );
  
  Y = Time;                              # Number of year for the simulation
  Yearlength=Y*31557600;
  times <- seq(0,Yearlength, by= 86400); # Time period in seconds converted in Y year per days 
  state <- c(A=Bx, B=By, C=Bz);          # Starting biomasses
  
  out <- as.data.frame(ode(times= times, func = ODEsyst, y = state, parms = parameters,
                           events = list(func = eventfun, root = TRUE), rootfun = rootfun));
  
  A_ext_time <- ifelse(length(out[which(out[,2] == 0), 1])>0, min(out[which(out[,2] == 0), 1]), NA)
  B_ext_time <- ifelse(length(out[which(out[,3] == 0), 1])>0, min(out[which(out[,3] == 0), 1]), NA)
  C_ext_time <- ifelse(length(out[which(out[,4] == 0), 1])>0, min(out[which(out[,4] == 0), 1]), NA)
  Time_ext <- c(A_ext_time, B_ext_time, C_ext_time)
  res<-append(tail(out,1), Time_ext);
  
  ## The last outcomes are used for final biomasses
  out_tail <- tail(out[,2:4],3653);             # take the last 10 years from out data
  Di <- seq(1, 3653, by=365.25);                # First day of the last years
  Df <- Di-1; Df <- Df[-1]; Df <- c(Df,3653);   # Last day of the last years
  Di <- ceiling(Di); Df <- ceiling(Df);
  
  pop <- c(); i <- 0;
  for(i in 1:length(Di)){
    j <- 0;
    out_sub <- out_tail[Di[i]:Df[i],];
    for(j in 1:3){
      min <- min(out_sub[,j]);   # Calculate min density of species i
      max <- max(out_sub[,j]);   # Calculate max density of species i
      mean <- mean(out_sub[,j]); # Calculate mean density of species i 
      sd <- sd(out_sub[,j]);     # Calculate sd of density of species i
      cv  <- sd/mean;            # Calculate the coeff. of variation of species density i
      measures <- c(min, max, mean, sd, cv);
      pop <- append(pop, measures);
      j<-j+1;
    }
    i <- i+1;
  }
  res<-append(res,pop);
}
########## Intraguild Predation
IGP<-function(Time,Temp,Car,Bx,By,Bz,Mx,My,Mz,Alpha,Beta,Gamma){
  ODEsyst <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fBA <- aba*A/(1+aba*A*thba);         # Type II Functional response of B on prey A  
      fCA <- aca*A/(1+aca*A*thca);         # Type II Functional response of C on prey A 
      fCB <- acb*B/(1+acb*B*thcb);         # Type II Functional response of C on prey B
      dA <- ra*A*(1-A/Ka)-(C*fCA+B*fBA);   # Biomass changes of basal species A
      dB <- e*B*(fBA-xb)-C*fCB;            # Biomass changes of top species B
      dC <- e*C*(fCA+fCB)-xc*C;            # Biomass changes of top species C
      list(c(dA,dB,dC));
    })
  }
  
  ## event triggered if state variable = Threshold value
  Tr<-10^(-12)                         # Threshold value for extinction
  rootfun <- function(t, state, pars) {
    return(state - Tr);
  }
  
  ## sets state variable = 0
  eventfun <- function(t, state, pars) {
    which ((state-Tr) <= 0);
    STATE <- state;
    STATE[which ((state-Tr) <=0)] <- 0;
    return(STATE);
  }
  
  parameters <- c(
    ra  <- r(Temp, Mx),                # Growth rate of A
    Ka  <- K(Temp, Mx, Car),           # Carrying capacity of A
    xb  <- x2(Temp, My),               # Metabolism of B
    xc  <- x2(Temp, Mz),               # Metabolism of C
    thba <- th(Temp, My, Mx),          # Handling time of B on A
    thca <- th(Temp, Mz, Mx),          # Handling time of C on A
    thcb <- th(Temp, Mz, My),          # Handling time of C on B
    aca  <- a(Temp, Mz, Mx),           # Attack rate of C on A
    acb  <- a(Temp, Mz, My),           # Attack rate of C on B
    aba  <- a(Temp, My, Mx)            # Attack rate of B on A
  );
  
  Y = Time;                            # Number of year for the simulation
  Yearlength=Y*31557600;
  times <- seq(0,Yearlength, by= 86400);  # Time period in seconds converted in Y year per days 
  state <- c(A=Bx, B=By, C=Bz);           # Starting biomasses
  
  out <- as.data.frame(ode(times= times, func = ODEsyst, y = state, parms = parameters,
                           events = list(func = eventfun, root = TRUE), rootfun = rootfun));
  
  A_ext_time <- ifelse(length(out[which(out[,2] == 0), 1])>0, min(out[which(out[,2] == 0), 1]), NA)
  B_ext_time <- ifelse(length(out[which(out[,3] == 0), 1])>0, min(out[which(out[,3] == 0), 1]), NA)
  C_ext_time <- ifelse(length(out[which(out[,4] == 0), 1])>0, min(out[which(out[,4] == 0), 1]), NA)
  Time_ext <- c(A_ext_time, B_ext_time, C_ext_time)
  res<-append(tail(out,1), Time_ext);
  
  ## The last outcomes are used for final biomasses
  out_tail <- tail(out[,2:4],3653);             # take the last 10 years from out data
  Di <- seq(1, 3653, by=365.25);                # First day of the last years
  Df <- Di-1; Df <- Df[-1]; Df <- c(Df,3653);   # Last day of the last years
  Di <- ceiling(Di); Df <- ceiling(Df);
  
  pop <- c(); i <- 0;
  for(i in 1:length(Di)){
    j <- 0;
    out_sub <- out_tail[Di[i]:Df[i],];
    for(j in 1:3){
      min <- min(out_sub[,j]);   # Calculate min density of species i
      max <- max(out_sub[,j]);   # Calculate max density of species i
      mean <- mean(out_sub[,j]); # Calculate mean density of species i 
      sd <- sd(out_sub[,j]);     # Calculate sd of density of species i
      cv  <- sd/mean;            # Calculate the coeff. of variation of species density i
      measures <- c(min, max, mean, sd, cv);
      pop <- append(pop, measures);
      j<-j+1;
    }
    i <- i+1;
  }
  res<-append(res,pop);
}

######################## Jacobian functions: Calculus of dominant eigenvalues in steady state  ====
### Consumer-Resource B-A (Cr-Rr)
lambda_2sp_AB <- function(Temp, Car, Bx, By, Mx, My, Alpha) {
  parameters <- c(
    ra  <- r(Temp, Mx),                # Growth rate of A
    Ka  <- K(Temp, Mx, Car),           # Carrying capacity of A
    xb  <- x2(Temp, My),               # Metabolism of B
    thba <- th(Temp, My, Mx),          # Handling time of B on A
    aba  <- a(Temp, My, Mx)            # Attack rate of B on A
  )
  ODEsyst_2sp_AB <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fB <- aba*A/(1+aba*A*thba);      # Type II Functional response of B on prey A
      dA <- ra*A*(1-A/Ka)-B*fB;        # Biomass changes of basal species A
      dB <- e*B*fB-xb*B;               # Biomass changes of top species B
      list(c(dA,dB));
    })
  }
  EIG      <- eigen(jacobian.full(y = c(A=Bx, B=By), func = ODEsyst_2sp_AB, parms = as.list(parameters)))$values
  return(max(Re(EIG)))
}
### Consumer-Resource C-A (Cr-Rr)
lambda_2sp_AC <- function(Temp, Car, Bx, Bz, Mx, Mz, Gamma) {
  parameters <- c(
    ra   <- r(Temp, Mx),                 # Growth rate of A
    Ka   <- K(Temp, Mx, Car),            # Carrying capacity of A
    xc  <- x2(Temp, Mz),                # Metabolism of C
    thca <- th(Temp, Mz, Mx),            # Handling time of C on A
    aca  <- a(Temp, Mz, Mx)              # Attack rate of C on A
  )
  ODEsyst_2sp_AC <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fC <- aca*A/(1+aca*A*thca);       # Type II Functional response of B on prey A
      dA <- ra*A*(1-A/Ka)-C*fC;         # Biomass changes of basal species A
      dC <- e*C*fC-xc*C;                # Biomass changes of top species C
      list(c(dA,dC));
    })
  }
  EIG      <- eigen(jacobian.full(y = c(A=Bx, C=Bz), func = ODEsyst_2sp_AC, parms = as.list(parameters)))$values
  return(max(Re(EIG)))
}
### Consumer-Resource C-B (Cr-Ri)
lambda_2sp_BC <- function(Temp, Car, By, Bz, My, Mz, Beta) {
  parameters <- c(
    rb  <- r(Temp, My),                # Growth rate of B
    Kb  <- K(Temp, My, Car),           # Carrying capacity of B
    xc  <- x2(Temp, Mz),               # Metabolism of C
    thcb <- th(Temp, Mz, My),          # Handling time of C on B
    acb  <- a(Temp, Mz, My)            # Attack rate of C on B
  )
  ODEsyst_2sp_BC <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fC <- acb*B/(1+acb*B*thcb);      # Type II Functional response of C on prey B
      dB <- rb*B*(1-B/Kb)-C*fC;        # Biomass changes of basal species B
      dC <- e*C*fC-xc*C;               # Biomass changes of top species C
      list(c(dB,dC));
    })
  }
  EIG      <- eigen(jacobian.full(y = c(B=By, C=Bz), func = ODEsyst_2sp_BC, parms = as.list(parameters)))$values
  return(max(Re(EIG)))
}
### 3 sp system, Apparent competition
lambda_3sp_AppComp <- function(Temp, Car, Bx, By, Bz, Mx, My, Mz, Alpha, Beta, Gamma) {
  parameters <- c(
    ra  <- r(Temp, Mx),                # Growth rate of A
    Ka  <- K(Temp, Mx, Car),           # Carrying capacity of A
    rb  <- r(Temp, My),                # Growth rate of B
    Kb  <- K(Temp, My, Car),           # Carrying capacity of B
    xc  <- x2(Temp, Mz),               # Metabolism of C
    thca <- th(Temp, Mz, Mx),          # Handling time of C on A
    aca  <- a(Temp, Mz, Mx),           # Attack rate of C on A
    thcb <- th(Temp, Mz, My),          # Handling time of C on B
    acb  <- a(Temp, Mz, My)            # Attack rate of C on B
  )
  ODEsyst_3sp_AppComp <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fCA <- aca*A/(1+aca*A*thca);      # Type II Functional response of C on prey A
      fCB <- acb*B/(1+acb*B*thcb);      # Type II Functional response of C on prey B
      dA <- ra*A*(1-A/Ka)-C*fCA;        # Biomass changes of resident basal species A
      dB <- rb*B*(1-B/Kb)-C*fCB;        # Biomass changes of inding basal species B
      dC <- e*C*(fCA+fCB)-xc*C;         # Biomass changes of consumer species C
      list(c(dA,dB,dC));
    })
  }
  EIG      <- eigen(jacobian.full(y = c(A=Bx, B=By, C=Bz), func = ODEsyst_3sp_AppComp, parms = as.list(parameters)))$values
  return(max(Re(EIG)))
}
### 3 sp system, Exploitative competition
lambda_3sp_ExpComp <- function(Temp, Car, Bx, By, Bz, Mx, My, Mz, Alpha, Beta, Gamma) {
  parameters <- c(
    ra  <- r(Temp, Mx),                # Growth rate of A
    Ka  <- K(Temp, Mx, Car),           # Carrying capacity of A
    xb  <- x2(Temp, My),               # Metabolism of B
    thba <- th(Temp, My, Mx),          # Handling time of B on A
    aba  <- a(Temp, My, Mx),           # Attack rate of B on A
    xc  <- x2(Temp, Mz),               # Metabolism of C
    thca <- th(Temp, Mz, Mx),          # Handling time of C on A
    aca  <- a(Temp, Mz, Mx)            # Attack rate of C on A
  )
  ODEsyst_3sp_ExpComp <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fC <- aca*A/(1+aca*A*thca);      # Type II Functional response of C on prey A
      fB <- aba*A/(1+aba*A*thba);      # Type II Functional response of B on prey A
      dA <- ra*A*(1-A/Ka)-(B*fB+C*fC); # Biomass changes of basal species A
      dB <- e*B*fB-xb*B;               # Biomass changes of invading consumer species B
      dC <- e*C*fC-xc*C;               # Biomass changes of resident consumer species C
      list(c(dA,dB,dC));
    })
  }
  EIG      <- eigen(jacobian.full(y = c(A=Bx, B=By, C=Bz), func = ODEsyst_3sp_ExpComp, parms = as.list(parameters)))$values
  return(max(Re(EIG)))
}
### 3 sp system, Trophic competition
lambda_3sp_TC <- function(Temp, Car, Bx, By, Bz, Mx, My, Mz, Alpha, Beta, Gamma) {
  parameters <- c(
    ra  <- r(Temp, Mx),                # Growth rate of A
    Ka  <- K(Temp, Mx, Car),           # Carrying capacity of A
    xb  <- x2(Temp, My),               # Metabolism of B
    thba <- th(Temp, My, Mx),          # Handling time of B on A
    aba  <- a(Temp, My, Mx),           # Attack rate of B on A
    xc  <- x2(Temp, Mz),               # Metabolism of C
    thcb <- th(Temp, Mz, My),          # Handling time of C on B
    acb  <- a(Temp, Mz, My)            # Attack rate of C on B
  )
  ODEsyst_3sp_TC <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fC <- acb*B/(1+acb*B*thcb);      # Type II Functional response of C on prey B
      fB <- aba*A/(1+aba*A*thba);      # Type II Functional response of B on prey A
      dA <- ra*A*(1-A/Ka)-B*fB;        # Biomass changes of basal species A
      dB <- e*B*fB-(C*fC+xb*B);        # Biomass changes of intermediate consumer species B
      dC <- e*C*fC-xc*C;               # Biomass changes of invading top species C
      list(c(dA,dB,dC));
    })
  }
  EIG      <- eigen(jacobian.full(y = c(A=Bx, B=By, C=Bz), func = ODEsyst_3sp_TC, parms = as.list(parameters)))$values
  return(max(Re(EIG)))
}
### 3 sp system, Intraguild predation
lambda_3sp_IGP <- function(Temp, Car, Bx, By, Bz, Mx, My, Mz, Alpha, Beta, Gamma) {
  parameters <- c(
    ra  <- r(Temp, Mx),                # Growth rate of A
    Ka  <- K(Temp, Mx, Car),           # Carrying capacity of A
    xb  <- x2(Temp, My),               # Metabolism of B
    thba <- th(Temp, My, Mx),          # Handling time of B on A
    aba  <- a(Temp, My, Mx),           # Attack rate of B on A
    xc  <- x2(Temp, Mz),               # Metabolism of C
    thcb <- th(Temp, Mz, My),          # Handling time of C on B
    thca <- th(Temp, Mz, Mx),          # Handling time of C on A
    acb  <- a(Temp, Mz, My),           # Attack rate of C on B
    aca  <- a(Temp, Mz, Mx)            # Attack rate of C on A
  )
  ODEsyst_3sp_IGP <- function(t,state,parameters){
    with(as.list(c(state,parameters)), {
      
      fCa <- aca*A/(1+aca*A*thca);     # Type II Functional response of C on prey A
      fCb <- acb*B/(1+acb*B*thcb);     # Type II Functional response of C on prey B
      fB <- aba*A/(1+aba*A*thba);      # Type II Functional response of B on prey A
      dA <- ra*A*(1-A/Ka)-(B*fB+C*fCb);# Biomass changes of basal species A
      dB <- e*B*fB-(C*fCb+xb*B);       # Biomass changes of intermediate consumer species B
      dC <- e*C*(fCa+fCb)-xc*C;        # Biomass changes of invading top species C
      list(c(dA,dB,dC));
    })
  }
  EIG      <- eigen(jacobian.full(y = c(A=Bx, B=By, C=Bz), func = ODEsyst_3sp_IGP, parms = as.list(parameters)))$values
  return(max(Re(EIG)))
}