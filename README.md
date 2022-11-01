# data

This directory houses all the datasets generated and used in our study. The datasets are organized into four folders, each described below:  
 * **1_generated_data**: raw data generated from the transient simulations,
 * **2_complete_data**: complete datasets gathering transient and equilibrium analyses;
 * **3_BMR_gradients_data**: data of invasion influences on resident communities over gradients of species body mass ratio;
 * and **4_POP_data**: data of species properties for the 10 last years of the transient simulations across trophic modules.

## 1. Data generated from the transient analyses.

The data stored in the **1_generated_data** folder are the raw data generated from the transient simulations. They consist in seven datasets of (non-invaded and invaded) communities structure after 5000 years.  
 * _2spAB.txt_ & _2spAC_: Two-species resident community composed of a consumer species (species B or C) preying on its basal resource (species A).  
 * _AC.data.txt_: Apparent competition module, in which a basal resource species (species B) invades the resident system and competes with resident basal species (species A).  
 * _EC.data.txt_: Exploitative competition module, in which a consumer species (species B) invades the resident system and competes with resident consumer (species C).  
 * _TC.data.txt_: Trophic chain module, in which an apex predator (species C) invades the resident system.  
 * _IGPB.data.txt_: Intraguild predation module with an invasion of intraguild prey species (species B).  
 * _IGPC.data.txt_: Intraguild predation module with an invasion of intraguild predator (species C).  
 
 Each of these datasets present the same structure, with 13 to 19 variables:
 
 | Variables | Units | Descriptions |
 | --------- | ----- | ------------ |
 | Time 		| s | Time duration of the simulation |
 | Temperature 	| °C | Temperature |
 | Carr | g.L<sup>-1</sup>| Intercept of the carrying capacity _I<sub>$\Phi$_ |
 | A, B, C | g.m<sup>-2</sup> | Species biomass densities |
 | Afin, Bfin, Cfin | - | Species presence/absence at the end of the transient dynamics |
 | Tot | - | Number of species present at the end of the transient dynamics |
 | ExtTimeA, B or C| s | Species extinction time during the simulation |
 | M<sub>A</sub>, M<sub>B</sub>, M<sub>C</sub>| mg | Species body masses |
 | Alpha, Beta, Gamma | - | Species mass ratio |
 
## 2. _Tot_Data_: Complete data sets (transient and equilibrium analyses).
The data stored in the **_2_complete_data_** folder extend the _generated data_ from the transient analyses with additional equilibrium analyses. They denote any changes in the stability regime within modules (before & end of the simulations) and between resident and invaded systems (under similar biotic and abiotic conditions). They also denote the influences of invading species on local communities through biodiversity change and involved invasion mechanisms.  
For a better readibility, we named each species by their trophic position (R= resource, C= consumer, P=predator) and their status (r=resident and i for invader) instead of their initial letters (A, B and C). The seven datasets are entitled as following:

* _CR_Tot_Alpha_ & _CR_Tot_Gamma_: respectively extend _2spAB.txt_ & _2spAC_, _alpha_ refering to the mass ratio between species A and B and _gamma_ the mass ratio between species A and C. 
* _AC_Tot_data.txt_, _EC_Tot_data.txt_, _TC_Tot_data.txt_, _IGPB_Tot_data.txt_ & _IGPC_Tot_data.txt_ extend the previous data sets of invaded communities (AC, EC, TC, IGPB and IGPC).

All datasets have 26 common variables. The datasets of invaded communities have 15 additional variables, highlighted as _<sup>*</sup>Variable_ in the following Table. We also used _(species)_ to ease the nomenclature. It actually accounts for 2-3 columns in the datasets depending on the number of interacting species (i.e., 2 species in resident community without invasion, and 3 species otherwise).

| Variables | Units | Descriptions |
| --------- | ----- | ------------ |
| Temp 			| °C 	| Temperature |
| Car 			| g.L<sup>-1</sup> | Intercept of the carrying capacity _I<sub>$\Phi$_ |
| <sup>*</sup>Fw| - | Food web module of study |
| <sup>*</sup>Inv.position | - | Trophic position of the invading species (1=basal, 2=consumer, 3=predator) |
| Res_BMR 		| - | Body mass ratio between resident species |
| Alpha, Beta, Gamma | - | Species body mass ratio |
| M_(species) 	| mg | Species body mass |
| (species)_i 	| g.m<sup>-2</sup> | Species biomass density before the transient dynamics ( at Cr-Rr equilibrium) |
| (species)_f 	| g.m<sup>-2</sup> | Species biomass density at the end of the transient dynamics |
| (species)_fin | -	| Species presence/absence at the end of the transient dynamics |
| Tot 			| - | Total number of species remaining at the end of the transient dynamics |
| (species)_extime | s | Species extinction time during the simulation |
| Sim_duration	| s | Time duration of the simulation |
| (species)_eq	| g.m<sup>-2</sup> | Species biomass density at equilibrium |
| Dom_Eigenvalue| - | Real part of the dominant eigenvalue of the Jacobian matrix |
| Tot_eq		| - | Number of species present in the studied module|
| Dyn_state_i | - | Structure of the community (composition and stability regime) before the transient dynamics |
| Dyn_state_f | - | Structure of the community (composition and stability regime) after transient & equilibrium analyses |
| Zone_i | - | Codification of _Dyn_state_i_ used for illustration/comparison in the analyses (from 1-3)|
| Zone_f | - | Codification of _Dyn_state_f_ used for illustration/comparison in the analyses (from 0-9),
| EQ_transit_InitFin | - | Change in the regime stability state (Eq=equilibrium, Cycles=cycles, Null=collapse) of a module before/after analyses |
| Zone_EQtransit_InitFin | - | Codification of _EQ_transit_InitFin_ used for illustration/comparison in the analyses (from 1-5)|
| <sup>*</sup>EQ_transit_ResInv | - | Change in the regime stability of resident community before/after invasion (_State<sub>r</sub>.State<sub>i</sub>_)|
| <sup>*</sup>Zone_EQtransit_ResInv | - | Codification of _EQ_transit_ResInv_ used for illustration/comparison in the analyses (from 1-9)|
| <sup>*</sup>BioD | - | Change in diversity resulting from invasion (from -2 to 3)|
| <sup>*</sup>State | - | Invasion mechanism involved|
| <sup>*</sup>Res_zonef | - | _Zone_f_ of resident data set for same abiotic conditions and species mass ratio |



## 3. Responses over BMR gradients

Tables S9 and S10 in our study result from the compilation of smaller datasets that are generated when using the script _BMR.MeanProp.R_ (see in **_code_** folder). These sets are stored in the **3_BMR_gradients_data** folder. It generated the following data:

 * The file extensions "_IS.txt" describe the averaged percentages for each invasion mechanisms and diversity change over gradient of the studied body mass ratio. The diversity changes variables _Bio_gain_, _Bio_neutral_ & _Bio_loss_ correspond respectively to the columns $\Delta$D>0,  $\Delta$D=0,  $\Delta$D<0 in Table S10.

 * The file extensions "RS-ResInv" describe the averaged percentages for each qualitative change in the stability regime of the community due to species invasion over gradients of species body mass ratio. The variables _DeltaS_gain_, _DeltaS_neutral_ & _DeltaS_loss_ correspond respectively to the columns $\Delta$S>0,  $\Delta$S=0,  $\Delta$S<0 in Table S11.

## 4. POP data

The datasets stored in the **_4. POP data_** folder consist of records of the population properties for each species in each invaded module over the last 10 years of the transient dynamics: species minimal, maximal and average biomass densities, its standard deviation and coefficient of variation for the sudied year. These data are only subsets of much larger datasets (80200 abiotic conditions* 16-25 bmr* 10 years). They represents populations properties across a gradient of one abiotic condition while fixing the other (either temperature or nutrient level), for all studied body mass ratio and food web modules:
* _N_data_: data over temperature gradient for fixed intercept of the carrying capacity at 5 g.L<sup>-1</sup>.
* _T_data_: data over gradient of nutrient for fixed temperature at 10°C.

