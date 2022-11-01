########################################################################
### R Script template to create/delete the folder directories containing the scripts & outputs of simulations
########################################################################

### Create a path directing where will be created all "BMRs" directory folders:
# example "./Template_2spAC";
Path <- "ZZ HERE ZZ";

### Creation of Directory for simulations
Gamma <- c(1, 2, 2.5, 3, 3.5, 3.75, 4, 5, 7.5, 10, 15, 20, 25, 50, 100);

for(gamma in Gamma){
  bmr <- paste("BMRs",gamma, sep="_");
  dir1 <- paste(Path,bmr,sep='/');
  dir2 <- paste(dir1,"out",sep='/');
  if(!file.exists(dir1)){
      dir.create(dir1);
      dir.create(dir2);
      }
}

### Automatised command to delete the folders
for(gamma in Gamma){
  bmr <- paste("BMRs",gamma, sep="_");
  setwd(Path); if(file.exists(bmr)){unlink(bmr, recursive=TRUE)};
}