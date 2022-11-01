########################################################################
### R Script template to create/delete the folder directories containing the scripts & outputs of simulations
########################################################################

### Create a path directing where will be created all "BMRs" directory folders:
# example "./Template_2spAB";
Path <- "ZZ HERE ZZ";

### Creation of Directory for simulations
Alpha <- c(1, 2, 5, 10);

for(alpha in Alpha){
  bmr <- paste("BMRs",alpha, sep="_");
  dir1 <- paste(Path,bmr,sep='/');
  dir2 <- paste(dir1,"out",sep='/');
  if(!file.exists(dir1)){
    dir.create(dir1);
    dir.create(dir2);
  }
}


### Automatised command to delete the folders
for(alpha in Alpha){
  bmr <- paste("BMRs",alpha, sep="_");
  dir1 <- paste(Path,bmr,sep='/');
  dir2 <- paste(dir1,"out",sep='/');
  if(file.exists(dir1)){
    dir.delete(dir2);
    dir.delete(dir1);
  }
}