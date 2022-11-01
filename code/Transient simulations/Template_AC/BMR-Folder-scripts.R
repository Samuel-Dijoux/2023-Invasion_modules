########################################################################
### R Script template to create/delete the folder directories containing the scripts & outputs of simulations
########################################################################

### Create a path directing where will be created all "BMRs" directory folders:
# example "./Template_AC";
Path <- "ZZ HERE ZZ";

### Creation of Directory for simulations
Beta <- c(1, 2, 5, 10);
Alpha <- c(0.1, 0.2, 0.5, 0.75, 1, 1.5, 2, 5, 10);

for(alpha in Alpha){
  for(beta in Beta){
    gamma <- alpha*beta;
    if(gamma<1){next
      }else{
      bmr <- paste("BMRs",paste(alpha, beta, sep='.'), sep="_");
      dir1 <- paste(Path,bmr,sep='/');
      dir2 <- paste(dir1,"out",sep='/');
      if(!file.exists(dir1)){
        dir.create(dir1);
        dir.create(dir2);
      }
      }
  }
}

### Automatised command to delete the folders
for(alpha in Alpha){
  for(beta in Beta){
    gamma <- alpha*beta;
    if(gamma<1){next
    }else{
    bmr <- paste("BMRs",paste(alpha, beta, sep='.'), sep="_");
    setwd(Path); if(file.exists(bmr)){unlink(bmr, recursive=TRUE)};
    }
  }
}