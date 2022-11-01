########################################################################
### R Script template to generate/delete the Functions.R
### for each species body mass combinations
########################################################################
library(readtext);

fun <- readtext("./Functions.txt");
fn <-"Functions.R";
Alpha <- c(1, 2, 5, 10);
Beta <- c(0.1, 0.2, 0.5, 0.75, 1, 1.5, 2, 5, 10);

# Create a path directing where will be created all "BMRs" directory folders. Example: "./Template_EC";
Path <- "ZZ HERE ZZ";

################# Automatised command to create Functions.R scripts for all body mass combinations
for(alpha in Alpha){
  for(beta in Beta){
    gamma <- (alpha * beta);
    if(gamma < 1){ next 
    }else{
      setwd(paste(Path, paste("BMRs", paste(alpha, beta, sep='.'), sep="_"), sep='/'));
      texte <- fun;
      cat(texte[,2], file='Functions.R');
    }
  }
} 

################# Automatised command to delete Functions.R scripts for all body mass combinations
for(alpha in Alpha){
  for(beta in Beta2){
    gamma <- (alpha * beta);
    if(gamma < 1){ next 
    }else{
      setwd(paste(Path, paste("BMRs", paste(alpha, beta, sep='.'), sep="_"), sep='/'));
      ## Check its existence
      if(file.exists(fn)){
        ## Delete file if it exists
        file.remove(fn)}
    }
  }
}
