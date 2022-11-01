########################################################################
### R Script template to generate/delete the Functions.R
### for each species body mass combinations
########################################################################
library(readtext);

fun <- readtext("./Functions.txt");
fn <-"Functions.R";
Beta <- c(1, 2, 5, 10);
Alpha <- c(0.1, 0.2, 0.5, 0.75, 1, 1.5, 2, 5, 10);

# Create a path directing where will be created all "BMRs" directory folders. Example: "./Template_AC";
Path <- "ZZ HERE ZZ";
################# Automatised command to create Functions.R scripts for all body mass combinations
for(beta in Beta){
  for(alpha in Alpha){
    gamma <- (alpha * beta);
    if(gamma < 1){ next 
    }else{
    ## Select the right directory Path
    dir1 <- paste(Path, paste("BMRs",gamma, sep="_"), sep='/')
    if(!exists(dir1)){dir.create(dir1) };
    setwd(paste(Path, paste("BMRs",gamma, sep="_"), sep='/'));
    texte <- fun;
    cat(texte[,2], file='Functions.R');
    }
  }
} 

################# Automatised command to delete Functions.R scripts for all body mass combinations
for(beta in Beta){
  for(alpha in Alpha){
    gamma <- (alpha * beta);
    if(gamma < 1){ next 
    }else{
    ## Select the right directory Path
    setwd(paste(Path, paste("BMRs",gamma, sep="_"), sep='/'));
    ## Check its existence
    if(file.exists(fn)){
    ## Delete file if it exists
    file.remove(fn)}
    }
  }
}
