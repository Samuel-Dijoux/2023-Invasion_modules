########################################################################
### R Script template to generate/delete the Functions.R
### for each species body mass combinations
########################################################################
library(readtext);

fun <- readtext("./Functions.txt");
fn <-"Functions.R";
Gamma <- c(1, 2, 1.5, 2.5, 3, 3.5, 3.75, 4, 5, 7.5, 10, 15, 20, 25, 50, 100);

# Create a path directing where will be created all "BMRs" directory folders. Example: "./Template_2spAC";
Path <- "ZZ HERE ZZ";

################# Automatised command to create Functions.R scripts for all body mass combinations
for(gamma in Gamma){
  ## Select the right directory Path
  dir1 <- paste(Path, paste("BMRs",gamma, sep="_"), sep='/')
  if(!exists(dir1)){dir.create(dir1) };
  setwd(paste(Path, paste("BMRs",gamma, sep="_"), sep='/'));
  texte <- fun;
  cat(texte[,2], file='Functions.R');
} 

################# Automatised command to delete Functions.R scripts for all body mass combinations
for(gamma in Gamma){
  ## Select the right directory Path
  setwd(paste(Path, paste("BMRs",gamma, sep="_"), sep='/'));
  ## Check its existence
  if(file.exists(fn)){
  ## Delete file if it exists
    file.remove(fn)}
}
