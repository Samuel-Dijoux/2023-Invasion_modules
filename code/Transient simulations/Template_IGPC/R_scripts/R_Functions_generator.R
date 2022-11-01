########################################################################
### R Script template to generate/delete the Functions.R
### for each species body mass combinations
########################################################################
library(readtext);

fun <- readtext("./Functions.txt");
fn <-"Functions.R";

# Create a path directing where will be created all "BMRs" directory folders. Example: "./Template_IGPC";
Path <- "ZZ HERE ZZ";

################# Automatised command to create Functions.R scripts for all body mass combinations
for(alpha in c(1, 2, 5, 10)){
  for(beta in c(1, 2, 5, 10)){
    setwd(paste(Path, paste("BMRs",paste(alpha, beta, sep='.'), sep="_"), sep='/'));
    texte <- fun;
    cat(texte[,2], file='Functions.R');
  }
}

################# Automatised command to delete Functions.R scripts for all body mass combinations
for(alpha in c(1, 2, 5, 10)){
  for(beta in c(1, 2, 5, 10)){
    setwd(paste(Path, paste("BMRs",paste(alpha, beta, sep='.'), sep="_"), sep='/'));
    ## Check its existence
    if(file.exists(fn)){
      ## Delete file if it exists
      file.remove(fn)}
  }
}
