########################################################################
### R Script template to generate/delete 40 bash subscripts for each species body mass combinations
### These scripts are necessary when performing large amounts of simulations 
### using computing infrastructure (such as MetaCentrum)
########################################################################

library(readtext);

Fw <- '2spAC';             # Food web nomenclature, consumer-resource system without invasion
# Create a path directing where will be created all "BMRs" directory folders. Example: "./Template_2spAB";
Path <- "ZZ HERE ZZ";

################# Automatised command to delete Bash scripts for all body mass combinations
bloc1 <- readtext("./BS1.txt");
bloc2 <- readtext("./BS2.txt");
bloc3 <- readtext("./BS3.txt");
bloc4 <- readtext("./BS4.txt");
bloc5 <- readtext("./BS5.txt");

for(gamma in Gamma){
  ## Select the right directory Path
  setwd(paste(Path, paste("BMRs", gamma, sep="_"), sep='/'));
  for(i in 1:40){
    texte <- paste(bloc1$text, paste( paste(Fw,gamma,i, sep='_'), 'sh', sep='.'),
                   bloc2$text, paste("BMRs", gamma, sep='_'),
                   bloc3$text, i,
                   bloc4$text, gamma,
                   bloc5$text, sep="")
    cat(texte, file= paste( paste(Fw,gamma,i, sep='_'), 'sh', sep='.'));
  }
} 
################# Automatised command to delete Bash scripts for all body mass combinations
for(gamma in Gamma){
  ## Select the right directory Path
  setwd(paste(Path, paste("BMRs", gamma, sep="_"),sep='/'));
  for(i in 1:40){
    fn <- paste( paste(Fw,gamma,i, sep='_'), 'sh', sep='.');
    ## Check its existence
    if (file.exists(fn)){ 
      ## Delete file if it exists
      file.remove(fn)}
  }
} 
