########################################################################
### R Script to generate or delete 40 subscripts for transient dynamics
### for each species body mass combinations
########################################################################
library(readtext);

Fw <-'2spAC'             # Food web nomenclature, no invasion in Consumer B-resource A system
Gamma <- c(1, 1.5, 2, 2.5, 3, 3.5, 3.75, 4, 5, 7.5, 10, 15, 20, 25, 50, 100);

# Create a path directing where will be created all "BMRs" directory folders. Example: "./Template_2spAC";
Path <- "ZZ HERE ZZ";
################# Automatised command to create R Scripts for all body mass combinations
bloc1 <- readtext("./R1.txt");
bloc2 <- readtext("./R2.txt");
bloc3 <- readtext("./R3.txt");
bloc4 <- readtext("./R4.txt");

for(gamma in Gamma){
  ## Select the right directory Path
   setwd(paste(Path, paste("BMRs", gamma, sep="_"), sep='/'));
  for(i in 1:40){
    texte <- paste(bloc1$text, paste( paste(Fw,gamma,i, sep='_'), 'R', sep='.'),
               bloc2$text, i,
               bloc3$text, gamma,
               bloc4$text, sep="")
    cat(texte, file=paste( paste(Fw,gamma,i, sep='_'), 'R', sep='.'))
  }
} 
################# Automatised command to delete R Scripts for all body mass combinations
for(gamma in Gamma){
    ## Select the right directory Path
    setwd(paste(Path, paste("BMRs", gamma, sep="_"), sep='/'));
    for(i in 1:40){
      fn <- paste( paste(Fw,gamma,i, sep='_'), 'R', sep='.');
      ## Check its existence
      if (file.exists(fn)){ 
        ## Delete file if it exists
        file.remove(fn)}
  }
}
