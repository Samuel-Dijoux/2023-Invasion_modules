########################################################################
### R Script to generate or delete 40 subscripts for transient dynamics
### for each species body mass combinations
########################################################################
library(readtext);

Fw <-'2spAB'             # Food web nomenclature, no invasion in Consumer B-resource A system
# Create a path directing where will be created all "BMRs" directory folders. Example: "./Template_2spAB";
Path <- "ZZ HERE ZZ";
################# Automatised command to create R Scripts for all body mass combinations
bloc1 <- readtext("./R1.txt");
bloc2 <- readtext("./R2.txt");
bloc3 <- readtext("./R3.txt");
bloc4 <- readtext("./R4.txt");

for(alpha in c(1, 2, 5, 10)){
  ## Select the right directory Path
   setwd(paste(Path, paste("BMRs", alpha, sep="_"), sep='/'));
  for(i in 1:40){
    texte <- paste(bloc1$text, paste( paste(Fw,alpha,i, sep='_'), 'R', sep='.'),
               bloc2$text, i,
               bloc3$text, alpha,
               bloc4$text, sep="")
    cat(texte, file=paste( paste(Fw,alpha,i, sep='_'), 'R', sep='.'))
  }
} 
################# Automatised command to delete R Scripts for all body mass combinations
for(alpha in c(1, 2, 5, 10)){
    ## Select the right directory Path
    setwd(paste(Path, paste("BMRs", alpha, sep="_"), sep='/'));
    for(i in 1:40){
      fn <- paste( paste(Fw,alpha,i, sep='_'), 'R', sep='.');
      ## Check its existence
      if (file.exists(fn)){ 
        ## Delete file if it exists
        file.remove(fn)}
  }
}
