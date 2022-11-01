########################################################################
### R Script to generate or delete 40 subscripts for transient dynamics
### for each species body mass combinations
########################################################################
library(readtext);

Fw <-'EC'             # Food web nomenclature, invasion of intermediate consumer species
Alpha <- c(1, 2, 5, 10);
Beta <- c(0.1, 0.2, 0.5, 0.75, 1, 1.5, 2, 5, 10);

# Create a path directing where will be created all "BMRs" directory folders. Example: "./Template_EC";
Path <- "ZZ HERE ZZ";
################# Automatised command to create R Scripts for all body mass combinations
bloc1 <- readtext("./R1.txt");
bloc2 <- readtext("./R2.txt");
bloc3 <- readtext("./R3.txt");
bloc4 <- readtext("./R4.txt");
bloc5 <- readtext("./R5.txt");

for(alpha in Alpha){
  for(beta in Beta2){
    gamma <- (alpha * beta);
    if(gamma < 1){ next 
    }else{
      setwd(paste(Path, paste("BMRs", paste(alpha, beta, sep='.'), sep="_"), sep='/'));
      for(i in 1:40){
        texte <- paste(bloc1$text, paste( paste(Fw,paste(alpha, beta, sep='.'),i, sep='_'), 'R', sep='.'),
                       bloc2$text, i,
                       bloc3$text, alpha,
                       bloc4$text, beta,
                       bloc5$text, sep="")
        cat(texte, file=paste( paste(Fw,paste(alpha, beta, sep='.'),i, sep='_'), 'R', sep='.'));
      }
    }
  }
}
################# Automatised command to delete R Scripts for all body mass combinations
for(alpha in Alpha){
  for(beta in Beta2){
    gamma <- (alpha * beta);
    if(gamma < 1){ next 
    }else{
      setwd(paste(Path, paste("BMRs", paste(alpha, beta, sep='.'), sep="_"), sep='/'));
      for(i in 1:40){
        fn <- paste( paste(Fw,paste(alpha, beta, sep='.'),i, sep='_'), 'R', sep='.');
        ## Check its existence
        if (file.exists(fn)){ 
          ## Delete file if it exists
          file.remove(fn)}
      }
    }
  }
}
