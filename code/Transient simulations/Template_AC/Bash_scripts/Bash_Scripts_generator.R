########################################################################
### R Script template to generate/delete 40 bash subscripts for each species body mass combinations
### These scripts are necessary when performing large amounts of simulations 
### using computing infrastructure (such as MetaCentrum)
########################################################################

library(readtext);

Fw <- 'AC';             # Food web nomenclature, invasion of basal species
Beta <- c(1, 2, 5, 10);
Alpha <- c(0.1, 0.2, 0.5, 0.75, 1, 1.5, 2, 5, 10);

# Create a path directing where will be created all "BMRs" directory folders. Example: "./Template_AC";
Path <- "ZZ HERE ZZ";

################# Automatised command to delete Bash scripts for all body mass combinations
bloc1 <- readtext("./BS1.txt");
bloc2 <- readtext("./BS2.txt");
bloc3 <- readtext("./BS3.txt");
bloc4 <- readtext("./BS4.txt");
bloc5 <- readtext("./BS5.txt");

for(beta in Beta){
  for(alpha in Alpha2){
    gamma <- (alpha * beta);
    if(gamma < 1){ next 
    }else{
      dir <- paste(Path, paste("BMRs",paste(alpha, beta, sep='.'), sep="_"), sep='/')
      setwd(dir);
      
      for(i in 1:40){
        texte <- paste(bloc1$text, paste( paste(Fw,paste(alpha, beta, sep='.'),i, sep='_'), 'sh', sep='.'),
                       bloc2$text, paste("BMRs", paste(alpha, beta, sep='.'), sep='_'),
                       bloc3$text, i,
                       bloc4$text, alpha,
                       bloc5$text, beta,
                       bloc6$text, sep="");
        cat(texte, file= paste( paste(Fw,paste(alpha, beta, sep='.'),i, sep='_'), 'sh', sep='.'));
      }
    }
  }
}
################# Automatised command to delete Bash scripts for all body mass combinations
for(beta in Beta){
  for(alpha in Alpha2){
    gamma <- (alpha * beta);
    if(gamma < 1){ next 
    }else{
      ## Select the right directory Path
      dir <- paste(Path, paste("BMRs",paste(alpha, beta, sep='.'), sep="_"), sep='/')
      setwd(Path);
      
      for(i in 1:40){
        fn <- paste( paste(Fw,paste(alpha, beta, sep='.'),i, sep='_'), 'sh', sep='.');
        ## Check its existence
        if (file.exists(fn)){ 
          ## Delete file if it exists
          file.remove(fn)}
      }
    }
  }
}