########################
## ChromHMM main call ##
########################

# Description: Function to call ChromHMM.jar
runChromHmmLoop <- function(in.dir, out.dir, chrom.dir,
                            min.nr.states = 2, 
                            max.nr.states = 35, 
                            by.states = 3,
                            max.iterations = 1000, 
                            cores) {
  for(nr.states in seq(from = min.nr.states, to = max.nr.states, by = by.states)) {
    cmd <- paste("java -jar", chrom.dir, "LearnModel -p", cores, 
                 "-holdcolumnorder -r", max.iterations, 
                 "-nobrowser -noenrich -printstatebyline ", 
                 in.dir, out.dir, nr.states, "mm10", sep = ' ')
    system(command = cmd)
  }
}

