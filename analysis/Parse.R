############################################################################
### Parse the output of the processing pipeline to extract readable data ###
############################################################################

# Packages

# Set directory
setwd("../processing/variant40")

# Intitialisation
file=list.files(pattern = "122_10")
sample = c()
pos = c()
ref = c()
var = c()
read_mut = c()
read_tot = c()
repet = c()

# Parse
for (i in file){
  tmp = read.csv2(i, sep = "\t")
  sample = c(sample, rep(substring(i, 1,11), nrow(tmp)))
  pos = c(pos, tmp$Position)
  ref = c(ref, tmp$Ref)
  var = c(var, tmp$VarAllele)
  read_mut = c(read_mut, tmp$Reads2)
  read_tot = c(read_tot, tmp$Reads1 + tmp$Reads2)
  repet = c(repet, rep(substring(i, 1,1), nrow(tmp)))
}

data = data.frame(Sample = sample, 
                  Position = pos, 
                  Reference = ref, 
                  Mutation = var,  
                  Occurences = read_mut, 
                  Depth = read_tot,
                  Batch = repet)

# Export
write.csv2(data, file = "donneees_ngs_variantcall.csv", row.names = FALSE) # csv
variants = data
save(variants, file = "../../analysis/variants.rds")
