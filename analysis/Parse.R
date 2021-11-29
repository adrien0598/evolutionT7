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

# Mutator definition
mutators = c("T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
             "TadA*-T7", "CGG", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
             "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
             "ABE8.20-m-CGG", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
             "ABE8.20-m-T7", "pSEVA221", "pSEVA471")
id = c("100", "101", "102", "103", "104", "109", "110", "111", "112",
       "113", "114", "115", "116", "117", "118", "119", "221", "471")
trans = data.frame(id = id, mut = mutators)

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
nom = c()
site = c()
for (i in 1:nrow(data)) {
  nom = c(nom, trans$mut[trans$id == substring(as.character(data$Sample[i]), 9, 11)])
  if (data$Position[i] < 3126 & data$Position[i] > 135){
    site = c(site, "targeted")
  } else {
    site = c(site, "non-targeted")
  }
}
data[["Mutator"]] = nom
data[["Region"]] = site

# Export
write.csv2(data, file = "donneees_ngs_variantcall.csv", row.names = FALSE) # csv
variants = data
save(variants, file = "../../analysis/variants.rds")
