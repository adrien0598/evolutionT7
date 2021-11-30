################################################################################
###  Tab1 : top 10 of mutator/polymerase among diverse criterias ###############
################################################################################

# Packages
library(stringr)

# Data
setwd("~/Documents/iGEM/analyse_ngs/evolutionT7/analysis") # set to source file location
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # (if using Rstudio)
load("variants.rds")

# Parse
type = c()
type_only = c()
sens = c()
pol = c()
mutator_type = c()
adenine = c("TadA*-T7", "ABE8.20-m-T7","TadA*-CGG", "ABE8.20-m-CGG")
non = c("pSEVA221", "T7", "pSEVA471", "CGG")
tmp <- data.frame(variants[(variants$Depth > 20),])
for (i in 1:nrow(tmp)) {
  if (substring(tmp$Mutator[i], nchar(tmp$Mutator[i])-2, nchar(tmp$Mutator[i])) == "CGG"){
    pol = c(pol, "CGG")
  }
  else if (substring(tmp$Mutator[i], nchar(tmp$Mutator[i])-1, nchar(tmp$Mutator[i])) == "T7"){
    pol = c(pol, "T7")
  }
  else {
    pol = c(pol, "no")
  }
  if (tmp$Reference[i] == "C" & tmp$Mutation[i] == "T"){
    type = c(type, "Cytosine_forward")
    type_only = c(type_only, "Cytosine")
    sens = c(sens, "Forward")
  }
  else if (tmp$Reference[i] == "G" & tmp$Mutation[i] == "A"){
    type = c(type, "Cytosine_reverse")
    type_only = c(type_only, "Cytosine")
    sens = c(sens, "Reverse")
  }
  else if (tmp$Reference[i] == "A" & tmp$Mutation[i] == "G"){
    type = c(type, "Adenine_forward")
    type_only = c(type_only, "Adenine")
    sens = c(sens, "Forward")
  }
  else if (tmp$Reference[i] == "T" & tmp$Mutation[i] == "C"){
    type = c(type, "Adenine_reverse")
    type_only = c(type_only, "Adenine")
    sens = c(sens, "Reverse")
  }
  else {
    type = c(type, "uncategorized")
    type_only = c(type_only, "uncategorized")
    sens = c(sens, "uncategorized")
  }
  if (tmp$Mutator[i] %in% non){
    mutator_type = c(mutator_type, "no")
  }
  else if (tmp$Mutator[i] %in% adenine){
    mutator_type = c(mutator_type, "Adenine")
  }
  else {
    mutator_type = c(mutator_type, "Cytosine")
  }
}
tmp[["Mutation_type"]] = type
tmp[["Orientation"]] = sens
tmp[["Type"]] = type_only
tmp[["Polymerase"]] = pol
tmp[['Mutator_type']] = mutator_type

tmp$Mutator_type = as.factor(tmp$Mutator_type)
tmp$Polymerase = as.factor(tmp$Polymerase)

# Target specificity
gg <- data.frame(Mutator = rep(c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", 
                                   "rAPOBEC1-T7", "TadA*-T7", 
                                   "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                                   "ABE8.20-m-T7", "pSEVA471", "CGG", "AID-CGG", 
                                   "pmCDA1-CGG", "rAPOBEC1-CGG",
                                   "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", 
                                   "evo-CDA1-BE4max-CGG", "ABE8.20-m-CGG"),2),
                 Region = c(rep("non-targeted", 18), rep("targeted", 18)))
mut = c()
for (i in 1:length(gg$Mutator)){
  mut = c(mut, sum(tmp$Occurences[(tmp$Mutator == gg$Mutator[i] & tmp$Region == gg$Region[i])]))
}
gg[["Mutation_number"]] = mut
tab = data.frame(Mutator = c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", 
                             "rAPOBEC1-T7", "TadA*-T7", 
                             "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                             "ABE8.20-m-T7", "pSEVA471", "CGG", "AID-CGG", 
                             "pmCDA1-CGG", "rAPOBEC1-CGG",
                             "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", 
                             "evo-CDA1-BE4max-CGG", "ABE8.20-m-CGG"))
a = c()
for (i in tab$Mutator){
  a = c(a, gg$Mutation_number[(gg$Mutator == i & gg$Region == "targeted")] /
          gg$Mutation_number[(gg$Mutator == i & gg$Region == "non-targeted")])
}
tab[['Target_specificity']] = a

# Strand and reaction specificity

trans = data.frame(mutator = c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", 
                               "rAPOBEC1-T7", "TadA*-T7", 
                               "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                               "ABE8.20-m-T7", "pSEVA471", "CGG", "AID-CGG", 
                               "pmCDA1-CGG", "rAPOBEC1-CGG",
                               "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", 
                               "evo-CDA1-BE4max-CGG", "ABE8.20-m-CGG"),
                   mutation = c(NA, NA, "Cytosine_forward", 
                                "Cytosine_forward", "Cytosine_forward",
                                "Adenine_forward", "Cytosine_forward",
                                "Cytosine_forward", "Adenine_forward",
                                NA, NA, "Cytosine_reverse",
                                "Cytosine_reverse", "Cytosine_reverse",
                                "Adenine_reverse", "Cytosine_reverse",
                                "Cytosine_reverse", "Adenine_reverse"),
                   strand = c(NA, rep("Forward", 8),
                              NA, rep("Reverse", 8)),
                   reaction = c(rep(c(NA, NA, "Cytosine", 
                                    "Cytosine", "Cytosine",
                                    "Adenine", "Cytosine",
                                    "Cytosine", "Adenine"), 2)))
mut = c()
for (i in tab$Mutator){
  mut = c(mut, 
          sum(tmp$Occurences[(tmp$Mutator == i & 
                                tmp$Mutation_type == trans$mutation[trans$mutator == i])])/
          sum(tmp$Occurences[(tmp$Mutator == i & 
                                tmp$Mutation_type != "uncategorized" & 
                                tmp$Mutation_type != trans$mutation[trans$mutator == i])]))
}
tab[["Mutation_specificity"]] = mut

# Reaction specificity
mut = c()
for (i in tab$Mutator){
  mut = c(mut, 
          sum(tmp$Occurences[(tmp$Mutator == i & 
                                tmp$Type == trans$reaction[trans$mutator == i])])/
            sum(tmp$Occurences[(tmp$Mutator == i & 
                                  tmp$Type != "uncategorized" & 
                                  tmp$Type != trans$reaction[trans$mutator == i])]))
}
tab[["Reaction_specificity"]] = mut

# Strand specificity
mut = c()
for (i in tab$Mutator){
  mut = c(mut, 
          sum(tmp$Occurences[(tmp$Mutator == i & 
                                tmp$Orientation == trans$strand[trans$mutator == i])])/
            sum(tmp$Occurences[(tmp$Mutator == i & 
                                  tmp$Orientation != "uncategorized" & 
                                  tmp$Orientation != trans$strand[trans$mutator == i])]))
}
tab[["Strand_specificity"]] = mut

# normalisation by fragment size
# map size : 4115
# targeted region size : 2991
# => non targeted region size : 1124
#norm_factor = 2991/1124
#tab[,2:5] = tab[,2:5]/norm_factor

# tot mutation number (fiability indicator)
a = c()
b = c()
for (i in tab$Mutator){
  a = c(a, sum(tmp$Occurences[(tmp$Mutator == i & tmp$Mutation_type != "uncategorized")]))
  b = c(b, length(tmp$Occurences[(tmp$Mutator == i & tmp$Mutation_type != "uncategorized")]))
}
tab[['Total mutation number']] = a
tab[['Total mutation count']] = b

# Print
summary(tab)

# Write
write.csv(tab, file = "Figure/Tab1.csv", quote = FALSE, row.names = FALSE, col.names = TRUE)
