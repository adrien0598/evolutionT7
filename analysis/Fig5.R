################################################################################
### Fig4 : orthogonality of direction and mutation type (mutations number) ######
################################################################################

# Packages
library(ggplot2)
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
tmp <- data.frame(variants[(variants$Depth > 20 & variants$Region == "targeted"),])
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

tmp$Mutation_type = as.factor(tmp$Mutation_type)
tmp$Orientation = as.factor(tmp$Orientation)
tmp$Mutator_type = as.factor(tmp$Mutator_type)
tmp$Polymerase = as.factor(tmp$Polymerase)

gg = data.frame(Mutator_type = c(rep("Adenine deaminase reverse (CGG)", 4), 
                                 rep("Cytosine deaminase reverse (CGG)", 4),
                                 rep("CGG without mutator", 4),
                                 rep("Control", 4),
                                 rep("Adenine deaminase forward (T7)", 4), 
                                 rep("Cytosine deaminase forward (T7)", 4),
                                 rep("T7 without mutator", 4)),
                Mutation_type = c(rep(c("Forward adenine deamination", 
                                        "Forward cytosine deamination", 
                                        "Reverse adenine deamination",
                                        "Reverse cytosine deamination"), 7)))
depth = c()
nb_mut = c()
nb_mut_dif = c()
for (i in levels(tmp$Polymerase)){
  for (j in levels(tmp$Mutator_type)){
    if (i == "no"){
      if (j == "no"){
        selection = tmp[(tmp$Polymerase == i & tmp$Mutator_type == j & tmp$Mutation_type == "Adenine_forward"),]
        depth = c(depth, sum(selection$Depth))
        nb_mut = c(nb_mut, sum(selection$Occurences))
        nb_mut_dif = c(nb_mut_dif, dim(selection)[1])
        
        selection = tmp[(tmp$Polymerase == i & tmp$Mutator_type == j & tmp$Mutation_type == "Cytosine_forward"),]
        depth = c(depth, sum(selection$Depth))
        nb_mut = c(nb_mut, sum(selection$Occurences))
        nb_mut_dif = c(nb_mut_dif, dim(selection)[1])
        
        selection = tmp[(tmp$Polymerase == i & tmp$Mutator_type == j & tmp$Mutation_type == "Adenine_reverse"),]
        depth = c(depth, sum(selection$Depth))
        nb_mut = c(nb_mut, sum(selection$Occurences))
        nb_mut_dif = c(nb_mut_dif, dim(selection)[1])
        
        selection = tmp[(tmp$Polymerase == i & tmp$Mutator_type == j & tmp$Mutation_type == "Cytosine_reverse"),]
        depth = c(depth, sum(selection$Depth))
        nb_mut = c(nb_mut, sum(selection$Occurences))
        nb_mut_dif = c(nb_mut_dif, dim(selection)[1])
      }
    }
    else {
      selection = tmp[(tmp$Polymerase == i & tmp$Mutator_type == j & tmp$Mutation_type == "Adenine_forward"),]
      depth = c(depth, sum(selection$Depth))
      nb_mut = c(nb_mut, sum(selection$Occurences))
      nb_mut_dif = c(nb_mut_dif, dim(selection)[1])
      
      selection = tmp[(tmp$Polymerase == i & tmp$Mutator_type == j & tmp$Mutation_type == "Cytosine_forward"),]
      depth = c(depth, sum(selection$Depth))
      nb_mut = c(nb_mut, sum(selection$Occurences))
      nb_mut_dif = c(nb_mut_dif, dim(selection)[1])
      
      selection = tmp[(tmp$Polymerase == i & tmp$Mutator_type == j & tmp$Mutation_type == "Adenine_reverse"),]
      depth = c(depth, sum(selection$Depth))
      nb_mut = c(nb_mut, sum(selection$Occurences))
      nb_mut_dif = c(nb_mut_dif, dim(selection)[1])
      
      selection = tmp[(tmp$Polymerase == i & tmp$Mutator_type == j & tmp$Mutation_type == "Cytosine_reverse"),]
      depth = c(depth, sum(selection$Depth))
      nb_mut = c(nb_mut, sum(selection$Occurences))
      nb_mut_dif = c(nb_mut_dif, dim(selection)[1])
    }
  }
}

gg[['Depth']] = depth
gg[['Mutation_number']] = nb_mut
gg[['Mutation_count']] = nb_mut_dif

# Graphs
gg1 <- ggplot(gg) +
  aes(x = Mutator_type, y = Mutation_number, fill = Mutation_type)+
  geom_col(position = "fill") +
  scale_fill_brewer(palette="Paired") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator type", 
                   limits=c("Adenine deaminase forward (T7)", 
                            "Cytosine deaminase forward (T7)", 
                            "T7 without mutator", 
                            "Control",
                            "Adenine deaminase reverse (CGG)", 
                            "Cytosine deaminase reverse (CGG)", 
                            "CGG without mutator")) +
  scale_y_continuous(name = "Mutation number") +
  labs(fill = "Mutation type")

# Plot
gg1
# Save

