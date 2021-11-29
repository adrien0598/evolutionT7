################################################################################
### Fig2 : Mutation in the targeted zone vs non targeted zone (for each pol) ###
################################################################################

# Packages
library(ggplot2)

# Data
setwd("~/analysis") # set to source file location
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # (if using Rstudio)
load("variants.rds")

# Some parsing
ortho1 = data.frame(variants[variants$Depth > 20,])
ortho1$Region = as.factor(ortho1$Region)
ortho1$Mutator = as.factor(ortho1$Mutator)

mutator = c()
zone = c()
mutations = c()
pol = c()
for (i in levels(ortho1$Mutator)){
  for (j in levels(ortho1$Region)){
    mutations = c(mutations, (sum(ortho1$Occurences[ortho1$Mutator == i & ortho1$Region == j])/sum(ortho1$Depth[ortho1$Mutator == i & ortho1$Region == j]))*1000)
    zone = c(zone, j)
    print(substring(i, nchar(i)-1, nchar(i)))
    if (substring(i, nchar(i)-1, nchar(i)) == "T7"){
      print('oui')
      if (nchar(substring(i, 1,nchar(i)-3)) > 0){
        mutator = c(mutator, substring(i, 1,nchar(i)-3))
      }
      else {
        mutator = c(mutator, "no")
      }
      pol = c(pol, "T7")
    }
    else if (substring(i, nchar(i)-2, nchar(i)) == "CGG"){
      print('non')
      if (nchar(substring(i, 1,nchar(i)-4)) > 0){
        mutator = c(mutator, substring(i, 1,nchar(i)-4))
      }
      else {
        mutator = c(mutator, "no")
      }
      pol = c(pol, "CGG")
    }
    else {
      pol = c(pol, "no")
      mutator = c(mutator, "no")
    }
  }
}

ortho = data.frame(mutator = mutator, zone = zone, mutations = mutations, polymerase = pol)

# Graphs

gg <- ggplot(ortho[ortho$polymerase != "no",]) +
  aes(x = polymerase, y = mutations, fill = zone) +
  xlab("Polymerase") +
  ylab("Mutations (/Kb)") +
  geom_boxplot()

gg2 <- ggplot(ortho[ortho$polymerase != "no" & ortho$mutator != "no",]) +
  aes(x = polymerase, y = mutations, fill = zone) +
  xlab("Polymerase") +
  ylab("Mutations (/Kb)") +
  geom_boxplot()

# Plot
gg
gg2

# Save
ggsave("Figure/Fig.2.png", 
       plot = gg,
       width = 25, 
       height = 20, 
       units = "cm",
       dpi = 300)
ggsave("Figure/Fig.2*.png", 
       plot = gg2,
       width = 25, 
       height = 20, 
       units = "cm",
       dpi = 300)
