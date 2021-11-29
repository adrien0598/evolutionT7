################################################################################
### Fig3 : Mutation (/Kb) for each mutator and each polymerase #################
################################################################################

# Packages
library(ggplot2)

# Data
setwd("~/analysis") # set to source file location
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # (if using Rstudio)
load("variants.rds")

# Calc mut
tmp = data.frame(variants[(variants$Region == "targeted" & variants$Depth > 30),])
tmp$Mutator = as.factor(tmp$Mutator)

mutations = c()
mutator = c()
for (i in levels(tmp$Mutator)){
  mutator = c(mutator, as.character(i))
  mutations = c(mutations, (sum(tmp$Occurences[tmp$Mutator == i & tmp$Region == j])/sum(tmp$Depth[tmp$Mutator == i & tmp$Region == j]))*1000)
}
tmp = data.frame(Mutator = mutator, Mutations = mutations)

# Graphs
gg <- ggplot(tmp) +
  aes(x = Mutator, y = Mutations) +
  xlab("Mutator") +
  ylab("Mutations (/Kb)") +
  geom_col() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                            "ABE8.20-m-T7", "pSEVA471", "CGG", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
                            "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
                            "ABE8.20-m-CGG"))

#Parse2
tmp = data.frame(variants[(variants$Region == "targeted" & variants$Depth > 30),])
tmp$Mutator = as.factor(tmp$Mutator)

mutations = c()
mutator = c()
batch = c()
for (i in levels(tmp$Mutator)){
  mutator = c(mutator, rep(as.character(i),2))
  mutations = c(mutations, (sum(tmp$Occurences[tmp$Mutator == i & tmp$Region == j & tmp$Batch == "A"])/sum(tmp$Depth[tmp$Mutator == i & tmp$Region == j & tmp$Batch == "A"]))*1000)
  batch = c(batch, "A")
  mutations = c(mutations, (sum(tmp$Occurences[tmp$Mutator == i & tmp$Region == j & tmp$Batch == "D"])/sum(tmp$Depth[tmp$Mutator == i & tmp$Region == j & tmp$Batch == "D"]))*1000)
  batch = c(batch, "D")
}
tmp = data.frame(Mutator = mutator, Mutations = mutations, Batch = batch)

#Graphs2
gg2 <- ggplot(tmp) +
  aes(x = Mutator, y = Mutations, fill = Batch) +
  xlab("Mutator") +
  ylab("Mutations (/Kb)") +
  geom_col(position = "dodge") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                            "ABE8.20-m-T7", "pSEVA471", "CGG", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
                            "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
                            "ABE8.20-m-CGG"))

# Plot
gg
gg2

# Save
ggsave("Figure/Fig.3.png", 
      plot = gg,
      width = 35, 
      height = 20, 
      units = "cm",
      dpi = 300)
ggsave("Figure/Fig.3*.png", 
       plot = gg2,
       width = 35, 
       height = 20, 
       units = "cm",
       dpi = 300)