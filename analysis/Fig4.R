################################################################################
### Fig4 : orthogonality of direction and mutation type ########################
################################################################################

# Packages
library(ggplot2)
library(stringr)

# Data
setwd("~/analysis") # set to source file location
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # (if using Rstudio)
load("variants.rds")

# Parse
type = c()
type_only = c()
sens = c()
pol = c()
tmp <- data.frame(variants[variants$Depth > 20,])
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
}
tmp[["Mutation_type"]] = type
tmp[["Orientation"]] = sens
tmp[["Type"]] = type_only
tmp[["Polymerase"]] = pol

# Graphs
gg <- ggplot(tmp[tmp$Batch == "A",]) + 
  aes(x = Mutator, fill = Mutation_type) +
  geom_bar() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                            "ABE8.20-m-T7", "pSEVA471", "CGG", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
                            "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
                            "ABE8.20-m-CGG"))

gg2 <- ggplot(tmp[tmp$Batch == "A",]) + 
  aes(x = Mutator, fill = Orientation) +
  geom_bar() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                            "ABE8.20-m-T7", "pSEVA471", "CGG", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
                            "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
                            "ABE8.20-m-CGG"))
gg3 <- ggplot(tmp[tmp$Batch == "A",]) + 
  aes(x = Mutator, fill = Type) +
  geom_bar() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                            "ABE8.20-m-T7", "pSEVA471", "CGG", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
                            "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
                            "ABE8.20-m-CGG"))

gg4 <- ggplot(tmp[tmp$Batch == "A",]) + 
  aes(x = Polymerase, fill = Orientation) +
  geom_bar()

gg5 <- ggplot(tmp[(tmp$Batch == "A" & tmp$Mutator != "T7" & tmp$Mutator != "CGG"),]) + 
  aes(x = Polymerase, fill = Orientation) +
  geom_bar()
# Plot
gg
gg2
gg3
gg4
gg5

# Save
ggsave("Figure/Fig.4-a.png", 
       plot = gg,
       width = 36, 
       height = 25, 
       units = "cm",
       dpi = 300)
ggsave("Figure/Fig.4-b.png", 
       plot = gg2,
       width = 36, 
       height = 25, 
       units = "cm",
       dpi = 300)
ggsave("Figure/Fig.4-c.png", 
       plot = gg3,
       width = 36, 
       height = 25, 
       units = "cm",
       dpi = 300)
ggsave("Figure/Fig.4-d.png", 
       plot = gg4,
       width = 25, 
       height = 20,  
       units = "cm",
       dpi = 300)
ggsave("Figure/Fig.4-d*.png", 
       plot = gg5,
       width = 25, 
       height = 20, 
       units = "cm",
       dpi = 300)



