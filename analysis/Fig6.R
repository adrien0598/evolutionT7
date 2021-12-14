################################################################################
### Fig6 :Mutation rate #####
################################################################################

# Packages
library(ggplot2)
library(stringr)

# Data
setwd("~/Documents/iGEM/analyse_ngs/evolutionT7/analysis") # set to source file location
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # (if using Rstudio)
load("variants.rds")
nb_C = 709
nb_T = 769
nb_A = 793
nb_G = 719

# Parse
gg = variants[(variants$Depth > 20 & variants$Batch == "A" & variants$Region == "targeted" & !(variants$Mutator %in% c("CGG", "pSEVA471",
                                                                          "pSEVA221", "T7"))),]
mutator = c()
mutrate = c()
reac = c()
gg$Mutator = as.factor(gg$Mutator)
for (i in levels(gg$Mutator)){
  mutator = c(mutator, rep(as.character(i),4))
  reac = c(reac, c("C -> T", "T -> C", "G -> A", "A -> G"))
  print(mean(gg$Depth[(gg$Reference == "C" & 
                         gg$Mutation == "T" & 
                         gg$Mutator == i)]))
  mutrate = c(mutrate, mean(gg$Occurences[(gg$Reference == "C" & 
                                            gg$Mutation == "T" & 
                                            gg$Mutator == i)]) /
                mean(gg$Depth[(gg$Reference == "C" & 
                                        gg$Mutation == "T" & 
                                        gg$Mutator == i)]))
  mutrate = c(mutrate, mean(gg$Occurences[(gg$Reference == "T" & 
                                            gg$Mutation == "C" & 
                                            gg$Mutator == i)]) /
                mean(gg$Depth[(gg$Reference == "T" & 
                                   gg$Mutation == "C" & 
                                   gg$Mutator == i)]))
  mutrate = c(mutrate, mean(gg$Occurences[(gg$Reference == "G" & 
                                            gg$Mutation == "A" & 
                                            gg$Mutator == i)]) /
                mean(gg$Depth[(gg$Reference == "G" & 
                                   gg$Mutation == "A" & 
                                   gg$Mutator == i)]))
  mutrate = c(mutrate, mean(gg$Occurences[(gg$Reference == "A" & 
                                            gg$Mutation == "G" & 
                                            gg$Mutator == i)]) /
                mean(gg$Depth[(gg$Reference == "A" & 
                                   gg$Mutation == "G" & 
                                   gg$Mutator == i)]))
}
gg = data.frame(Mutator = mutator, Mutation_rate = mutrate, Activity = reac)
gg = na.omit(gg)
gg$Mutation_rate = gg$Mutation_rate # per generations

# Graph
## modeling plot (/mutation site)
gg1 <- ggplot(gg) +
  aes(x = Mutator, y = (Mutation_rate/1.23)-0.0065) + 
  geom_col() +
  xlab("Mutator") +
  ylab("Mutation rate (*number of sites)") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                            "ABE8.20-m-T7", "", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
                            "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
                            "ABE8.20-m-CGG"))

ggsave("Figure/Fig.6-a.png", 
       plot = gg1,
       width = 30, 
       height = 21, 
       units = "cm",
       dpi = 300)

## engineering plot  (/kb)
gg$Mutation_rate[gg$Activity == "C -> T"] = (gg$Mutation_rate[gg$Activity == "C -> T"]/nb_C)*1000
gg$Mutation_rate[gg$Activity == "T -> C"] = (gg$Mutation_rate[gg$Activity == "T -> C"]/nb_T)*1000
gg$Mutation_rate[gg$Activity == "A -> G"] = (gg$Mutation_rate[gg$Activity == "A -> G"]/nb_A)*1000
gg$Mutation_rate[gg$Activity == "G -> A"] = (gg$Mutation_rate[gg$Activity == "G -> A"]/nb_G)*1000

gg2 <- ggplot(gg) +
  aes(x = Mutator, y = Mutation_rate, fill = Activity) + 
  geom_col() +
  xlab("Mutator") +
  ylab("Mutation rate") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                            "ABE8.20-m-T7", "", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
                            "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
                            "ABE8.20-m-CGG"))

gg$Mutation_rate[((gg$Activity == "C -> T" | gg$Activity == "G -> A") & 
                    gg$Mutator %in% c("TadA*-T7", "ABE8.20-m-T7",
                                      "TadA*-CGG", "ABE8.20-m-CGG"))] = NA
gg$Mutation_rate[((gg$Activity == "T -> C" | gg$Activity == "A -> G") & 
                    !(gg$Mutator %in% c("TadA*-T7", "ABE8.20-m-T7",
                                      "TadA*-CGG", "ABE8.20-m-CGG")))] = NA
gg = na.omit(gg)
gg3 <- ggplot(gg) +
  aes(x = Mutator, y = Mutation_rate, fill = Activity) + 
  geom_col() +
  xlab("Mutator") +
  ylab("Mutation rate") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                            "ABE8.20-m-T7", "", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
                            "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
                            "ABE8.20-m-CGG"))
# Save
ggsave("Figure/Fig.6-b.png", 
       plot = gg3,
       width = 30, 
       height = 21, 
       units = "cm",
       dpi = 300)