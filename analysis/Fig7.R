################################################################################
### Fig7 : Modelisation figures #####
################################################################################

# Packages
library(ggplot2)
library(stringr)

# Data
setwd("~/Documents/iGEM/analyse_ngs/evolutionT7/analysis") # set to source file location
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # (if using Rstudio)
load("variants.rds")

# Parse
gg = variants[(variants$Depth > 20 & variants$Batch == "A" & variants$Region == "targeted" & !(variants$Mutator %in% c("CGG", "pSEVA471",
                                                                                                                       "pSEVA221", "T7"))),]
mutator = c()
mutrate = c()
reac = c()
mutrate2 = c()
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
gg$Mutation_rate = gg$Mutation_rate/40 # per generations
gg$Mutator = as.factor(gg$Mutator)

tmp = c()
for (i in levels(gg$Mutator)){
  if (i %in% c("AID-T7", "pmCDAI-T7", "rAPOBEC1-T7")){
    tmp = c(tmp, sum(gg$Mutation_rate[gg$Mutator == i & gg$Activity == "C -> T"]))
    tmp = c(tmp, sum(gg$Mutation_rate[gg$Mutator == i & gg$Activity == "G -> A"]))
  }
  if (i == "TadA*-T7"){
    
    tmp = c(tmp, sum(gg$Mutation_rate[gg$Mutator == i & gg$Activity == "T -> C"]))
    tmp = c(tmp, sum(gg$Mutation_rate[gg$Mutator == i & gg$Activity == "A -> G"]))
  }
}
esti = c(0.00354, 0.00360, 0.0564, 0.00664, 0.00553, 0.00972, 0.00452, 0.000709) #estimated values (spanich publication)
miniplot = data.frame(x = esti, y = tmp)

# Graph
## fig 1
gg1 <- ggplot(gg) +
  aes(x = Mutator, y = Mutation_rate) + 
  geom_col() +
  xlab("Mutator") +
  ylab("Mutation rate") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7"))

gg2 <- ggplot(miniplot) +
  aes(x = x, y = y) +
  xlab("Estimation at g = 3 (literature)") +
  ylab("Estimation at g = 45") +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype="dashed", color = "red")

gg3 <- ggplot(miniplot) +
  aes(x = x, y = y*32 -0.013) +
  xlab("Estimation at g = 3 (literature)") +
  ylab("Corrected estimation at g = 45") +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype="dashed", color = "red")


# Save
ggsave("Figure/Fig.7-a.png", 
       plot = gg1,
       width = 25, 
       height = 20, 
       units = "cm",
       dpi = 300)
ggsave("Figure/Fig.7-b.png", 
       plot = gg2,
       width = 20, 
       height = 20, 
       units = "cm",
       dpi = 300)
ggsave("Figure/Fig.7-c.png", 
       plot = gg3,
       width = 20, 
       height = 20, 
       units = "cm",
       dpi = 300)