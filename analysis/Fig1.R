################################################################################
### Fig1 : barplot + jitter of variant read / total for each mutator (A or D) ##
################################################################################

# Packages
library(ggplot2)

# Data
setwd("~/analysis") # set to source file location
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # (if using Rstudio)
load("variants.rds")

# Select the data
variants1 = variants[((variants$Reference == "C" | variants$Reference == "G") & 
                        (variants$Mutation == "T" | variants$Mutation == "A")),]
variants2 = variants[((variants$Reference == "T" | variants$Reference == "A") & 
                        (variants$Mutation == "C" | variants$Mutation == "G")),]

# Graphs
ggA = ggplot(variants1[(variants1$Batch == "A" & variants1$Region == "targeted" & variants1$Depth > 20),]) +
  aes(x = Mutator, y = Occurences/Depth) +
  ggtitle("C:G -> T:A") +
  xlab("Mutator") +
  ylab("variant read / total") +
  geom_jitter() +
  geom_boxplot(alpha = 0.5, outlier.size=0, outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                            "ABE8.20-m-T7", "pSEVA471", "CGG", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
                            "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
                            "ABE8.20-m-CGG"))

ggB = ggplot(variants1[(variants1$Batch == "D" & variants1$Region == "targeted" & variants1$Depth > 20),]) +
  aes(x = Mutator, y = Occurences/Depth) +
  ggtitle("C:G -> T:A") +
  xlab("Mutator") +
  ylab("variant read / total") +
  geom_jitter() +
  geom_boxplot(alpha = 0.5, outlier.size=0, outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                            "ABE8.20-m-T7", "pSEVA471", "CGG", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
                            "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
                            "ABE8.20-m-CGG"))

ggA2 = ggplot(variants2[(variants2$Batch == "A" & variants2$Region == "targeted" & variants2$Depth > 20),]) +
  aes(x = Mutator, y = Occurences/Depth) +
  ggtitle("A:T -> G:C") +
  xlab("Mutator") +
  ylab("variant read / total") +
  geom_jitter() +
  geom_boxplot(alpha = 0.5, outlier.size=0, outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                            "ABE8.20-m-T7", "pSEVA471", "CGG", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
                            "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
                            "ABE8.20-m-CGG"))

ggB2 = ggplot(variants1[(variants2$Batch == "D" & variants2$Region == "targeted" & variants2$Depth > 20),]) +
  aes(x = Mutator, y = Occurences/Depth) +
  ggtitle("A:T -> G:C") +
  xlab("Mutator") +
  ylab("variant read / total") +
  geom_jitter() +
  geom_boxplot(alpha = 0.5, outlier.size=0, outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                            "ABE8.20-m-T7", "pSEVA471", "CGG", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
                            "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
                            "ABE8.20-m-CGG"))
# Plot
ggA
ggB
ggA2
ggB2

# Save
ggsave("Figure/Fig.1-b.png", 
       plot = ggB,
       width = 36, 
       height = 25, 
       units = "cm",
       dpi = 300)
ggsave("Figure/Fig.1-a.png", 
       plot = ggA,
       width = 36, 
       height = 25, 
       units = "cm",
       dpi = 300)
ggsave("Figure/Fig.1-c.png", 
       plot = ggB2,
       width = 36, 
       height = 25, 
       units = "cm",
       dpi = 300)
ggsave("Figure/Fig.1-d.png", 
       plot = ggA2,
       width = 36, 
       height = 25, 
       units = "cm",
       dpi = 300)
