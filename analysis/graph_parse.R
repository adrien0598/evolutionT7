library(ggplot2)

setwd("~/Documents/iGEM/analyse_ngs/variant40")
name = c()
pos = c()
ref = c()
var = c()
read_mut = c()
read_tot = c()
file=list.files(pattern = "122_10")

for (i in file){
  tmp = read.csv2(i, sep = "\t")
  name = c(name, rep(substring(i, 1,11), nrow(tmp)))
  pos = c(pos, tmp$Position)
  ref = c(ref, tmp$Ref)
  var = c(var, tmp$VarAllele)
  read_mut = c(read_mut, tmp$Reads2)
  read_tot = c(read_tot, tmp$Reads1 + tmp$Reads2)
}

data = data.frame(Sample = name, Position = pos, Ref = ref, Mutation = var,  Occurences = read_mut, Depth = read_tot)
data = data[data$Ref %in% c("A", "T", "C", "G"),]
write.csv2(data, file = "donneees_ngs_variantcall.csv", sep = ";", dec = ",")

data$Sample = as.factor(data$Sample)
gg = data.frame(data[((data$Ref %in% c('C', 'G')) & data$Mutation %in% c("T", "A") & data$Depth > 500),])
#gg = data.frame(data[(data$Depth > 500),])

mutators = c("T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
             "TadA*-T7", "CGG", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
             "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
             "ABE8.20-m-CGG", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
             "ABE8.20-m-T7", "pSEVA221", "pSEVA471")
id = c("100", "101", "102", "103", "104", "109", "110", "111", "112",
       "113", "114", "115", "116", "117", "118", "119", "221", "471")
trans = data.frame(id = id, mut = mutators)

nom = c()
rep = c()
for (i in 1:nrow(gg)) {
  nom = c(nom, trans$mut[trans$id == substring(as.character(gg$Sample[i]), 9, 11)])
  rep = c(rep, substring(as.character(gg$Sample[i]), 1, 1))
}
gg[["nom"]] = nom
gg[['rep']] = rep

ggA = ggplot(gg[gg$rep == "A",]) +
  aes(x = nom, y = Occurences/Depth) +
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

ggB = ggplot(gg[gg$rep == "D",]) +
  aes(x = nom, y = Occurences/Depth) +
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

gg2 = data
nom = c()
type = c()
for (i in 1:nrow(gg2)) {
  nom = c(nom, trans$mut[trans$id == substring(as.character(gg2$Sample[i]), 9, 11)])
  if (((data$Ref[i] %in% c('C', 'G')) & data$Mutation[i] %in% c("T", "A"))){
    type = c(type, "C:G -> T:A")
  }
  else {
    type = c(type, "A:T -> G:C")
  }
}
gg2[["nom"]] = nom
gg2[["type"]] = type

ggplot(gg2) +
  aes(x = nom, fill = type) +
  ggtitle("total mutations") +
  xlab("Mutator") +
  ylab("mutations") +
  geom_bar() +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

id = levels(as.factor(gg2$nom))
mut = c()
idd = c()
type = c()
for (i in id){
  mut = c(mut, sum(gg2$Occurences[gg2$type == "C:G -> T:A" & gg2$nom == i])/sum(gg2$Depth[gg2$type == "C:G -> T:A" & gg2$nom == i]))
  mut = c(mut, sum(gg2$Occurences[gg2$type == "A:T -> G:C" & gg2$nom == i])/sum(gg2$Depth[gg2$type == "A:T -> G:C" & gg2$nom == i]))
  type = c(type, "C:G -> T:A", "A:T -> G:C")
  idd = c(idd, i, i)
}
gg3 = data.frame(id =idd, mut = mut, type = type)

ggplot(gg3) +
  aes(x = id, y = mut, fill = type) + 
  ggtitle("Mutations of each type") +
  xlab("Mutator") +
  ylab("Proportion") +
  geom_col() +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

ggplot(gg2[gg2$type == "A:T -> G:C",]) +
  aes(x = nom, y = Occurences/Depth) +
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

# mutations rate (validation)
gg4 = gg2[gg2$nom %in% c("AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", "TadA*-T7"),]
id = c("AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", "TadA*-T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", "TadA*-T7")
activ = c("C -> T", "C -> T", "C -> T", "T -> C", "G -> A", "G -> A", "G -> A", "A -> G")
alpha = c()
for (i in c("AID-T7", "pmCDAI-T7", "rAPOBEC1-T7")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "C" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "C" & gg4$nom == i]))
}
for (i in c("TadA*-T7")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "T" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "T" & gg4$nom == i]))
}
for (i in c("AID-T7", "pmCDAI-T7", "rAPOBEC1-T7")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "G" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "G" & gg4$nom == i]))
}
for (i in c("TadA*-T7")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "A" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "A" & gg4$nom == i]))
}

#plot
gg4 = data.frame(Mutator = id, Activity = activ, Mutation_rate = alpha)
ggplot(gg4) +
  aes(x = Mutator, y = Mutation_rate, fill = Activity) + 
  geom_col(position = "dodge") +
  xlab("Mutator") +
  ylab("Mutation rate") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(name ="Mutator", 
                   limits=c("AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
                            "TadA*-T7"))

#graphe robustesse
esti = c(0.00354, 0.0564, 0.00553, 0.00452, 0.00360, 0.00664, 0.00972, 0.000709) #estimated values (spanich publication)
gg4[["Estimation"]] = esti
ggplot(gg4) +
  aes(x = Estimation, y = Mutation_rate/1.23 - 0.0065) +
  xlab("Estimation at g = 3 (literature)") +
  ylab("Corrected estimation at g = 30") +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype="dashed", color = "red")
# wo 
esti = c(0.00354, 0.00553, 0.00452, 0.00360, 0.00664, 0.00972, 0.000709) #estimated values (spanich publication)
gg4 = gg4[-2,]
gg4[["Estimation"]] = esti
ggplot(gg4) +
  aes(x = Estimation, y = Mutation_rate) +
  ggtitle("Model accuracy") +
  xlab("Estimation at 3g (literature)") +
  ylab("Estimation at 30g (our NGS results)") +
  geom_point() +
  geom_smooth(method = "lm")


# other estimationslimits=c("pSEVA221", "T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 

gg4 = gg2[gg2$nom %in% c("evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
                         "ABE8.20-m-T7"),]
id = c("evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7","ABE8.20-m-T7", 
       "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7","ABE8.20-m-T7")
activ = c("C -> T", "C -> T", "T -> C", 
          "G -> A", "G -> A", "A -> G")
alpha = c()
for (i in c("evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "C" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "C" & gg4$nom == i]))
}
for (i in c("ABE8.20-m-T7")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "T" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "T" & gg4$nom == i]))
}
for (i in c("evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "G" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "G" & gg4$nom == i]))
}
for (i in c("ABE8.20-m-T7")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "A" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "A" & gg4$nom == i]))
}
alpha = alpha/1.23 -0.0065 # esti
alpha

# final
gg4 = gg2[gg2$nom %in% c("AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG", "TadA*-CGG", 
                         "evoAPOBEC1-BE4max-CGG","evo-CDA1-BE4max-CGG","ABE8.20-m-CGG"),]
id = c("AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
       "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
       "ABE8.20-m-CGG")
activ = c("C -> T", "C -> T","C -> T", "T -> C", "C -> T","C -> T", "T -> C", 
          "G -> A", "G -> A", "G -> A","A -> G", "G -> A", "G -> A","A -> G")
alpha = c()
for (i in c("AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "C" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "C" & gg4$nom == i]))
}
for (i in c("TadA*-CGG")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "T" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "T" & gg4$nom == i]))
}
for (i in c("evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "C" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "C" & gg4$nom == i]))
}
for (i in c("ABE8.20-m-CGG")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "T" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "T" & gg4$nom == i]))
}
for (i in c("AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "G" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "G" & gg4$nom == i]))
}
for (i in c("TadA*-CGG")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "A" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "A" & gg4$nom == i]))
}
for (i in c("evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "G" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "G" & gg4$nom == i]))
}
for (i in c("ABE8.20-m-CGG")){
  alpha = c(alpha, sum(gg4$Occurences[gg4$Ref == "A" & gg4$nom == i])/sum(gg4$Depth[gg4$Ref == "A" & gg4$nom == i]))
}

# the graph
gg0 = data.frame(Mutator = c("AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", "TadA*-T7", 
                             "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", "TadA*-T7",
                             "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7","ABE8.20-m-T7", 
                             "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7","ABE8.20-m-T7",
                             "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG", "TadA*-CGG", 
                             "evoAPOBEC1-BE4max-CGG","evo-CDA1-BE4max-CGG","ABE8.20-m-CGG",
                             "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG", "TadA*-CGG", 
                             "evoAPOBEC1-BE4max-CGG","evo-CDA1-BE4max-CGG","ABE8.20-m-CGG"),
                 Activity = c("C -> T", "C -> T", "C -> T", "T -> C", 
                              "G -> A", "G -> A", "G -> A", "A -> G",
                              "C -> T", "C -> T", "T -> C", 
                              "G -> A", "G -> A", "A -> G",
                              "C -> T", "C -> T","C -> T", "T -> C", "C -> T","C -> T", "T -> C", 
                              "G -> A", "G -> A", "G -> A","A -> G", "G -> A", "G -> A","A -> G"),
                 Mutation_rate = alpha)

ggplot(gg0) +
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


# orthogonalité
ortho1 = data.frame(data[data$Depth > 20,])

nom = c()
rep = c()
zone = c()
for (i in 1:nrow(ortho1)) {
  nom = c(nom, trans$mut[trans$id == substring(as.character(ortho1$Sample[i]), 9, 11)])
  rep = c(rep, substring(as.character(ortho1$Sample[i]), 1, 1))
  if (ortho1$Position[i] < 3126 & ortho1$Position[i] > 135){
    zone = c(zone, "targeted")
  }
  else {
    zone = c(zone, "non-targeted")
  }
}
ortho1[["nom"]] = nom
ortho1[['rep']] = rep
ortho1[['zone']] = zone
ortho1$nom = as.factor(ortho1$nom)
ortho1$zone = as.factor(ortho1$zone)
ortho1$rep = as.factor(ortho1$rep)

mutator = c()
zone = c()
mutations = c()
pol = c()
for (i in levels(ortho1$nom)){
  for (j in levels(ortho1$zone)){
    mutations = c(mutations, (sum(ortho1$Occurences[ortho1$nom == i & ortho1$zone == j])/sum(ortho1$Depth[ortho1$nom == i & ortho1$zone == j]))*1000)
    zone = c(zone, j)
    if (substring(i, nchar(i)-2, nchar(i)) == "-T7"){
      if (nchar(substring(i, 1,nchar(i)-3)) > 0){
        mutator = c(mutator, substring(i, 1,nchar(i)-3))
      }
      else {
        mutator = c(mutator, "no")
      }
      pol = c(pol, "T7")
    }
    else{
      if (nchar(substring(i, 1,nchar(i)-4)) > 0){
        mutator = c(mutator, substring(i, 1,nchar(i)-4))
      }
      else {
        mutator = c(mutator, "no")
      }
      pol = c(pol, "CGG")
    }
  }
}

ortho = data.frame(mutator = mutator, zone = zone, mutations = mutations, polymerase = pol)

ggplot(ortho[ortho$polymerase != "no",]) +
  aes(x = polymerase, y = mutations, fill = zone) +
  xlab("Polymerase") +
  ylab("Mutations (/Kb)") +
  geom_boxplot()
