################################################################################
####  Identification of each mutations impact on IC50 #########################
################################################################################

# Packages
library(tidyr)

# Data
tab <- read.csv("mutations_azt.csv", sep = ",", header = TRUE, dec = ",")
row.names(tab) = tab[,1]
tab = tab[,-1]
tab <- data.frame(apply(tab, 2, function(x) replace_na(x, 0)))
tab <- data.frame(apply(tab, 2, function(x) as.factor(x)))

# mod
mod <- lm( IC50 ~ L49F * E61K * E102K * E108K * V117I * R118K * R162C * R162H * 
             D207N * G236S * E237K * G263R * S264N * R273K, data = tab)
mod2 <- lm( IC50 ~ L49F + E61K + E102K + E108K + V117I + R118K + R162C + R162H + 
             D207N + G236S + E237K + G263R + S264N + R273K, data = tab)

mod3 <- lm(IC50 ~ P12S + A16T + L19F + D36N + L49F + E61K + S80F + E102K + E108K + V117I + 
R118K + R162C + R162H + V182I + G194S + E195K + D207N + G236S + E237K + D250N + 
G251S + G263R + S264N + E270K + R273K + E277K, data = tab)

# only known mutations
tab_main = tab[(tab$P12S == 0 &
                  tab$A16T == 0 &
                  tab$L19F == 0 &
                  tab$D36N == 0 &
                  tab$L49F == 0 &
                  tab$E61K == 0 &
                  tab$S80F == 0 &
                  tab$E108K == 0 &
                  tab$V117I == 0 &
                  tab$R118K == 0 &
                  tab$V182I == 0 &
                  tab$G194S == 0 &
                  tab$E195K == 0 &
                  tab$D207N == 0 &
                  tab$D250N == 0 &
                  tab$G251S == 0 &
                  tab$G263R == 0 &
                  tab$S264N == 0 &
                  tab$E270K == 0 &
                  tab$R273K == 0 &
                  tab$E277K == 0),]

tab_main = data.frame(as.factor(tab_main$E102K), 
                      as.factor(tab_main$R162C), 
                      as.factor(tab_main$R162H), 
                      as.factor(tab_main$G236S), 
                      as.factor(tab_main$E237K), 
                      tab_main$IC50)
colnames(tab_main) = c("E102K", "R162C", "R162H", "G236S", "E237K", "IC50")
tab_main = tab_main[(apply(tab_main[1:5], 1, function(x) sum(as.integer(x))) > 0),]
tab_main$IC50 = as.numeric(tab_main$IC50)


mod = lm(IC50 ~ E102K * R162C * G236S * E237K, data = tab_main)
mod = lm(IC50 ~ E102K * R162H * R162C * G236S * E237K, data = tab_main) #best
mod = lm(IC50 ~ E102K + R162H + R162C + G236S + E237K, data = tab_main)

# Predictions
predic_mod = predict(mod, newdata = tab, interval = "prediction")
ecart = round(as.numeric(tab$IC50) - predic_mod[,1] , digits = 4)
dedu = data.frame(tab[,-c(8,12,13,18,19)])
dedu[["PredIC50"]] = ecart
dedu = dedu[(apply(dedu[,1:21], 1, function(x) sum(as.numeric(x))) > 0),]
write.csv(dedu, file = "New_mutations.csv")
