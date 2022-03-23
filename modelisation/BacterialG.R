setwd("D:/adrie/Documents/igem")
data = read.csv2("DO_NEB_R.csv")

library(ggplot2)

gg = rbind(data.frame(time = data$time, Particles = data$E01),
           data.frame(time = data$time, Particles = data$F01),
           data.frame(time = data$time, Particles = data$G01),
           data.frame(time = data$time, Particles = data$H01))
gg[["Rep"]] = c(rep("E01", 127), rep("F01", 127), rep("G01", 127), rep("H01", 127))

ggplot(gg) + 
  aes(x = time, y = Particles, color = Rep)+
  geom_point()

data[["Mean_particles"]] = apply(data[,-5], 1, function(x) mean(x))

ggplot(data)+
  aes(x = time, y = Mean_particles)+
  geom_point() +
  geom_hline(yintercept = 4.15) +
  geom_abline(slope = 0.0348, intercept = -0.3)
# A = 4.15
# lambda = 10 min
# mu = 0.0348 min^-1

## growth model
A = 4.1
l = 60
mu = 0.00031
logistic <- function(x){
  tmp = A/(1+ exp((4*mu*(l-x))/A +2))
  return(tmp)
}

data[["Logistic"]] = sapply(data$time *60, function(x) logistic(x))

gompertz <- function(x){
  tmp = A*exp(-exp((exp(1)*mu*(l-x)/A)+1))
  return(tmp)
}

data[["Gompertz"]] = sapply(data$time *60, function(x) gompertz(x))


gg = rbind(data.frame(time = data$time, Particles = data$Mean_particles, Model = rep("Truth", 127)),
           data.frame(time = data$time, Particles = data$Logistic, Model = rep("Logistic", 127)),
           data.frame(time = data$time, Particles = data$Gompertz, Model = rep("Gompertz", 127)))

ggplot(gg[gg$Model != "Logistic",])+
  aes(x = time, y = Particles, color = Model)+
  theme_classic()+
  geom_point()

corel = data.frame(Model = c("Logistc", "Gompertz"), 
                Pearson.cor = c(cor.test(gg$Particles[gg$Model == "Logistic"], 
                                       gg$Particles[gg$Model == "Truth"],
                                       method = "pearson")$estimate[[1]],
                                cor.test(gg$Particles[gg$Model == "Gompertz"], 
                                         gg$Particles[gg$Model == "Truth"],
                                         method = "pearson")$estimate[[1]]))
                
                
# derivate
data[['g']] = data$time/20
derive = data.frame(Diff = diff(data$Gompertz) /diff(data$g), g = data$g[-1])
ggplot(derive)+
  aes(x = g, y = Diff) +
  geom_point()
                
                
                
