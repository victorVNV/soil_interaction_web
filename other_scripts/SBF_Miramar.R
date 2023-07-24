
#devtools::install_github('lsaravia/EcoIndicators')


library(tidyverse)
library(readxl)
library(ade4)
library(vegan)
library(cluster)
library(BiodiversityR)
library(Ecoindicators)

miramar <- read_excel("soil_networks.xlsx", sheet = "SBF_Miramar")

spc_miramar<-miramar[,1]

t_miramar<-t(miramar[,-1])

dimnames(t_miramar)[2]<-spc_miramar

matriz_similitud<- vegdist(t_miramar,method="bray")
cluster_gradiente <- hclust(matriz_similitud,method = "average")
plot(cluster_gradiente, main="Similitud de los embudos de Berleese de las réplicas en Miramar", sub="Coeficiente Bray Curtis", xlab = "Réplicas en Miramar")

indice_simpson<-diversity(t_miramar, index="simpson")
plot(indice_simpson)

indice_shanon<-diversity(t_miramar, index="shannon")
plot(indice_shanon)

riqueza<-specnumber(t_miramar)
plot(riqueza)

ind_wilson<-betadiver(t_miramar, method=8)
cluster_Beta_Wilson <- hclust(ind_wilson,method = "average")
plot(cluster_Beta_Wilson)

whitaker_ANT<-rankabundance(t_miramar[c(1:5),])
grafico_whitakerANT<-rankabunplot(whitaker_ANT, scale='proportion', addit=FALSE, specnames=c(1:15),srt=90,ylim=c(0, 5), pch=15,col= c(1:15))

whitaker_NAT<-rankabundance(t_miramar[c(6:10),])

grafico_whitakerNAT<-rankabunplot(whitaker_NAT, scale='proportion', addit=FALSE, specnames=c(1:15),srt=90,ylim=c(0, 1.5), pch=15,col= c(1:15))

#Prueba de especies indicadoras
  #Prueba con réplicas

indicator_species_miramar<-Ecoindicators::select_indicator_species(com=as.data.frame(t_miramar),group = as.vector(row.names(t_miramar)))
indval_miramar<-labdsv::indval(t_miramar,clustering = as.vector(row.names(t_miramar)))
summary(indval_miramar)
#No se hallaron taxas indicadoras del sitio.
#Ok, si cambio alfa, no cambian los resultados.
#Si quito los animales que dominan, o no conformn el grupo, tambíen quedan los mismos resultados

  #Prueba por sitio

ANT<-colSums(t_miramar[c(1:5),])
NAT<-colSums(t_miramar[c(6:10),])

t_miramar2<-rbind(ANT,NAT)

indicator_species_miramar2<-Ecoindicators::select_indicator_species(com=as.data.frame(t_miramar2),group = as.vector(row.names(t_miramar2)))
indval_miramar2<-labdsv::indval(t_miramar2,clustering = as.vector(row.names(t_miramar2)))
summary(indval_miramar2)

summary(indice_simpson[c(1:5)])
summary(indice_simpson[c(6:10)])

mean(indice_simpson[c(1:5)])
mean(indice_simpson[c(6:10)])

sd(indice_simpson[c(1:5)])
sd(indice_simpson[c(6:10)])

###################

library(boot)

# Función para calcular la media
mean_func <- function(data, indices) {
  mean(data[indices])
}

# Función para calcular el desvío estándar
sd_func <- function(data, indices) {
  sd(data[indices])
}

# Número de repeticiones de bootstrap
n_boot <- 1000

# Realizar bootstrap para calcular la media y el desvío estándar
boot_mean_A <- boot(indice_simpson[c(1:5)], mean_func, R = n_boot)
boot_sd_A <- boot(indice_simpson[c(1:5)], sd_func, R = n_boot)

# Obtener resultados
bootstrap_mean_A <- boot.ci(boot_mean_A, type = "basic")$basic
bootstrap_sd_A <- boot.ci(boot_sd_A, type = "basic")$basic

# Imprimir resultados
print("Bootstrap Mean:")
print(bootstrap_mean_A)
print("Bootstrap Standard Deviation:")
print(bootstrap_sd_A)

# Realizar bootstrap para calcular la media y el desvío estándar
boot_mean_N<- boot(indice_simpson[c(6:10)], mean_func, R = n_boot)
boot_sd_N <- boot(indice_simpson[c(6:10)], sd_func, R = n_boot)

# Obtener resultados
bootstrap_mean_N <- boot.ci(boot_mean_N, type = "basic")$basic
bootstrap_sd_N <- boot.ci(boot_sd_N, type = "basic")$basic

# Imprimir resultados
print("Bootstrap Mean:")
print(bootstrap_mean_N)
print("Bootstrap Standard Deviation:")
print(bootstrap_sd_N)
