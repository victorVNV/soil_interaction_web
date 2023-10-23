
#Carga paquetes
{library(tidyverse);library(readxl);library(ade4)
library(vegan);library(cluster);library(BiodiversityR)}

#Anális de datos de ROSANA
#Clásicos multivariados

comunity <- read_excel("soil_networks.xlsx", sheet = "Ro_web")

resumen_bruto<-summary(comunity);resumen_bruto

#NOT FILTER, macrofauna, y anélidos

full_db<-comunity %>% group_by(system,site,year,season,
                                   order,guild,specie) %>%
  summarise(abundance=n(), peso.ind=mean(biomass), 
            bmicro=(mean(bm)), POM=mean(POM),
            res.mec=mean(RM5),DAP=mean(DAP),HR=mean(HR),
            pH=mean(Ph),CE=mean(CE),
            MO=mean(MO),N=mean(N),P=mean(P),K=mean(K),
            Ca=mean(Ca),Mg=mean(Mg),Na=mean(Sodio))%>%
  distinct(specie,.keep_all=TRUE)

#POM, correlacionado con MO, bm biomasa microbiana correlacionada con Respiración

#Generar matriz de abundancia, peso corporal, Biomasa
# de especies en: estación x año x sitio x uso

full_ab_temporal<-full_db %>% 
  pivot_wider(names_from = specie, 
              values_from = abundance, 
              values_fill = 0)

full_peso_temporal<-full_db %>% 
  pivot_wider(names_from = specie, 
              values_from = peso.ind, 
              values_fill = 0)

full_B_temporalx<-full_db %>% 
  mutate(B=(peso.ind*abundance))%>%
  pivot_wider(names_from = specie, 
              values_from = B, 
              values_fill = 0)

#Generar matriz de abundancia, peso corporal, Biomasa
# de especies en: año x sitio x uso

year_db<-comunity %>% group_by(system,site,year,
                               order,guild,specie) %>%
  summarise(abundance=n(), peso.ind=mean(biomass), 
            bmicro=(mean(bm)), POM=mean(POM),
            res.mec=mean(RM5),DAP=mean(DAP),HR=mean(HR),
            pH=mean(Ph),CE=mean(CE),
            MO=mean(MO),N=mean(N),P=mean(P),K=mean(K),
            Ca=mean(Ca),Mg=mean(Mg),Na=mean(Sodio))%>%
  distinct(specie,.keep_all=TRUE)

year_ab_temporal<-year_db %>% 
  pivot_wider(names_from = specie, 
              values_from = abundance, 
              values_fill = 0)

year_peso_temporal<-year_db %>% 
  pivot_wider(names_from = specie, 
              values_from = peso.ind, 
              values_fill = 0)

year_B_temporalx<-year_db %>% 
  mutate(B=(peso.ind*abundance))%>%
  pivot_wider(names_from = specie, 
              values_from = B, 
              values_fill = 0)

#Generar matriz de abundancia, peso corporal, Biomasa
# de especies en: sitio x uso

sitio_db<-comunity %>% group_by(system,site,
                               order,guild,specie) %>%
  summarise(abundance=n(), peso.ind=mean(biomass), 
            bmicro=(mean(bm)), POM=mean(POM),
            res.mec=mean(RM5),DAP=mean(DAP),HR=mean(HR),
            pH=mean(Ph),CE=mean(CE),
            MO=mean(MO),N=mean(N),P=mean(P),K=mean(K),
            Ca=mean(Ca),Mg=mean(Mg),Na=mean(Sodio))%>%
  distinct(specie,.keep_all=TRUE)

sitio_ab_temporal<-sitio_db %>% 
  pivot_wider(names_from = specie, 
              values_from = abundance, 
              values_fill = 0)

sitio_peso_temporal<-sitio_db %>% 
  pivot_wider(names_from = specie, 
              values_from = peso.ind, 
              values_fill = 0)

sitio_B_temporalx<-sitio_db %>% 
  mutate(B=(peso.ind*abundance))%>%
  pivot_wider(names_from = specie, 
              values_from = B, 
              values_fill = 0)

#Generar matriz de abundancia, peso corporal, Biomasa
# de especies en: x uso

uso_db<-comunity %>% group_by(system,
                                order,guild,specie) %>%
  summarise(abundance=n(), peso.ind=mean(biomass), 
            bmicro=(mean(bm)), POM=mean(POM),
            res.mec=mean(RM5),DAP=mean(DAP),HR=mean(HR),
            pH=mean(Ph),CE=mean(CE),
            MO=mean(MO),N=mean(N),P=mean(P),K=mean(K),
            Ca=mean(Ca),Mg=mean(Mg),Na=mean(Sodio))%>%
  distinct(specie,.keep_all=TRUE)

uso_ab_temporal<-uso_db %>% 
  pivot_wider(names_from = specie, 
              values_from = abundance, 
              values_fill = 0)

uso_peso_temporal<-uso_db %>% 
  pivot_wider(names_from = specie, 
              values_from = peso.ind, 
              values_fill = 0)

uso_B_temporalx<-uso_db %>% 
  mutate(B=(peso.ind*abundance))%>%
  pivot_wider(names_from = specie, 
              values_from = B, 
              values_fill = 0)
