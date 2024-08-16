#
# Ajuste con fitdistrplus
#
require(fitdistrplus)
require(actuar)
plotdist(da_G$biomass, histo = TRUE, demp = TRUE)
fitW <- fitdist(da_N$biomass, "pareto", start = list(scale =11, shape = 1.2))
fitln <- fitdist(da_N$biomass, "lnorm")
fitg <- fitdist(da_N$biomass, "gamma")
denscomp(list(fitW, fitln, fitg),xlogscale=TRUE,ylogscale=TRUE)
gofstat(list(fitW, fitln, fitg), fitnames = c("Pareto", "lnorm", "Gamma"))

dcomp <- cdfcomp(list(fitW, fitln, fitg), legendtext = c("Pareto", "lognormal", "gamma"), xlogscale = TRUE, ylogscale = TRUE,
    xlab = "body biomass", 
    fitcol = c("red", "green", "orange"), fitlty = 1, fitlwd = .5, 
    xlegend = "topright", plotstyle = "ggplot", addlegend = FALSE)
dcomp + ggplot2::theme_minimal() 
