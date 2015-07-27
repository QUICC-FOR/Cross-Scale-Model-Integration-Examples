library(rgdal)
library(sp)

P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
P4S.albers = CRS("+init=epsg:2150")
P4S.albers = CRS("+init=epsg:3087")

ocean = readOGR(dsn="dat/figure/ne_50m_ocean", layer="ne_50m_ocean")
lakes = readOGR(dsn="dat/figure/ne_50m_lakes", layer="ne_50m_lakes")
lkNames = c("Huron", "Michigan", "Superior", "Ontario", "Erie", "St. Clair")
grLakes = lakes[as.integer(sapply(lkNames, grep, lakes$name)),]

presDat = readRDS("dat/raw/AceSac_pres.rds")
coordinates(presDat) = c('lon', 'lat')
proj4string(presDat) = P4S.latlon
presDat.albers = spTransform(presDat, P4S.albers)
ocean.albers = spTransform(ocean, P4S.albers)
grLakes.albers = spTransform(grLakes, P4S.albers)

# par(mfrow=c(1,2))
# plot(presDat, pch=4, cex=0.8)
# plot(ocean, col="white", add=T)
# plot(grLakes, col="white", add=T)
# 
pdf(file="plotmap.pdf", width=4.5, height=4.5)
par(mar=c(0.1,0.1,0.1,0.1))
plot(presDat.albers, pch=4, cex=0.1, axes=F, ylim=c(400000, 3500000), bty='o')
plot(ocean.albers, col="white", add=T)
plot(grLakes.albers, col="white", add=T)
dev.off()