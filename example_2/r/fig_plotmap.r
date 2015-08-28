library(rgdal)
library(sp)

P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
P4S.albers <- CRS("+init=epsg:5070") # albers equal area conic NAD-83 north america

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
paperwidth = 4.5
dpi = 600
hToWRatio = 1
width = as.integer(dpi*paperwidth)
height = as.integer(width * hToWRatio)
fontsize = 12
png(w=width, h=height, file="plotmap.png", pointsize=fontsize, res = dpi)

par(mar=c(2,2,2,2))
plot(presDat.albers, pch=4, cex=0.1, ylim=c(400000, 3500000), bty='o', xaxt='n', yaxt='n', axes=T)
plot(ocean.albers, col="#9CE3FF", add=T)
plot(grLakes.albers, col="#9CE3FF", add=T)

# plot grid lines and labels
llgridlines(ocean.albers, easts=seq(-100, -70, 10), norths=seq(30,50,10), ndiscr=100, 
		side="ES", offset=1e6, col='#000000FF')
par(xpd=TRUE)
#	westText = c("30°N", "40°N", "50°N")
eastText = c("30°N", "40°N", "50°N")
#	westCoords = list(rep(-3e5,3), c(0.75e6, 1.9e6, 3e6))
eastCoords = list(rep(3.12e6,3), c(1.3e6, 2.45e6, 3.6e6))
southText = c("90°W", "80°W", "70°W")
southCoords = list(c(0.65e6, 1.7e6, 2.8e6), rep(0.2e6,3))
cex.ll = 0.8
text(eastCoords[[1]], eastCoords[[2]], eastText, cex=cex.ll)
#	text(westCoords[[1]], westCoords[[2]], westText, cex=cex.ll)
text(southCoords[[1]], southCoords[[2]], southText, cex=cex.ll)
par(xpd=FALSE)


dev.off()