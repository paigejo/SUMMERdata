require(rgdal)
# Shapefiles obtained from https://gadm.org/download_country_v3.html
adm0 = readOGR(dsn = "~/git/SUMMERdata/gadm36_KEN_shp", layer = "gadm36_KEN_0")
adm1 = readOGR(dsn = "~/git/SUMMERdata/gadm36_KEN_shp", layer = "gadm36_KEN_1")
adm2 = readOGR(dsn = "~/git/SUMMERdata/gadm36_KEN_shp", layer = "gadm36_KEN_2")
kenyaPoly = adm0@polygons[[1]]@Polygons[[77]]@coords

# cluster locations from the 2014 Kenya Demographic Health Survey (KDHS).
library(maptools)
gpsDat = readShapePoints("~/git/U5MR/Kenya2014gps/KEGE71FL.shp")
coords = attr(gpsDat, "coords")
coords = coords[coords[,1] > 20,] # remove bad coordinates set to (0,0)

plot(coords, pch=19, cex=.1, col="blue")
polygon(kenyaPoly)

# project from lat/lon to UTM northing/easting in kilometers.  Use epsg=21097
# either pass lon/east and lat/north, or a matrix with 2 columns: first being lon/east, second being lat/north
# inverse: if FALSE, projects from lon/lat to easting/northing.  Else from easting/northing to lon/lat
library(sp)
projKenya = function(lon, lat=NULL, inverse=FALSE) {
  if(is.null(lat)) {
    lat = lon[,2]
    lon = lon[,1]
  }
  
  if(!inverse) {
    # from lon/lat coords to easting/northing
    # rgdal:::showEPSG("+proj=longlat")
    lonLatCoords = SpatialPoints(cbind(lon, lat), proj4string=CRS("+proj=longlat"))
    coordsUTM = spTransform(lonLatCoords, CRS("+init=epsg:21097 +units=km"))
    out = attr(coordsUTM, "coords")
  }
  else {
    # from easting/northing coords to lon/lat
    east = lon
    north = lat
    coordsUTM = SpatialPoints(cbind(east, north), proj4string=CRS("+init=epsg:21097 +units=km"))
    lonLatCoords = spTransform(coordsUTM, CRS("+proj=longlat"))
    out = attr(lonLatCoords, "coords")
  }
  
  out
}

# project 2014 KDHS data to easting/northing in km
eastNorth = projKenya(coords)

# construct INLA mesh grid for SPDE models
# get a reasonable default mesh triangulation for the SPDE model for the Kenya data
getSPDEMeshKenya = function(locs, n=5000, max.n=5000, doPlot=FALSE, max.edge=c(7, 200), 
                            offset=-.08, cutoff=4, jitterAmount=max.edge[1]/4, seed=123) {
  
  if(is.null(locs)) {
    # jitter the locations used to create the mesh so that they do not always lie on mesh points
    locs=cbind(jitter(mort$east, amount=jitterAmount), jitter(mort$north, amount=jitterAmount))
  }
  
  # generate mesh on R2
  mesh = inla.mesh.2d(loc=locs, n=n, max.n=max.n, offset=offset, cutoff=cutoff, max.edge=max.edge)
  
  mesh
}

# half the size of the maximum shorter edge length. This makes sure peaks of the basis functions 
# are not always at the observation locations (which would underestimate uncertainty elsewhere), 
# but are also more dense near the observations
jitterAmount = 7/2
jitteredEastNorth = cbind(jitter(eastNorth[,1], amount=jitterAmount), jitter(eastNorth[,2], amount=jitterAmount))

library(INLA)
kenyaMesh = inla.mesh.2d(loc=jitteredEastNorth, n=5000, max.n=5000, max.edge=c(7, 200), 
                         offset=-.08, cutoff=4)

save(adm2, adm1, adm0, kenyaPoly, kenyaMesh, file="~/git/SUMMERdata/kenyaMaps.rda")
