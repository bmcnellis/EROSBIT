
desert <- "CH"

library(raster)

beginCluster()
all_files <- list.files(pattern = ".img")
img_1 <- raster(all_files[1])
desert_mask <- raster(paste0(desert, "_mask_30m.grd"))
desert_mask <- projectRaster(desert_mask, crs = crs(img_1), method = "ngb")

for (i in seq_along(all_files)) {
    ifl <- all_files[i]
    inm <- paste0(unlist(strsplit(ifl, ".img")), "_masked_", desert, ".grd")
    i_rast <- raster(ifl)
    i_rast <- crop(i_rast, desert_mask)

    extent(desert_mask)@xmin <- ifelse(abs(extent(i_rast)@xmin - extent(desert_mask)@xmin) < 10, extent(i_rast)@xmin, extent(desert_mask)@xmin)
    extent(desert_mask)@xmax <- ifelse(abs(extent(i_rast)@xmax - extent(desert_mask)@xmax) < 10, extent(i_rast)@xmax, extent(desert_mask)@xmax)
    extent(desert_mask)@ymin <- ifelse(abs(extent(i_rast)@ymin - extent(desert_mask)@ymin) < 10, extent(i_rast)@ymin, extent(desert_mask)@ymin)
    extent(desert_mask)@ymax <- ifelse(abs(extent(i_rast)@ymax - extent(desert_mask)@ymax) < 10, extent(i_rast)@ymax, extent(desert_mask)@ymax)

    stopifnot(all.equal(extent(i_rast), extent(desert_mask)))
    i_rast <- calc(i_rast, function(x) ifelse(x > 100, NA, x))
    i_rast <- mask(i_rast, desert_mask, maskvalue = NA)
    writeRaster(i_rast, file = inm, format = "raster", datatype = "FLT4S", overwrite = T)
}

endCluster()


