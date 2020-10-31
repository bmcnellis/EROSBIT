#' @title Collects BIT data specified by a BITraster-class object
#'
#' @description  NA
#'
#' @param BITraster Object of class BITraster
#'
#' @name GetBIT
#' @author Brandon McNellis
#' @export
GetBIT <- function(BITraster) {
  require(raster)

  stopifnot(validObject(BITraster))

  out_BIT <- BITraster

  if (length(BITraster@mask) > 0) {
    stopifnot(CheckMaskValidity(BITraster))
    mask0 <- raster::raster(BITraster@mask)
  }

  for (i in seq_along(BITraster@variables)) {
    i_var <- BITraster@variables[i]

    for (j in seq_along(BITraster@years)) {
      j_year <- BITraster@years[j]

      j_fname <- EROSBIT::GetFilename(BITraster, i_var, j_year)
      j_img <- raster::raster(j_fname)

      # apply mask
      if (length(BITraster@mask) > 0) {
        browser()
        j_img <- raster::crop(j_img, mask0)
        raster::compareRaster(j_img, mask0, values = F, orig = TRUE, stopiffalse = T)
        j_img <- raster::mask(j_img, mask0)
      }

      # add to out_BIT
      # addLayer needs a method for BITraster or it will just return a RasterStack
      browser()
      # addLayer(test0, j_img)


    }
  }


}
