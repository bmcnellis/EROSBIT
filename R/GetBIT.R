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
  cat('\nstart GetBIT')

  stopifnot(validObject(BITraster))

  out_BIT <- BITraster

  if (length(BITraster@mask) > 0) {
    stopifnot(CheckMaskValidity(BITraster))
    mask0 <- raster::raster(BITraster@mask)
  }

  for (i in seq_along(BITraster@variables)) {
    i_var <- BITraster@variables[i]
    cat('\nvariable:', i_var)

    for (j in seq_along(BITraster@years)) {
      j_year <- BITraster@years[j]
      cat('\nyear:', j_year)

      j_fname <- EROSBIT::GetFilename(BITraster, i_var, j_year)
      j_img <- raster::raster(j_fname)

      # apply mask
      if (length(BITraster@mask) > 0) {
        j_img <- raster::crop(j_img, mask0)
        raster::compareRaster(j_img, mask0, values = F, orig = TRUE, stopiffalse = T)
        j_img <- raster::mask(j_img, mask0)
      }

      # add to out_BIT
      out_BIT <- EROSBIT::addLayer(out_BIT, j_img)
      names(out_BIT)[nlayers(out_BIT)] <- paste0(i_var, '_', j_year)

      if (all(is.na(values(BITraster$layer)))) {
        # Not sure if need this, need to make method for dropLayer if so
        #BITraster <- EROSBIT::dropLayer(BITraster, 1)
      }

    } # end j
  } # end i

  names(out_BIT)[1] <- 'template'

  stopifnot(validObject(out_BIT))
  return(out_BIT)


}
