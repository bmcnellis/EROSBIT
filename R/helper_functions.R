#' @title Miscellaneous supporting functions specific to this package
#'
#' @description  NA
#'
#' @name helper_functions
#' @author Brandon McNellis
NULL
#' @rdname helper_functions
#' @export
VarNames <- function() {

  vars <- c('Annual_Herb', 'Herb', 'Sage', 'Shrub', 'Bare', 'Litter')
  return(vars)
}
#' @rdname helper_functions
#' @export
AllYears <- function() {

  yrs <- c(1985:2018)
  return(yrs)
}
#' @rdname helper_functions
#' @export
matchResolution <- function(x,
                            ref,
                            method="bilinear",
                            filename="",
                            ...){

  # matchResolution is from mqueinnec/foster
  require(raster)

  if (!class(x)[1] %in% c("RasterLayer", "RasterBrick", "RasterStack")) {
    stop("x must be a Raster object")
  }

  if (!class(ref)[1] %in% c("RasterLayer", "RasterBrick", "RasterStack")) {
    stop("ref must be a Raster object")
  }

  #Check CRS
  if (is.na(raster::crs(x)) | is.na(raster::crs(ref))) {
    stop("CRS of x or ref is not defined")
  } else if (!raster::compareCRS(crs(x), crs(ref))) {
    warning("x and ref don't have the same CRS. x is projected to ref CRS before
            resampling")
    x <- raster::projectRaster(x, crs = raster::crs(ref))
  }

  if (raster::extent(ref) > raster::extent(x)) {
    #We crop ref to x extent. It avoids creating a large resampled x if ref
    #extent is much larger than x
    ref_crop <- raster::crop(ref, x, filename = "")
  } else {
    ref_crop <- ref
  }

  #Resampling
  out <- raster::resample(x = x, y = ref_crop, method = method, filename =
                            filename, ...)

  return(out)
}
