#' @title Miscellaneous supporting functions specific to this package
#'
#' @description  NA
#'
#' @name helper_functions
#' @author Brandon McNellis
NULL
#' @rdname helper_functions
#' @export
GetFilename <- function(BITraster, var, year) {

  stopifnot(
    validObject(BITraster),
    var %in% BITraster@variables,
    year %in% BITraster@years
  )

  l1 <- list.files(BITraster@data_directory, pattern = paste0('^', var))
  l2 <- grep(EROSBIT::GetYearSuffix(year), x = l1, value = T)
  l3 <- l2[dir.exists(paste0(BITraster@data_directory, '/', l2))]
  l4 <- paste0(BITraster@data_directory, '/', l3)
  l5 <- list.files(l4, pattern = '.img$', full.names = T)
  l6 <- grep(pattern = paste0(year, '_mos'), x = l5, value = T)

  return(l6)

}
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
GetYearSuffix <- function(year) {

  if (year %in% c(1985:1995)) {
    return('1985_1995')
  } else if (year %in% c(1996:2006)) {
    return('1996_2006')
  } else if (year %in% c(2007:2018)) {
    return('2007:2018')
  } else {
    stop('bad year input, out of range of EROSBIT')
  }

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
#' @rdname helper_functions
#' @export
CheckMaskValidity <- function(BITraster) {
  stopifnot(inherits(BITraster, 'BITraster'))

  if (length(BITraster@mask) == 0) {
    return(TRUE)
  } else {
    require(raster)
    require(sp)

    mask0 <- raster::raster(BITraster@mask)

    t0 <- raster::compareRaster(BITraster, mask0)
    return(t0)

  }

  stop("failed to check mask validity")
}
