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
    return('2007_2018')
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
#' @rdname helper_functions
#' @export
ExpandToMask <- function(BITraster) {

  require(raster)
  require(rgdal)
  stopifnot(validObject(BITraster))

  n0 <- parallel::detectCores() - 1
  s0 <- 'snow' %in% installed.packages()

  if (n0 > 0 && s0) {
    beginCluster(n0)
  }

  # Load the mask raster
  mask00 <- raster(BITraster@mask)

  # Reproject mask raster to BITraster CRS
  mask0 <- projectRaster(
    mask00,
    res = c(30, 30),
    crs = crs(BITraster),
    method = 'ngb'
  )

  # Shift mask raster origin
  dxdy <- 15 - origin(mask0)
  mask0 <- raster::shift(mask0, dxdy[1], dxdy[2])
  res(mask0) <- round(res(mask0), 4)

  # Convert mask to BITraster
  expanded_BIT <- as(mask0, 'BITraster')
  # fix resolution rounding - needs to happen AFTER coercion
  res(expanded_BIT) <- round(res(expanded_BIT), 4)

  if (validObject(expanded_BIT)) {
    f1 <- tools::file_path_sans_ext(BITraster@mask)
    f2 <- tools::file_ext(BITraster@mask)

    new_mask_fname <- paste0(f1, '_extended.', f2)
    if (f2 == 'tif') {
      new_mask_format <- 'GTiff'
    } else {
      stop('dont know the mask file format')
    }

    writeRaster(mask0, filename = new_mask_fname, format = new_mask_format, overwrite = T)
    expanded_BIT@mask <- new_mask_fname

  } else {
    stop('expanding function failed')
  }

  # Transfer other slots
  expanded_BIT@variables <- BITraster@variables
  expanded_BIT@data_directory <- BITraster@data_directory
  expanded_BIT@temp_directory <- BITraster@temp_directory
  expanded_BIT@years <- BITraster@years

  stopifnot(
    validObject(expanded_BIT),
    CheckMaskValidity(expanded_BIT)
  )

  if (n0 > 0 && s0) {
    endCluster()
  }

  return(expanded_BIT)

}
