#' @title S4 class for extending Raster* to EROS BIT dataset
#'
#' @description NA
#'
#' @details NA
#'
#' @include generics.R
#' @import sp
#' @import raster
#'
#' @author Brandon McNellis
#' @name BITraster
#' @rdname BITraster
NULL
#' An S4 class for extending Raster* to EROS BIT
#'
#' @slot variables Character vector describing the contained BIT variables
#' @slot data_directory Full path to un-zipped BIT data, defaults to current working directory
#' @slot temp_directory Optional place to put temporary raster files (they get big)
#' @slot years Years of the analysis
#' @slot character Character vector with filepath to a mask. See details
#'
#' @details
#'
#' `mask` must be a RasterLayer with the same extent and resolution as the BITraster
#'
#' @rdname BITraster
#' @export
BITraster <- setClass(
  'BITraster',
  contains = "RasterBrick",
  slots = list(
    variables = "character",
    data_directory = "character",
    temp_directory = "character",
    years = "numeric",
    mask = "character"
  )
)
#' @export
setValidity('BITraster', function(object) {
  require(raster)

  errors <- character()

  if (!dir.exists(object@data_directory)) {
    msg <- paste0('data directory does not exist')
    errors <- c(errors, msg)
  }
  #if (!dir.exists(object@temp_directory)) {
  #  msg <- paste0('temp directory does not exist or has been moved')
  #  errors <- c(errors, msg)
  #}
  if (!all(object@variables %in% VarNames())) {
    msg <- paste0('malformed variable names')
    errors <- c(errors, msg)
  }
  if (!all(object@years %in% AllYears())) {
    msg <- paste0('invalid year range')
    errors <- c(errors, msg)
  }

  # check origins
  orig1 <- raster::origin(object)
  orig1 <- paste0(orig1[1], ', ', orig1[2])
  orig2 <- raster::origin(EROSBIT::sample_BIT)
  orig2 <- paste0(orig2[1], ', ', orig2[2])

  if (!identical(orig1, orig2)) {
    msg <- paste0('bad origin, got ', orig1, ' should be ', orig2)
    errors <- c(errors, msg)
  }

  # check resolution
  res1 <- raster::res(object)
  res1 <- paste0(res1[1], ', ', res1[2])
  res2 <- raster::res(EROSBIT::sample_BIT)
  res2 <- paste0(res2[1], ', ', res2[2])

  if (!identical(res1, res2)) {
    msg <- paste0('bad resolution, got ', res1, ' should be ', res2)
    errors <- c(errors, msg)
  }

  mask <- object@mask
  if (length(mask) > 1) {
    msg <- paste0('mask must be length-1 character')
    errors <- c(errors, msg)
  } else if (length(mask) == 1) {
    if (!file.exists(mask)) {
      msg <- paste0('mask file does not exist')
      errors <- c(errors, msg)
    }
  }

  # returns
  if (length(errors) == 0) {
    TRUE
  } else {
    errors
  }
})
#' @rdname BITraster
#' @export
setMethod('initialize',
          signature(.Object = 'BITraster'),
          function (.Object, ...) {
            require(raster)
            params <- list(...)

            if ('temp_directory' %in% names(params)) {
              .Object@temp_directory <- params$temp_directory
            } else {
              .Object@temp_directory <- tempdir()
            }

            if ('data_directory' %in% names(params)) {
              .Object@data_directory <- params$data_directory
            } else {
              .Object@data_directory <- getwd()
            }

            # BITraster defaults
            .Object@crs <- EROSBIT::EROSBIT_CRS
            .Object@mask <- character()
            # Raster* defaults
            raster::res(.Object) <- raster::res(EROSBIT::sample_BIT)
            raster::crs(.Object) <- raster::crs(EROSBIT::sample_BIT)
            # origin MUST come last here, for some reason...
            raster::origin(.Object) <- raster::origin(EROSBIT::sample_BIT)

            # returns
            .Object <- callNextMethod()
            mt <- validObject(.Object)
            if (isTRUE(mt)) {
              return(.Object)
            } else {
              return(mt)
            }
          }
)
setAs('RasterLayer', 'BITraster', def = function(from) {
  require(raster)
  to <- as(from, 'RasterBrick')
  to <- as(to, 'BITraster')
  to
})
setAs('RasterStack', 'BITraster', def = function(from) {
  require(raster)
  to <- as(from, 'RasterBrick')
  to <- as(to, 'BITraster')
  to
})
#' @rdname BITraster
#' @export
.make_mask <- function(BITraster, mask, filename) {
  stopifnot(validObject(BITraster))
  require(raster)

  raster::values(mask) <- ifelse(is.na(raster::values(mask)), NA, 1)
  mask <- raster::resample(mask, BITraster, 'ngb')
  raster::writeRaster(x = mask, filename = filename, format = 'GTiff', overwrite = T)

  BITraster@mask <- filename
  stopifnot(EROSBIT::CheckMaskValidity(BITraster), validObject(BITraster))
  return(BITraster)

}
#' @rdname BITraster
#' @export
setMethod('MakeMask',
          signature(mask = 'SpatialPolygonsDataFrame'),
          function(BITraster, mask, filename) {
            require(raster)
            require(sp)

            # change mask to RasterLayer format
            mask0 <- sp::spTransform(mask, EROSBIT::EROSBIT_CRS)

            # define target raster
            mask_blank <- raster::raster(mask0)
            # set resolution to 30m x 30m
            raster::res(mask_blank) <- raster::res(EROSBIT::sample_BIT)
            # adjust origin
            s0 <- raster::origin(EROSBIT::sample_BIT) - raster::origin(mask_blank)
            mask_blank <- raster::shift(mask_blank, dx = s0[1], dy = s0[2])
            # convert template to raster
            mask0 <- raster::rasterize(mask0, mask_blank)

            .make_mask(
              BITraster = BITraster,
              mask = mask0,
              filename = filename
            )
          })
#' @rdname BITraster
#' @export
.template_BIT <- function(template, dir, variables, years) {
  require(raster) # implicit

  # Input checks:
  stopifnot(
    dir.exists(dir),
    is.numeric(years),
    is.character(variables),
    all(variables %in% EROSBIT::VarNames()),
    all(years %in% EROSBIT::AllYears()),
    inherits(template, 'RasterLayer')
  )

  # Reset origin, pad the extent, and resample to new raster
  blank_template <- template
  values(blank_template) <- NA
  s0 <- origin(EROSBIT::sample_BIT) - origin(blank_template)
  blank_template <- shift(blank_template, dx = s0[1], dy = s0[2])
  blank_template <- extend(blank_template, 35) # ca. 1km if res = 30, 30
  m0 <- ifelse(template@data@isfactor, 'ngb', 'bilinear')
  #attr0 <- template@data@attributes # preserve attributes
  template <- resample(template, blank_template, m0)
  #template@data@attributes <- attr0

  # Convert to BITraster and add meta-data
  out_BIT <- as(template, 'BITraster')
  out_BIT@data_directory <- dir
  out_BIT@variables <- variables
  out_BIT@years <- years

  # Output checks:

  compareRaster(out_BIT, template, res = T, orig = T, stopiffalse = T)

  gc()
  return(out_BIT)

}
#' @rdname BITraster
#' @export
setMethod('TemplateBIT',
          signature(template = 'SpatialPolygonsDataFrame'),
          function(template, dir, variables, years) {
            require(raster)
            require(sp)

            # TODO: parallelize the rasterize? seems slow

            # Conver SPDF to RasterLayer

            # convert CRS
            template_new <- sp::spTransform(template, crs(EROSBIT::sample_BIT))

            # define target raster
            template_rast <- raster::raster(template_new)
            # set resolution to 30m x 30m
            raster::res(template_rast) <- raster::res(EROSBIT::sample_BIT)
            # convert template to raster
            template_rast <- raster::rasterize(template_new, template_rast)

            .template_BIT(
              template = template_rast,
              dir = dir,
              variables = variables,
              years = years
            )
          })
#' @rdname BITraster
#' @export
setMethod('TemplateBIT',
          signature(template = 'RasterLayer'),
          function(template, dir, variables, years) {
            .template_BIT(
              template = template,
              dir = dir,
              variables = variables,
              years = years
            )
          })
#' @rdname BITraster
#' @export
setMethod('addLayer',
          signature(x = 'BITraster'),
          function(x, ...) {
            require(raster)
            stopifnot(validObject(x))

            y <- callNextMethod()
            y <- as(y, 'BITraster')

            y@variables <- x@variables
            y@data_directory <- x@data_directory
            y@temp_directory <- x@temp_directory
            y@years <- x@years
            y@mask <- x@mask

            return(y)

          })
#' @rdname BITraster
#' @export
.subset_BIT_points <- function(BITraster, points, variables, years, parallel = FALSE) {

  require(sp)
  require(raster)

  # This function subsets 30-m EROSBIT pixels from a SpatialPoints object
  # Adapted from FDF project code

  # Returns a SpatialPointsDataFrame
  # Future options to return a raster?

  # This could be changed to a method for BITraster if need be

  # sanity checks
  stopifnot(
    inherits(BITraster, 'BITraster'),
    #inherits(points, 'SpatialPoints'),
    is.character(variables),
    is.numeric(years),
    validObject(BITraster),
    all(years %in% BITraster@years),
    all(variables %in% BITraster@variables)
  )

  # main loop
  out_df <- points[, 2, drop = FALSE]

  if (parallel) {
    # implement a foreach here with foreach package
    require(foreach)
    stop('not implemented')

  } else if (!parallel) {

    for (i in seq_along(years)) {
      ii <- years[i]

      ilayer <- subset(BITraster, which(BITraster@years == ii))
      ivals <- extract(ilayer, points[, 1])
      out_df <- cbind(out_df, ivals)
      colnames(out_df)[ncol(out_df)] <- as.character(ii)

    }

    return(out_df)

  }

}
#' @rdname BITraster
#' @export
setMethod('SubsetPoints',
          signature(points = 'matrix'),
          function(BITraster, points, variables, years, parallel = FALSE) {
            .subset_BIT_points(
              BITraster = BITraster,
              points = points,
              variables = variables,
              years = years,
              parallel = parallel)

          })
#' @rdname BITraster
#' @export
setMethod('SubsetPoints',
          signature(points = 'data.frame'),
          function(BITraster, points, variables, years, parallel = FALSE) {
            .subset_BIT_points(
              BITraster = BITraster,
              points = as.data.frame(points),
              variables = variables,
              years = years,
              parallel = parallel)

          })
#' @rdname BITraster
#' @export
setMethod('SubsetPoints',
          signature(points = 'SpatialPoints'),
          function(BITraster, points, variables, years, parallel = FALSE) {

            require(sp)
            # xyFromCell from the parent BITraster gets coordinates

            stopifnot(
              identical(crs(BITraster) == crs(points)),
              identical(bbox(BITraster), bbox(points))
            )

            points0 <- BIT_spatial@coords

            .subset_BIT_points(
              BITraster = BITraster,
              points = points0,
              variables = variables,
              years = years,
              parallel = parallel)

          })
#' @rdname BITraster
#' @export
setMethod('SubsetPoints',
          signature(points = 'SpatialPolygonsDataFrame'),
          function(BITraster, points, variables, years, parallel = FALSE) {

            require(sp)
            # xyFromCell from the parent BITraster gets coordinates

            stopifnot(
              identical(crs(BITraster) == crs(points)),
              identical(bbox(BITraster), bbox(points))
            )

            stop('not implemented')

            .subset_BIT_points(
              BITraster = BITraster,
              points = points0,
              variables = variables,
              years = years,
              parallel = parallel)

          })
