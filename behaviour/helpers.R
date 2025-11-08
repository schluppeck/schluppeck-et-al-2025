# ---------------------------------
# from stat_ellipse code! ggplot
# broken out, so we can inspect maths / scaling, etc.
# ds, 2025-01-15
# ---------------------------------

calculate_ellipse <- function(data, vars, type = "t", level = 0.95, segments = 51){
  dfn <- 2
  dfd <- nrow(data) - 1

  if (!type %in% c("t", "norm", "euclid")) {
    cli::cli_inform("Unrecognized ellipse type")
    ellipse <- matrix(NA_real_, ncol = 2)
  } else if (dfd < 3) {
    cli::cli_inform("Too few points to calculate an ellipse")
    ellipse <- matrix(NA_real_, ncol = 2)
  } else {
    if (type == "t") {
      v <- MASS::cov.trob(data[,vars])
    } else if (type == "norm") {
      v <- stats::cov.wt(data[,vars])
    } else if (type == "euclid") {
      v <- stats::cov.wt(data[,vars])
      v$cov <- diag(rep(min(diag(v$cov)), 2))
    }
    shape <- v$cov
    center <- v$center
    chol_decomp <- chol(shape)
    if (type == "euclid") {
      radius <- level/max(chol_decomp)
    } else {
      radius <- sqrt(dfn * stats::qf(level, dfn, dfd))
    }
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    ellipse <- t(center + radius * t(unit.circle %*% chol_decomp))
  }

  colnames(ellipse) <- vars
  mat_2_df(ellipse)
}


# make a function that takes a CSV file and returns a dataframe

#' read data and parse amblyope, control + subject
#'
#' @param fname A valid CSV file name
#' @returns A tidy dataframe with behavioural data
#' @examples
#' readPlus("./fmri-behaviour/AbmDichoptic.csv")
readPlus <- function(fname) {
  read_csv(fname) |>
    mutate(filename = basename(fname)) |>
    mutate(kind = as_factor(str_extract(fname, "[AC]"))) |>
    mutate(which_eye = as_factor(str_to_lower(
      str_extract(filename, ".*(Good|Bad).csv", group = 1)
    ))) |>
    mutate(sub = as_factor(str_extract(fname, ".*[AC]([a-z]*).*", group = 1))) |>
    select(-filename)
}


readPlusDichoptic <- function(fname) {
  read_csv(fname) |>
    mutate(filename = basename(fname)) |>
    mutate(kind = as_factor(str_extract(fname, "[AC]"))) |>
    mutate(which_eye = as_factor(str_to_lower(
      str_extract(filename, ".*(Good|Bad|Dichoptic).csv", group = 1)
    ))) |>
    mutate(sub = as_factor(str_extract(fname, ".*[AC]([a-z]*).*", group = 1))) |>
    select(-filename)
}


