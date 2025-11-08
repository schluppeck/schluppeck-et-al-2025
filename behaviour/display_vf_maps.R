# quick look at distortion maps in Amblyopic and Control subjects
#
# this code shows the logic of data aggregation and plotting
#
# ds 2018-06-22
# ds 2025-01-14 picked up again 

## See readme for detailed notes.

library(tidyverse)
library(janitor)
library(cetcolor)
library(cowplot)
library(measurements)
library(broom)
library(ggridges)
library(gtsummary)


# source helper functions
source("./helpers.R")
source("./performance.R")

#good orange, bad blueish
theColors <- c("#ff9922","#25baff")

theme_set(theme_cowplot())

# map it over all files that match the pattern. In this case Amblyopbia data
# - through good eye

theData <- list.files(
  path = "./fmri-behaviour/",
  pattern = "[AC]*(Good|Bad).csv",
  full.names = TRUE
)

# map over a bunch of the CSV files, rearrange, clean names, rename
data <- map_df(theData, ~ readPlus(.)) |>
  select(sub, kind, everything()) |>
  clean_names() |>
  select(sub, kind,
    target_x0 = target_x_512,
    target_y0 = target_y_384,
    response_x0 = response_x_512,
    response_y0 = response_y_384,
    target_orientation_deg,
    target_radius_pixels,
    which_eye
  ) |> 
mutate(sub_num = as.numeric(sub))

glimpse(data)

# clean up column names and also add polar coordinates
# also adding a "match_x, match_y" location, which is on equal and opposite side
# of fixation. for origin-aligned data, negative 1 is to flip x and y coords, which achieves this point symmetry
cart2pol_r <- function(x, y) {
  sqrt(x^2 + y^2)
}

cart2pol_theta <- function(x, y) {
  180 / pi * atan2(y, x)
}

#' convert from pixels to degree
#' 
#'  in the study, the stimuli were placed a 1,3,5,7º ecc
#'  conversion is by this factor
#'  take into account small rounding errors, but only for STIMULUS locations
pix2deg <- function(r) {
  DEG_PER_PIX <- 0.01566
  # return this:
  r * DEG_PER_PIX
}

data <- data |>
  # convert everything to degrees
  mutate(
    target_x0 = pix2deg(target_x0),
    target_y0 = pix2deg(target_y0),
    match_x0 = -target_x0,
    match_y0 = -target_y0,
    response_x0 = pix2deg(response_x0),
    response_y0 = pix2deg(response_y0),
    response_r = cart2pol_r(response_x0, response_y0),
    response_theta = cart2pol_theta(response_x0, response_y0),
    match_r = cart2pol_r(match_x0, match_y0),
    match_theta = cart2pol_theta(match_x0, match_y0),
    target_radius_deg = round(pix2deg(target_radius_pixels))
  )

# make target id... inner ring... 1:6... next ring 7:12 ...
# label targets with ids per individual target, ori(theta) and radius(r)
error_data <- data |>
  mutate(
    target_ori_id = as.numeric(as_factor(target_orientation_deg)),
    target_rad_id = as.numeric(as_factor(target_radius_deg)),
    target_id = target_ori_id + (target_rad_id - 1) * max(target_ori_id)
  ) |>
  group_by(sub, kind, target_id) |>
  mutate(
    error_x = response_x0 - match_x0,
    error_y = response_y0 - match_y0,
    error_r = response_r - match_r,
    error_theta = response_theta - match_theta,
    error_r_mu = mean(error_r),
    error_r_sd = sd(error_r)
  )


error_data_lm <- error_data |>
  mutate(target_ori_id = as_factor(target_ori_id),
)

# pre-JB comments
the_formula = error_r ~ kind + match_r + which_eye  + match_theta

# with interaction
the_formula = error_r ~ kind * match_r + which_eye  + match_theta

m <- lm(the_formula , error_data)

summary(m)

d <- tidy(m, conf.int = TRUE)

ggplot(d, aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0)) +
  geom_point() +
  geom_vline(xintercept = 0, lty = 4) +
  geom_errorbarh() +
  theme(aspect.ratio = 1)

amb_data <- error_data |> 
  filter(kind == "A") |> 
  filter(sub_num %in% c(4))

control_data <- error_data |> 
  filter(kind == "C") |> 
  filter(sub_num %in% c(18))

# example plot of responses that subjects gave
  
plot_distortion_data <- function(data){
  data |> 
      ggplot(aes(
    x = response_x0,
    y = response_y0,
    colour = target_rad_id,
    group = target_id
  )) +
  geom_point(aes(x = target_x0, y = target_y0), color = "black", size = 0.5) +
  stat_ellipse() +
  geom_point(aes(x = response_x0, y = response_y0), alpha = 1, size = 0.1) +
  coord_fixed() +
  scale_colour_gradientn(colours = cet_pal(8, name = "rainbow")) +
  guides(color="none") +
  facet_wrap(~sub) +
  theme_minimal() + 
  scale_x_continuous(breaks = c(-7, 0, 7), labels = c(-7, 0, 7)) +  
  scale_y_continuous(breaks = c(-7, 0, 7), labels = c(-7, 0, 7)) +  
    theme(#axis.text.x = element_blank(), 
        #axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(x = "x Position (deg)", y = "y Position (deg)")
}

BEH_examples <-  error_data |> 
  filter(sub_num %in% c(4,18)) |> 
  plot_distortion_data() 


BEH_examples


#---------- do some calculations

# Group data and fit Gaussian distribution to error_theta
gaussian_fits <- error_data |>
  group_by(kind, which_eye, match_r, match_theta) |>
  summarise(
    fit_r = list(broom::tidy(MASS::fitdistr(error_r, "normal"))),
    fit_theta = list(broom::tidy(MASS::fitdistr(error_theta, "normal")))) |>
  unnest(c(fit_r, fit_theta), names_sep = "_")

# View the results
glimpse(gaussian_fits)


#' function for ploting eccentricity error as a ridge plot
#'
#' @param d The dataframe with error data
#' 
#' 
make_error_ridge_plot <- function(d) {
  d |> 
  ggplot(aes(x = error_r, y = match_r, group = match_r)) +
  geom_density_ridges() +
  #geom_function(fun = dnorm, colour = "red") +
  geom_vline(xintercept = 0, color="red", alpha = 0.5) +
  facet_wrap(~which_eye) +
  theme_minimal() +
  #theme(aspect.ratio = 1.618) +
  scale_x_continuous(limits = c(-2,2)) +
  scale_y_reverse(breaks = c(1, 3, 5, 7),labels = c(1, 3, 5, 7)) +
  guides(colour = "none") +
  labs(x = "Eccentricity error (deg)", y = "Target eccentricity (deg)") 
}


#' function for ploting angular error as a ridge plot
#'
#' @param d The dataframe with error data
#' 
#' 
make_error_ridge_plot_theta <- function(d) {
  d |> 
    select(-kind) |> 
  ggplot(aes(x = error_theta, y = match_r, group=interaction(which_eye, match_r))) +
  geom_density_ridges() +
  geom_vline(xintercept = 0, color="red", alpha = 0.5) +
  facet_wrap(~which_eye) +
  theme_minimal() +
  theme(aspect.ratio = 1.618) +
  scale_x_continuous(limits = c(-25,25)) +
  scale_y_reverse(breaks = c(1, 3, 5, 7), labels = c(1, 3, 5, 7)) +
  guides(colour = "none") +
  labs(x = "Angular error (º)", y = "Target ecentricity (deg)") 
}


# reorg for this plot: AMB / FE / controls
reorg_data <- error_data |>
  group_by(kind) |> 
    mutate(which_eye = if(any(kind == "C")) "control" else which_eye) |> 
      ungroup() |> 
    mutate(which_eye = fct_recode(which_eye, FE="good", AE = "bad")) |> 
    mutate(which_eye = fct_relevel(which_eye, "AE", "FE", "control")) |> 
    ungroup()
  
ECC_plot <- reorg_data |> make_error_ridge_plot() 
RAD_plot <- reorg_data |> make_error_ridge_plot_theta() 

ECC_plot
RAD_plot 
