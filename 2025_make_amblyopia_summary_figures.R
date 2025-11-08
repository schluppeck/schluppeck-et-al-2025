# code required to make minimal version of stats and figures
# for 2025 Amblyopia pRF paper.
#
# ds 2025-02, jb, zh 2025-03

library(tidyverse)
library(cowplot)
library(janitor)
library(broom)
library(lme4)
library(buildmer)
library(lmerTest)
library(effects)
library(effectsize)
library(officer)
library(flextable)

# source helper functions
source('2025_helpers.R')

# 2025 reprise... control / checks
theme_set(theme_cowplot(font_size=8))

# Options
doFig3 <-1 # group-average radial error
doFigS2 <-1 # group-average angular/tangential error
doFig5 <-1 # group-average pRF estimates
doFigS4 <-1 # individual amblyopes' pRF estimates
doStatsFig3 <- 1
doStatsFigS2 <- 1
doStatsFig5 <- 1
saveFiguresAndTables <- 0 # if 1, save figures without plotting them

# Figure options
theColors <- c(AE = "#25baff", FE = "#ff9922",control="#888888") #FE orange, AE: blueish
theColorsLabels <- c(AE = 'Amblyopic', FE = 'Fellow', control = 'Control')
theColorScale <- scale_colour_manual(name = "Viewing eye",values = theColors, labels = theColorsLabels) # (instead of changing the labels here,
theFillScale <- scale_fill_manual(name = "Viewing eye",values = theColors, labels = theColorsLabels)    # it should be done just after
theLineTypeScale <- scale_linetype_manual(name = "Hemisphere",values=c(RH = "solid", LH = "dashed"),    # just after reading in the data)
                                          labels = c(RH = "Right", LH = "Left"))
pointsSizePRF = 0.2
pointsSizeBehaviour = 0.5
pointsAlpha = 0.25
loessAlpha = 0.1
pointsJitterBehaviour = 0.1
pointsJitterPRF = 0.002
lineWidthBehaviour = 0.6
lineWidthPRF = 0.6

# data file from matlab code / model fitting
fname.behaviour <- "behavioural-data.csv"
fname.brain <- "2025-03-fit-data-summary.csv"
fname.clinical <- "clinical-data.csv"

# participants to keep (this will also defined the participants' number)
amblyopes.initials = c("bm","jh","kr","pb","sm","kg","ii","mp","cg","hb","ra","rc")
controls.initials = c("ac","ca","ez","jb","jp","ms","pm","rj","rs")

data.brain <- read_cort_mag_data(fname.brain);
data.brain <- prepare_data(data.brain,amblyopes.initials,controls.initials)
data.brain <- attach_clinical_data(data.brain,fname.clinical)

# glimpse(data.brain)

data.behaviour <- read_behaviour_data(fname.behaviour)
data.behaviour <- prepare_data(data.behaviour,amblyopes.initials,controls.initials)
data.behaviour <- attach_clinical_data(data.behaviour,fname.clinical)

#### Figure 3 (radial error by eccentricity)----
if (doFig3 || doFigS2 || doStatsFig3 || doStatsFigS2){
  # scale/unscale radial and angular errors
  data.behaviour <- data.behaviour |>
    mutate(error_rnorm = error_r/eccentricity) |> # radial error, ecc-scaled
    mutate(error_thetadva = error_theta/180*pi * eccentricity) |>  # theta error in deg visual angle
    mutate(error_thetarad = data.behaviour$error_theta/180*pi) # theta error in radians (naturally ecc-scaled)

  ## calculate dependent measures: mean and standard deviation of each error type for each subject at each eccentricity
  errors <- data.behaviour %>%
    group_by(subject, group, eye, eyes, eccentricity) %>%
    summarise(radial_spread = sd(error_r),
              radial_bias = mean(error_r),
              radial_spread_scaled = sd(error_rnorm),
              radial_bias_scaled = mean(error_rnorm),
              angular_spread = sd(error_theta),
              angular_bias = mean(error_theta),
              angular_spread_dva = sd(error_thetadva),
              angular_bias_dva = mean(error_thetadva),
              .groups = "keep")
  errors<-as.data.frame(errors) # because zh prefers dataframes

  # Y axis limits for each panel
  sd_lims = c(0,1)
  sd_norm_lims = c(0,0.4)
  mean_lims = c(-1,0.5)
  mean_norm_lims = c(-0.3,0.3)

  # Function for plot each panel, will be applied to each of the 4 error measures
  plot.behaviour.fig3.panel <- function(error.measure, error.label, ylims){
    panel <- errors |>
      ggplot(aes(x = eccentricity,
                 y = get(error.measure), color = eyes)) +
      geom_jitter(size = pointsSizeBehaviour, alpha = pointsAlpha,
                  width = pointsJitterBehaviour) +
      geom_smooth(method = "lm",
                  formula = y ~ x,  # Corrected formula
                  se = TRUE, linewidth = lineWidthBehaviour, aes(fill=eyes)) +
      theColorScale +
      theFillScale +
      scale_x_continuous(limits = c(0.5, 7.5), breaks = c(1, 3, 5, 7)) +
      theme(aspect.ratio = 1,
            strip.background = NULL,
            legend.key.width = unit(1, "line")) +
      labs(x = 'Target eccentricity (º)',
           y = error.label) +
      scale_y_continuous(limits = ylims)
    panel
  }

  if (doFig3){
    # Panel A
    r_sd_plot <- plot.behaviour.fig3.panel("radial_spread",
                                           "Radial error SD (deg)",
                                           sd_lims)

    # Panel B
    r_mean_plot <- plot.behaviour.fig3.panel("radial_bias",
                                             "Mean radial error (deg)",
                                             mean_lims) +
      geom_hline(yintercept=0, linetype = "dotted", linewidth = 0.2)

    # Panel C
    rnorm_sd_plot <-  plot.behaviour.fig3.panel("radial_spread_scaled",
                                                "Radial error SD scaled",
                                                sd_norm_lims)

    # Panel D
    rnorm_mean_plot <-  plot.behaviour.fig3.panel("radial_bias_scaled",
                                                  "Mean radial error scaled",
                                                  mean_norm_lims) +
      geom_hline(yintercept=0, linetype = "dotted", linewidth = 0.2)

    # Combine into a single figure without legend
    FIG_behaviour <- plot_grid(
      r_sd_plot + theme(legend.position="none"),
      r_mean_plot + theme(legend.position="none"),
      rnorm_sd_plot + theme(legend.position="none"),
      rnorm_mean_plot + theme(legend.position="none"),
      align = 'vh',
      labels = c("A", "B", "C", "D"),
      nrow = 2
    )
    # Add legend back
    # extract the legend from one of the plots
    legend <- get_legend(
      # create some space to the left of the legend
      r_sd_plot #+ theme(legend.box.margin = margin(0, 0, 0, 12))
    )
    # add the legend to the row we made earlier.
    FIG_behaviour <- plot_grid(FIG_behaviour, legend, rel_widths = c(2, .4))


    if (saveFiguresAndTables==1){
      FIG_behaviour
      ggsave(plot = FIG_behaviour, file = "figure-3-sd_mean_radial_error.png",
             type = "cairo-png",  bg = "white",
             width = 13, height = 10.5, units = "cm", dpi = 800)
    }else{
      print(FIG_behaviour)
    }
  }

 #### Figure S2 (include angular error)----
  if (doFigS2){
    # Panel A
    rnorm_sd_plot <-  plot.behaviour.fig3.panel("angular_spread","Angular error SD (º)",
                                                c(0,17))

    # Panel B
    rnorm_mean_plot <-  plot.behaviour.fig3.panel("angular_bias",
                                                  "Mean angular error (º)",
                                                  c(-5,8)) +
      geom_hline(yintercept=0, linetype = "dotted", linewidth = 0.2)

    # Panel C
    r_sd_plot <- plot.behaviour.fig3.panel("angular_spread_dva",
                                           "Angular error SD (deg)",
                                           sd_lims)

    # Panel D
    r_mean_plot <- plot.behaviour.fig3.panel("angular_bias_dva",
                                             "Mean angular error (deg)",
                                             c(-0.5,0.5)) +
      geom_hline(yintercept=0, linetype = "dotted", linewidth = 0.2)

    # Combine into a single figure without legend
    FIG_behaviour <- plot_grid(
      rnorm_sd_plot + theme(legend.position="none"),
      rnorm_mean_plot + theme(legend.position="none"),
      r_sd_plot + theme(legend.position="none"),
      r_mean_plot + theme(legend.position="none"),
      align = 'vh',
      labels = c("A", "B", "C", "D"),
      nrow = 2
    )
    # Add legend back
    # extract the legend from one of the plots
    legend <- get_legend(
      # create some space to the left of the legend
      r_sd_plot #+ theme(legend.box.margin = margin(0, 0, 0, 12))
    )
    # add the legend to the row we made earlier.
    FIG_behaviour <- plot_grid(FIG_behaviour, legend, rel_widths = c(2, .4))


    if (saveFiguresAndTables==1){
      FIG_behaviour
      # 1.5 column figure: 4.86 inches by 9.19 inches
      ggsave(plot = FIG_behaviour, file = "figure-S2-sd_mean_angular_error.png",
             type = "cairo-png",  bg = "white",
             width = 13, height = 11, units = "cm", dpi = 800)
    }else{
      print(FIG_behaviour)
    }
  }

  #### Figure 3 models & tables----
  if(doStatsFig3==1 || doStatsFigS2==1 ){


    vars <- NULL
    var.names <- NULL
    if (doFig3){
      vars<-c(vars,"radial_spread",
              "radial_bias",
              "radial_spread_scaled",
              "radial_bias_scaled")
      var.names<-c(var.names,"Radial error SD (spread)",
                   "Mean radial error (bias)",
                   "Radial error SD scaled",
                   "Mean radial error scaled")
    }
    if (doFigS2){
      vars<-c(vars,"angular_spread_dva",
              "angular_bias_dva",
              "angular_spread",
              "angular_bias")
      var.names<-c(var.names,"Angular error spread (degrees of visual angle)",
                   "Mean angular error (degrees of visual angle)",
                   "Angular error SD (spread)",
                   "Mean angular error (bias)")
    }
    groups<-c("Controls", "Amblyopes")
    amb.eyes<-c("AE", "FE")

  ##### 1. compare eyes within groups----
    # cycles through dependent measures
    ## no difference between the left and right control eyes for any measure
    ## significant difference between amblyopic and fellow eyes for all measures except thetadegva (mean and sd)
    for(i in 1:length(vars)){
      tables.behaviour = list()
      for(j in 1:length(groups)){
        print("")
        print(var.names[i])
        print(groups[j])
        print("")

        lm.eye.within <- errors |> filter(group==substring(groups[j], 1, 1)) |>
          reduced.maximal.feasible.lm(as.formula(paste(
            vars[i], "~ eccentricity*eye+( 1 + eccentricity*eye|subject)") ) )
        print("")
        print(lm.eye.within[[1]])
        if (!is.null(lm.eye.within[[2]])){
          print(anova(lm.eye.within[[2]]))
          if (saveFiguresAndTables==1){
            #} && (has_fixed_effects(lm.eye.within[[2]])
            # || has_random_effects(lm.eye.within[[2]])) ){
            # need to re-fit model to calculate eta_squared
            if (has_random_effects(lm.eye.within[[1]])){
              lm.eye.within[[2]] <- lmer(lm.eye.within[[1]],
                                         data = errors[errors$group==
                                                         substring(groups[j], 1, 1),])
            }else if (has_fixed_effects(lm.eye.within[[1]])){
              lm.eye.within[[2]] <- lm(lm.eye.within[[1]],
                                       data = errors[errors$group==
                                                       substring(groups[j], 1, 1),])
            }
            tables.behaviour[[j]] <- publication.ready.anova.table(
              lm.eye.within[[2]],groups[j])
          }
        }

      }

      ## 2. compare AE, FE against control ----
      ## amblyopic vs. control: significant main effect of eyes & eyes*ecc interaction
      ## for all measures except:
      ## - mean(theta degva): nothing
      ## - sd(theta degva): only eyes*ecc interaction
      ## fellow vs. control: effect sig more often for radial than theta, mean than sd
      ## bias: sig diff for radial error only (main or interaction); sd: main effect
      ## only
      for(k in 1:2){
        print("")
        print(var.names[i])
        print(paste(amb.eyes[k], "vs. controls"))
        print("")
        lm.eye.between <- errors |> filter(eyes %in% c("control", amb.eyes[k])) |>
          reduced.maximal.feasible.lm(as.formula(paste(
            vars[i], "~eccentricity*group + (1 + eccentricity|subject)") ) )
        print("")
        print(lm.eye.between[[1]])
        if (!is.null(lm.eye.between[[2]])){
          print(anova(lm.eye.between[[2]]))
          if (saveFiguresAndTables==1){
            #} && (has_fixed_effects(lm.eye.between[[2]]) ||
            # has_random_effects(lm.eye.within[[2]])) ){
            # need to re-fit model to calculate eta_squared
            if (has_random_effects(lm.eye.between[[1]])){
              lm.eye.between[[2]] <- lmer(lm.eye.between[[1]],
                                          data = errors[errors$eyes %in%
                                                          c("control", amb.eyes[k]),])
            }else if (has_fixed_effects(lm.eye.within[[1]])){
              lm.eye.between[[2]] <- lm(lm.eye.between[[1]],
                                        data = errors[errors$eyes %in%
                                                        c("control", amb.eyes[k]),])
            }
            tables.behaviour[[length(groups)+k]] <- publication.ready.anova.table(
              lm.eye.between[[2]],paste(amb.eyes[k], "vs. controls"))
          }
        }

      }

      if (saveFiguresAndTables==1){
        # Combine tables and save
        combine.save.table(rbind(tables.behaviour[[1]],
                                 tables.behaviour[[2]],
                                 tables.behaviour[[3]],
                                 tables.behaviour[[4]]),
                           var.names[i], sprintf("Effects_%s.docx",vars[i]),
                           cumsum(c(nrow(tables.behaviour[[1]]),
                                    nrow(tables.behaviour[[2]]),
                                    nrow(tables.behaviour[[3]]))))
      }
    }
  }

}

#### Figure 5 (pRF measures)----

if(doFig5==1){

  the_mapping_ecc <- aes(x = x_coords, y = pRF_ecc,
                     color=eyes,
                     linetype=hemisphere,
                     group=interaction(eyes, hemisphere))  # eyes in group

  ECC_plot <- data.brain |>
    ggplot(the_mapping_ecc) +
    geom_jitter(size = pointsSizePRF, alpha = pointsAlpha, width = pointsJitterPRF) +
    geom_smooth(method="nls",
                formula = y ~ bensonEcc(x,Q),
                method.args = list(start=c(Q=3)), se=FALSE,
                mapping=the_mapping_ecc,
                linewidth=lineWidthPRF) +
    geom_smooth(method="loess", formula = y ~ x,
                mapping=the_mapping_ecc,
                se=FALSE, linewidth=lineWidthPRF/2, alpha=loessAlpha)+
    theColorScale +
    theLineTypeScale +
    scale_x_continuous(limits=c(0.0,0.51), breaks = seq(0, 1, by=0.2)) +
    theme(aspect.ratio = 1,
          strip.background = NULL,
          strip.text = element_blank(),
          legend.key.width = unit(1, "line")) +
    labs(x = 'Normalised cortical coordinates',
         y = "Eccentricity (º)") +
    facet_wrap(~group)

  ECC_plot

  # what does r2 look like? remembering this is the (nan)mean across the same x coords
  data.brain |>
    ggplot(aes(x = x_coords,
               y = r2data, color=eyes, linetype = hemisphere,
               group=interaction(eyes, hemisphere))) +
    geom_point(size = pointsSizePRF, alpha = pointsAlpha) +
    geom_smooth(method = "loess", formula = y ~ x, se=FALSE) +
    theColorScale +
    labs(title="Mean r2 value vs cortical distance")


  # pRF size plot
  the_mapping_rf <- aes(x = x_coords, y = pRF_size,
                        color=eyes, fill=eyes, linetype=hemisphere,
                        group=interaction(eyes, hemisphere))

  RF_plot <-  data.brain |>
    ggplot(the_mapping_rf) +
    geom_jitter(size = pointsSizePRF, alpha = pointsAlpha, width = pointsJitterPRF) +
    geom_smooth(method = "lm",
                formula = y ~ x, mapping=the_mapping_rf, linewidth=lineWidthPRF) +
    theColorScale +
    theFillScale +
    theLineTypeScale +
    scale_x_continuous(limits=c(0,0.5), breaks = seq(0, 0.6, by=0.2)) +
    scale_y_continuous(limits = c(0,2.2)) +
    theme(aspect.ratio = 1,
          strip.background = NULL,
          strip.text = element_blank(),
          legend.key.width = unit(1, "line")) +
    labs(x = 'Normalised cortical coordinates',
         y = "pRF size, sigma (º)") +
    facet_wrap(~group)

  RF_plot

  # Combine into a single figure without legend
  FIG_cortex <- plot_grid(
    ECC_plot + theme(legend.position="none"),
    RF_plot + theme(legend.position="none"),
    align = 'vh',
    vjust = 1.1,
    labels = c("A", "B"),
    nrow = 2
  )
  # Add legend back
  # extract the legend from one of the plots
  legend <- get_legend(
            # create some space to the left of the legend
    RF_plot #+ theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  # add the legend to the row we made earlier.
  FIG_cortex <- plot_grid(FIG_cortex, legend, rel_widths = c(1, .23))

  if (saveFiguresAndTables==1){
    FIG_cortex
    # 1.5 column figure: 4.86 inches by 9.19 inches
    ggsave(plot = FIG_cortex, file = "figure-5-pRF_ecc_width.png",
           type = "cairo-png",  bg = "white",
           width = 11, height = 9, units = "cm", dpi = 800)
  }else{
    print(FIG_cortex)
  }
}

