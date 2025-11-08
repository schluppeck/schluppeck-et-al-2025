# helper code
#
# in conjunction with 2025_make_amblyopia_summary


#' function for reading in data from the pRF CSV file and do the following:
#'    - select reasonable cortical coordinates
#'    - create log(e) variable
#' 
#' @param fname.brain The filename of the pRF csv file
#' @param minCoord A cutoff value for reasonable min x (cortical coordinates)
#'                 The default is 0 (ie include all data, even around areas corresponding
#'                 to fixation)
#' @param maxCoord A cutoff value for reasonable max x (cortical coordinates)
#'                 The default is 0.5 - beyond this the stimuli in the fMRI experiment
#'                 led to biased pRF fits (as seen by eccentricity values that drop again)
#'
read_cort_mag_data <- function(fname.brain, minCoord = 0, maxCoord = 0.5){
  message(      '(read_cort_mag_data) read pRF data and restricted x - cortical distance - ')
  message(paste('                             to between', minCoord, 'and', maxCoord))
  
  data <- read_csv(fname.brain, show_col_types = FALSE)|>
    clean_names() |>
    arrange(which_eye) |> # make sure data are sorted according to whihc_eye before converting to factors
    mutate_if(is_character, as_factor) |>
    select(-sub)
  
  
  # filter based on min and max coordinates and add log(e) variable
  data <- data |>
    filter(e_coords >= minCoord & e_coords <= maxCoord)|>
    mutate(log_pRF_ecc = log(e_data)) |>
    filter(is.finite(log_pRF_ecc))
  
  # change names of pRF data
  names(data)[names(data) == "rf_data"] <- "pRF_size"
  names(data)[names(data) == "e_data"] <- "pRF_ecc"
  names(data)[names(data) == "e_coords"] <- "x_coords"
  names(data)[names(data) == "hemi"] <- "hemisphere"
  names(data)[names(data) == "which_eye"] <- "eye"
  names(data)[names(data) == "kind"] <- "group"
  
  data # return the dataset
}

#' function for reading in data from the pRF CSV files and do the following:
#'    - rename some variable names and levels names so they match the brain data
#' 
#' @param fname.behaviour The filename of the pRF csv file
#' 
read_behaviour_data <- function(fname.behaviour){

  message(      '(read_behaviour_data) read behaviour data')

  data <- read_csv(fname.behaviour, show_col_types = FALSE)|>
    clean_names() |>
    arrange(which_eye) |> # make sure data are sorted according to whihc_eye before converting to factors
    mutate_if(is_character, as_factor) |>
    select(-sub_num)
  
  # change variable names
  names(data)[names(data) == "sub"] <- "sub_name" # so this matches the brain data
  names(data)[names(data) == "which_eye"] <- "eye"
  names(data)[names(data) == "kind"] <- "group"
  names(data)[names(data) == "target_radius_deg"] <- "eccentricity"
  
  # change level names of eye
  levels(data$eye) <- c("AE","FE")

  data # return the dataset
}


#' function to prepare behavioural or brain data by:
#'    - adding numbered subject names (subject)
#'    - remove participants who are not in the predefined list of amblyopes or controls
#'    - correcting stimulated eye information for controls
#'    - creating new variable eyes with three levels (AE, FE and control)
#' 
#' @param data The brain or behaviour dataframe
#' @param amblyopes.initials The list of amblyopes to keep
#' @param controls.initials The list of controls to keep
#' 
prepare_data <- function(data,amblyopes.initials,controls.initials){
  
  # create numbered subject names
  data <- data |>
    mutate(subject = ifelse(group == "A", match(sub_name,amblyopes.initials),match(sub_name,controls.initials)))
  
  participants.to.exclude <- as.character(unique(data$sub_name[is.na(data$subject)]))
  if (length(participants.to.exclude) > 0){
    message(      sprintf('(prepare_data) Excluded participants %s', paste(participants.to.exclude, collapse = " ")))
    
    data<- data |>
      filter(!is.na(subject)) |> # delete any participant without subject (not in either list)
      mutate(subject = as_factor(sprintf('%s%02d', group, subject)))
  }
  
  # create new variable "eyes" with three levels (AE, FE and control)
  data$eyes <- data$eye
  levels(data$eyes) <- c(levels(data$eyes),"control")
  data$eyes[data$group == "C"] <- "control"
  
  # correct eye labels for controls (FE is always left eye and AE is always right eye)
  # (instead this should be corrected directly in the file)
  actual.eyes <- c("L","R")
  incorrect.eyes <- levels(data$eye)
  levels(data$eye) <- c(levels(data$eye),actual.eyes)
  levels(incorrect.eyes) <- levels(data$eye)
  for (i in 1:length(actual.eyes)){
    data$eye[data$eye==incorrect.eyes[i] & data$group=="C"]<-actual.eyes[i]
  }
  
  # re-order variables
  data <- data |>
    select(subject, eye, group, eyes, everything())
  
  data # return the dataset
}


#' function for adding clinical data to dataframe (for amblyopes)
#' 
#' @param fname.brain The filename of the pRF csv file
#' @param fname.clinical The filename of the associated clinical csv file
#' 
attach_clinical_data <- function(data, fname.clinical){
  
  # join in clinical data, including which eye (L/R) was affected
  eye_data <- read_csv(fname.clinical, show_col_types = FALSE) |>
    clean_names() |>
    mutate(sub_name = as_factor(str_to_lower(subject)),
           affected_eye = eye, .keep="unused") |>
    select(sub_name, everything())
  
  # now left_join (keeps all brain or behaviour data, but does not include clinical data from irrelevant participants)
  data <- left_join(data, eye_data, by="sub_name")
  
  data # return the dataset
}

#' function for finding the maximal feasible model and reduce it by backward elimination of non-significant fixed-effects
#' 
#' @param data dataframe
#' @param formula linear model formula for full model with all fixed and random effects
#' 
reduced.maximal.feasible.lm <- function(data,formula){
  
  print("(reduced.maximal.feasible.lm) Find maximal feasible model using buildmer...") # (i.e. including the largest number of random effects and still converging)
  m <- buildmer(formula,data=data,
                buildmerControl=buildmerControl(direction='order',
                                                quiet = TRUE,
                                                args=list(control=lmerControl(optimizer='bobyqa'))))
  f.max <- formula(m@model) # formula of maximal feasible model
  
  lm.reduced = NULL
  if (has_random_effects(f.max) | has_fixed_effects(f.max)){ # if there are some effects left
    print("(reduced.maximal.feasible.lm) Backward elimination of non-significant effects...")
    m <- buildmer(f.max,data=data,buildmerControl=list(direction='backward',
                                                       quiet = TRUE,
                                                       args=list(control=lmerControl(optimizer='bobyqa'))))
    f.reduced <- formula(m@model)  # reduced model
    if (has_random_effects(f.reduced)){
      lm.reduced <-lmer(f.reduced,data=data) # fit and test fixed effects
    }else if (has_fixed_effects(f.max)){
      lm.reduced <-lm(f.reduced,data=data) # fit and test fixed effects
      print("(reduced.maximal.feasible.lm) No random effects left after backward elimination, using lm for fit.")
    }else{
      f.reduced = NULL
      print("(reduced.maximal.feasible.lm) No random or fixed effects left after backward elimination.")
    }
    
  }else{
    f.reduced = f.max
    print("(reduced.maximal.feasible.lm) No random or fixed effects left, aborting.")
  }
  
  return(list(f.reduced,lm.reduced))
}

#' function testing whether a linear model formula contains fixed effects
#' 
#' @param data formula
#' 
has_fixed_effects <- function(formula) {
  # Extract terms from formula
  terms_info <- terms(formula)
  # Count fixed effects (excluding intercept)
  length(attr(terms_info, "term.labels")) > 0
}

#' function testing whether a linear model formula contains random effects
#' 
#' @param data formula
#' 
has_random_effects <- function(formula) {
  length(findbars(formula))>0
}

#' function formatting an anova table and adding eta-squared
#' 
#' @param model dataframe
#' @param model.name the name of the model (will be pre-prended to the formula)
#' 
publication.ready.anova.table <- function(model,model.name){

  # Get Type III ANOVA tables
  anova_results <- anova(model)
  n.effects = ifelse(class(model) == "lm",nrow(anova_results)-1,nrow(anova_results)) # number of effects
  
  p.thresholds <- c(0, .001, .01, 0.05, 1)
  labs <- c("***","**", "*", NA)
  
  # Compute effect sizes
  if (n.effects > 0){
    eta2 <- eta_squared(model, partial = TRUE)
    eta2.partial <- coalesce(eta2$Eta2_partial, eta2$Eta2) # take Eta2_partial if it exists, Eta2 is not
    cohen_f <- sqrt(eta2.partial / (1 - eta2.partial))
    
    Effects <- rownames(head(anova_results,n=n.effects))
    SS <- head(anova_results,n=n.effects)$"Sum Sq"
    if (class(model) == "lm"){
      NumDF <- head(anova_results,n=n.effects)$Df
      DenDF <- tail(anova_results$Df,n=1)
    }else{
      NumDF <- anova_results$NumDF
      DenDF <- anova_results$DenDF
    }
    F <- head(anova_results,n=n.effects)$"F value"
    p <- head(anova_results,n=n.effects)$"Pr(>F)"
    signif. <- labs[findInterval(p, p.thresholds)]
    eta.2 <- head(eta2.partial,n=n.effects)
    Cohens_f <- head(cohen_f,n=n.effects)
    
  }else{
    Effects <- "None"
    SS <- NA
    NumDF <- NA
    DenDF <- NA
    F <- NA
    p <- NA
    signif. <- NA
    eta.2 <- NA
    Cohens_f <- NA
  }

  # Convert to dataframe
  table_data <- data.frame(Model = sub(" \\+ \\(1.*", "", paste(model.name,deparse1(formula(model)),sep=": ")), # remove random effects from formula and prepend model name
                           Effects,SS,NumDF,DenDF,F,p,signif.,Cohens_f) |>
    mutate(across(where(is.numeric)& !any_of(c("p","NumDF")), ~ formatC(.x, format = "f", digits = 2))) # keep 2 decimals for all numeric variables except p
  if (n.effects > 0){
    table_data <- table_data |> mutate(p = formatC(p, format = "g", digits = 2)) # keep 2 significant digits for p
  }
  table_data

}


#' function saving a table to a Word document
#' 
#' @param input.table formatted table
#' @param table.title table title
#' @param file.name output file name
#' @param line.rows rows after which to add an horizontal line
#' 
combine.save.table <- function(input.table,table.title,file.name,line.rows){

  anova_flextable <- flextable(input.table) |>
    set_caption(table.title) |>
    merge_v(j = 1) |> # Merge model names vertically
    hline(i = line.rows, border = fp_border(width = 1.5, color = "gray")) |>  # Lines between each model
    set_header_labels(Model = "Reduced model", Cohens_f = "Cohen's f") |>
    align(align = "center", part = "header") |> # center-align all headers
    align(j = 2:9, align = "center", part = "body") |> # center-align columns 2:9
    width(j = 1, width = 1.8) |> # Adjust width of the model column narrow to force wrapping
    width(j = 2, width = max(0.5,max(nchar(input.table$Effects))/12) ) |> # Adjust width of the effects column according to max string length
    width(j = c(3,7,9), width = 0.8) |> # Adjust the SS, p value and cohen's f column wide to avoid wrapping
    set_table_properties(layout = "fixed") # Ensure word wrap works (do not use "autofit")
  
  ### Save the Merged Table to a Word Document
  doc <- read_docx() |>
    body_add_flextable(anova_flextable) |>
    body_add_par("") |>  # Add space
    body_end_section_landscape()  # Set landscape orientation
  
  # Save the document
  print(doc, target = file.name)

}

#' simple exponential fit to eccentricity data
#' 
#' takes in a cortical coordinate on [0,1] and returns eccentricity value in degrees 
#' 
#' @param x cortical distance, normalised units
#' @param Q exponent
#' 
bensonEcc <- function(x, Q) {
  90*exp(Q * (x-1))
}

# E = 90*exp(Q*(x-1))
# log(e) = log(90) + Q*(x-1)
# 1/Q(log(e) - log(90)) + 1 = x

#' inverse of that function, to allow for some additional calculations
#' 
#' takes in an Eccentricity value in degrees and returns cortical coordinate on [0,1]
#' 
#' @param E cortical distance, normalised units
#' @param Q exponent
#' 
invBenson <- function(E, Q) {
  # return E in cortical coords
  (log(E) - log(90))/Q + 1
}
