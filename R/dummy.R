# This script contains some functions to automate things:

# cat("Last updated: 23/03/2021 - 15:14")

# Function to automize dummy encoding:
RIVMgrowth_dummy <- function(
  data, # Dataset
  # design, # Design matrix, as in psychonetrics
  type = c("non-linear","linear"), # Type of analysis to do
  independent, # Variable name of predictor (only one supported at the moment)
  time, # time vector encoding distances in time
  ordered = FALSE,
  ...
){
  # library("dplyr")
  
  # Check data:
  stopifnot(is.data.frame(data))
  
  
  # Form design matrix:
  dots <- list(...)
  
  # FIXME: Still no support for multiple dependent variables...
  if (length(dots) > 1){
    stop("Currently only one dependent variable is supported.")
  }
  design <- t(dots[[1]])
  rownames(design) <- names(dots)
  
  # Check ordered:
  if (length(ordered) == 1){
    ordered <- rep(ordered,nrow(design))
  }
  
  # Check design matrix:
  allVars <- unlist(design)
  allVars <- allVars[!is.na(allVars)]
  if (!all(allVars %in% names(data))){
    stop("Not all names in the design matrix correspond to variables in the dataset.")
  }
  
  # Check ordered:
  if (length(ordered)!=nrow(design)){
    stop("'ordered' must be a logical vector with TRUE or FALSE for each row in the design matrix (dependent variable).")
  }
  for (i in seq_along(ordered)){
    if (ordered[i]){
      for (j in seq_along(ncol(design))){
        data[[design[i,j]]] <- as.ordered(data[[design[i,j]]])
      }
    } else {
      for (j in seq_along(ncol(design))){
        if (is.ordered(data[[design[i,j]]])||is.factor(data[[design[i,j]]])||!is.numeric(data[[design[i,j]]])){
          stop("Variable is not continuous.")
        }
      }
    }
  }
  
  # Extract from design:
  if (!is.matrix(design)){
    design <- t(design)
  }
  nVar <- nrow(design)
  nWave <- ncol(design)
  
  # Names of design:
  if (is.null(rownames(design))){
    rownames(design) <- paste0("Var_",seq_len(nVar))
  }
  variables <- rownames(design)
  
  # Time vector:
  if (missing(time)){
    time <- seq_len(nWave)
  }
  
  # Type of analysis:
  type <- match.arg(type)
  
  # Check if predictor is present:
  use_predictor <- !missing(independent)
  
  # Setup dummy:
  if (use_predictor){
    # Check data:
    if (!all(independent %in% names(data))){
      stop("independent does not correspond to a variable(s) in the dataset")
    }
    
    
    # Make factors if needed:
    for (v in independent){
      if (!is.factor(data[[v]])){
        data[[v]] <- as.factor(data[[v]])
      }      
    }
    
    # If there are multiple predictors, we need to collapse these in one:
    
    # first count all cases:
    predictor_summary <- data %>% group_by(across({{independent}})) %>% tally
    
    # Now make a dummy predictor:
    predictor_summary$PREDICTOR <- as.factor(seq_len(nrow(predictor_summary)))
    
    # And left join this to the dataset:
    data <- data %>% left_join(
      predictor_summary %>% select(-n),
      by = independent
    )
    
    # Now we can just continue as if we have one predictor:
    original_predictors <- independent
    independent <- "PREDICTOR"

    # Make dummy variables:
    predictor_levels <- levels(data[[independent]])
    nlevels_predictor <- length(predictor_levels)
    nDummy <- nlevels_predictor-1
    
    # Dummy coding:
    dummyvars <- paste0("DUMMY_",seq_len(nDummy))
    for (d in seq_len(nDummy)){
      data[[dummyvars[[d]]]] <- 1* (data[[independent]] == predictor_levels[d+1])
    }
    
    # Also compute marginal dummy variables for marginal distributions:
    marginal_dummies <- lapply(seq_along(original_predictors),function(i){
      data %>% group_by_(original_predictors[i]) %>% summarize(
        across(dummyvars, ~ mean(.x, na.rm=TRUE)), .groups = "drop_last"
      ) %>% ungroup
    }) 
    names(marginal_dummies) <- original_predictors
    
  } else {
    nDummy <- 0
    dummyvars <- character(0)
  }

  # Names of parameters needed in the model:
  intercepts <- paste0("int_",variables)
  slopes <- paste0("slope_",variables)
  residual <- paste0("resid_",variables)
  
  # Slope parameters:
  if (type == "non-linear"){
    slope_pars <- lapply(variables,function(x){
      c(0,1,paste0("slope_",x,"_",3:nWave))
    })  
    names(slope_pars) <- variables
  } else {
    slope_pars <- lapply(variables,function(x){
      time - time[1]
    })
    names(slope_pars) <- variables
  }
  
  # Intercept mean parameters:
  int_mean_pars <- lapply(seq_along(variables),function(i){
    res <- paste0("intercept_model_",variables[i],"_",seq_len(nDummy+1)-1)
    if (ordered[i]){
      res[1] <- 0
    }
    return(res)
  })
  names(int_mean_pars) <- variables
  
  # Slope mean parameters:
  slope_mean_pars <- lapply(variables,function(x){
    paste0("slope_model_",x,"_",seq_len(nDummy+1)-1)
  })
  names(slope_mean_pars) <- variables
  
  # Threshold parameters:
  thresholds <- lapply(seq_along(variables),function(i){
    if (!ordered[i]){
      return(NA)
    } else {
      # How many levels?
      nLevel <- length(na.omit(unique(unlist(data[,design[i,]]))))
      paste0("t_",i,"_",seq_len(nLevel-1))
    }
  })
  
  # Start forming the model:
  mod <- character(0)
  cur_line <- 0
  
  # loop over variables:
  for (v in seq_len(nVar)){
    # Variables from the design matrix:
    curvars <- design[v,]
    
    # Which are not NA?
    incl <- !is.na(curvars)
    inclvars <- which(incl)
    
    # Cut out NA:
    curvars <- curvars[!is.na(curvars)]
    
    # Update the line:
    cur_line <- cur_line + 2
    
    # Header:
    mod[cur_line] <- paste0("# Intercepts for variable ",variables[v])
    
    # Update the line:
    cur_line <- cur_line + 1
    
    # Write intercepts:
    mod[cur_line] <- paste0(intercepts[v], " =~ ", paste0("1 * ",curvars, collapse = " + "))
    
    # Update the line:
    cur_line <- cur_line + 1
    
    # Now model the intercept:
    mod[cur_line] <- paste0(intercepts[v], " ~ ", paste0(int_mean_pars[[v]],"*",c("1", dummyvars), collapse = " + "))
    
    # Update the line:
    cur_line <- cur_line + 2
    
    # Header:
    mod[cur_line] <- paste0("# Slopes for variable ",variables[v])
    
    # Update the line:
    cur_line <- cur_line + 1
    
    # Write slopes:
    mod[cur_line] <- paste0(slopes[v], " =~ ", paste0(slope_pars[[v]][incl]," * ",curvars, collapse = " + "))
    
    # Update the line:
    cur_line <- cur_line + 1
    
    # Now model the slopes:
    mod[cur_line] <- paste0(slopes[v], " ~ ", paste0(slope_mean_pars[[v]],"*",c("1", dummyvars), collapse = " + "))
    
    # Update the line:
    cur_line <- cur_line + 2
    
    # Header:
    mod[cur_line] <- paste0("# Residuals for variable ",variables[v])
    
 
    
    # If the variable is ordered, constrain the thresholds:
    if (ordered[v]){
      # Equal thresholds for every variable:
      for (i in seq_along(curvars)){
        # Update the line:
        cur_line <- cur_line + 1
        
        # Write the thresholds: 
        mod[cur_line] <- paste0(curvars[i], " | ", paste0(thresholds[[v]]," * t",seq_along(thresholds[[v]]),collapse=" + "))
      }
      
    } else {
      # Else write residuals:
      for (i in seq_along(curvars)){
        # Update the line:
        cur_line <- cur_line + 1
        
        mod[cur_line] <- paste0(curvars[i], " ~~ ", residual[v]," * ",curvars[i])
      }
    }
    
  }
  
  
  # Now add all expected values of all groups at all waves for all variables:
  if (use_predictor){
    ev_df <- expand.grid(
      predictor = seq_along(predictor_levels),
      wave = seq_len(nWave),
      var = variables,
      stringsAsFactors = FALSE
    )
    names(ev_df) <- c("predictor","wave","var")
    
    # Add dummies:
    for (d in seq_len(nDummy)){
      ev_df[[dummyvars[[d]]]] <- 1* (predictor_levels[ev_df$predictor] == predictor_levels[d+1])
    }
    
    # Parameter name:
    ev_df$par <-  paste0(ev_df$var,"_wave",ev_df$wave,"_",ev_df$predictor)
    
    # Make PREDICTOR column to match:
    ev_df$PREDICTOR <- predictor_levels[ev_df$predictor]
    
    # Left join the predictor matrix:
    ev_df <- ev_df %>% left_join(
      predictor_summary %>% select(-n),
      by = "PREDICTOR"
    ) %>% select(
      -predictor, -PREDICTOR
    )
    
    # Reorder cols:
    nPred <- ncol(predictor_summary) - 2
    ev_df <- ev_df[,c((ncol(ev_df)-nPred+1):(ncol(ev_df)),1:(ncol(ev_df)-nPred))]
    
    # Also make a data frame for the marginal expected values:
    for (i in seq_along(marginal_dummies)){
      marginal_dummies[[i]]$predictor_var <- original_predictors[i]
      marginal_dummies[[i]]$predictor_level <- seq_len(nrow(marginal_dummies[[i]]))
      marginal_dummies[[i]]$predictor_label <- marginal_dummies[[i]][[1]]
    }
    ev_df_marginal_sub <- bind_rows(marginal_dummies)
    ev_df_marginal_sub <- ev_df_marginal_sub[,c("predictor_var","predictor_level","predictor_label", dummyvars)]
    ev_df_marginal_sub$ID <- seq_len(nrow(ev_df_marginal_sub))
    
    ev_df_marginal <- expand.grid(
      ID = unique(ev_df_marginal_sub$ID),
      wave = seq_len(nWave),
      var = variables,
      stringsAsFactors = FALSE
    ) %>% left_join(ev_df_marginal_sub, by = c("ID")) %>% select(-ID)
    
    # Parameter name:
    ev_df_marginal$par <-  paste0(ev_df_marginal$var,"_wave",ev_df_marginal$wave,"_",ev_df_marginal$predictor_var,"_",ev_df_marginal$predictor_level)
    
  } else {
    ev_df <- expand.grid(
      wave = seq_len(nWave),
      var = variables,
      stringsAsFactors = FALSE
    )
    
    # Parameter name:
    ev_df$par <-  paste0(ev_df$var,"_wave",ev_df$wave)
    
    nPred <- 0
    
    ev_df_marginal <- NULL
  }

  # Update the line:
  cur_line <- cur_line + 2
  
  # Header:
  mod[cur_line] <- "# Dummy encoding"
  

  # Loop over these to add:
  for (i in seq_len(nrow(ev_df))){
    var <- ev_df$var[i]
    wave <- ev_df$wave[i]
    dummies <- unlist(ev_df[i,dummyvars])
    
    
    # Current parameter:
    curpar <- ev_df$par[i]
    
    # Update the line:
    cur_line <- cur_line + 1

    # Write the line:
    mod[cur_line] <-  paste0(
      curpar, " := ", paste0(
        int_mean_pars[[var]], "*", c("1",dummies),
        collapse = " + "
      ), " + ",
      paste0(
        slope_mean_pars[[var]], "*", c("1",dummies), " * ", slope_pars[[v]][wave],
        collapse = " + "
      )
    )
  }
  
  # Also use dummies for marginal predictions:
  if (!is.null(ev_df_marginal)){
    for (i in seq_len(nrow(ev_df_marginal))){
      var <- ev_df_marginal$var[i]
      wave <- ev_df_marginal$wave[i]
      dummies <- unlist(ev_df_marginal[i,dummyvars])
      
      
      # Current parameter:
      curpar <- ev_df_marginal$par[i]
      
      # Update the line:
      cur_line <- cur_line + 1
      
      # Write the line:
      mod[cur_line] <-  paste0(
        curpar, " := ", paste0(
          int_mean_pars[[var]], "*", c("1",dummies),
          collapse = " + "
        ), " + ",
        paste0(
          slope_mean_pars[[var]], "*", c("1",dummies), " * ", slope_pars[[v]][wave],
          collapse = " + "
        )
      )
    }
    
    # for (i in seq_along(marginal_dummies)){
    #   
    #   # Write header:
    #   cur_line <- cur_line + 3
    #   mod[cur_line] <- paste0("# Marginal Dummy encoding - ",original_predictors[i])
    #   
    #   # Loop over waves:
    #   for (j in seq_len(nrow(marginal_dummies[[i]]))){
    #     for (v in variables){
    #       for (w in seq_len(nWave)){
    #         
    #         var <- v
    #         wave <- w
    #         dummies <- unlist(marginal_dummies[[i]][j,dummyvars])
    #         
    #         # Current parameter:
    #         curpar <- paste0(original_predictors[i],"_",j,"_",var,"_wave",wave)
    #         
    #         # Update the line:
    #         cur_line <- cur_line + 1
    #         
    #         # Write the line:
    #         mod[cur_line] <-  paste0(
    #           curpar, " := ", paste0(
    #             int_mean_pars[[var]], "*", c("1",dummies),
    #             collapse = " + "
    #           ), " + ",
    #           paste0(
    #             slope_mean_pars[[var]], "*", c("1",dummies), " * ", slope_pars[[v]][wave],
    #             collapse = " + "
    #           )
    #         )
    #         
    #         
    #       }
    #     }
    #   }
    #  
    # 
    #   
    #   
    # }
    
    
  }
  
  # Fix empty lines:
  mod[is.na(mod)] <- ""
  
  # Aggregate:
  model <- paste0(mod, collapse = "\n")

  # No FIML for ordered categorical data:
  missing <- ifelse(any(ordered),"pairwise","fiml")
  
  # Run lavaan model:
  fit <- growth(model, data, missing = missing)

  # Connect the parameter estimates to the expected value table:
  parTable <- parameterEstimates(fit, standardized = TRUE)
  
  # Join with the expected values:
  ev_df <- ev_df %>% left_join(parTable, by = c("par" = "label"))
  
  # Add time:
  ev_df$time_of_wave <- time[ev_df$wave]
  
  # Also join marginal expected values:
  if (!is.null(ev_df_marginal)){
    # Join with the expected values:
    ev_df_marginal <- ev_df_marginal %>% left_join(parTable, by = c("par" = "label"))
    
    # Add time:
    ev_df_marginal$time_of_wave <- time[ev_df_marginal$wave]
  }
    
  # Fix predictor:
  # if (use_predictor){
  #   ev_df$predictor <- predictor_levels[ev_df$predictor]
  # }
  
  
  # Select columns:
  if (use_predictor){
    ev_df <- ev_df[,c( names(predictor_summary)[1:(ncol(predictor_summary)-2)],"var","wave","time_of_wave","est","se","ci.lower","ci.upper")]
    ev_df_marginal <- ev_df_marginal[,c("predictor_var","predictor_level","predictor_label","var","wave","time_of_wave","est","se","ci.lower","ci.upper")]
    
  } else {
    ev_df <- ev_df[,c("var","wave","time_of_wave","est","se","ci.lower","ci.upper")]
  }
  

  
  # Significance of the regressions:
  reg_df <- parTable %>% filter(op %in% c("~1", "~"), lhs %in% c(intercepts, slopes)) %>% 
    select(-label)
  
  # Correlations between slopes:
  cor_df <- parTable %>% filter(op == "~~", lhs %in% c(intercepts, slopes), rhs %in% c(intercepts, slopes)) %>% 
    select(lhs, op, rhs, est, se, z, pvalue, ci.lower, ci.upper, std.all)

  
  # Make list:
  results <- list(
    lavresult = fit,
    expected_values = ev_df,
    marginal_expected_values = ev_df_marginal,
    regressions = reg_df,
    correlations = cor_df,
    time = time,
    use_predictor = use_predictor,
    n_predictor = nPred,
    ordered = ordered,
    lav_model = model
  )
  class(results) <- "RIVMgrowth_dummy"
  
  
  # Return:
  return(results)
}

# Print method:
print.RIVMgrowth_dummy <- function(x){
  cat("=== RIVM Latent Growth Curve Analysis ===","\n\n")
  
  cat("Lavaan fit result: \n")
  print(x$lavresult)
  
  cat("\n\nLavaan fit measures: \n")
  fit <- round(fitMeasures(x$lavresult)[c("chisq", "df", "pvalue", "cfi", "tli", "nnfi", "rfi", 
                            "nfi", "pnfi", "ifi", "rni", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", 
                            "rmsea.pvalue", "srmr", "gfi", "agfi")], 2)
  
  df <- data.frame(
    index = names(fit),
    value = fit
  )
  rownames(df) <- NULL
  print(df, digits = 2)
  
  cat("\n\nRegression coefficients: \n")
  print(x$regressions, digits = 2)
  
  cat("\n\n(co)variance and correlation coefficients: \n")
  print(x$correlations, digits = 2)
  
  cat("\n\nExpected values: \n")
  print(x$expected_values, digits = 2)
  
  cat("\n\nMarginal Expected values: \n")
  print(x$marginal_expected_values, digits = 2)
  
  cat("\n\n=========================================")
}


# Plot method:
plot.RIVMgrowth_dummy <- function(x, dependent, independent, wave_label = "Wave "){
  # library("ggplot2")
  
  # Select variable:
  if (missing(dependent)){
    if (length(unique(x$expected_values$var))==1){
      dependent <- x$expected_values$var[1] 
    } else {
      stop(paste0("Dependent variable needs to be assigned using the 'dependent' argument. Set this to one of: ",paste0('"',unique(x$expected_values$var),'"',collapse=", ")))
    }
  }
  
  if (x$use_predictor){
    
    
    
    if (missing(independent)){
      predVars <- unique(x$marginal_expected_values$predictor_var)
      if (length(predVars) > 1){
        stop(paste0("Independent variable needs to be assigned using the 'independent' argument. Set this to one of: ",paste0('"',predVars,'"',collapse=", ")))
      } else {
        independent <- predVars[1]
      }
    }
    
    # Select data to use:
    useData <- x$marginal_expected_values %>% filter(
      predictor_var == independent, var == dependent
    )
    
    # How many levels:
    nLevel <- length(unique(useData$predictor_level))
    
    # Plot:
    p <- ggplot(useData, 
                
                # With aesthetics:
                aes(
                  x = time_of_wave, # This is the variable that enocdes the time of the wave.
                  y = est, # The estimate
                  ymin = ci.lower, # If we want CIs
                  ymax = ci.upper, # If we want CIs
                  colour = predictor_label,
                  fill = predictor_label
                )) +  
      
      # Nicer legend (and colorblind theme):
      scale_colour_manual(independent, 
                          values = qgraph:::colorblind(nLevel)) +
      scale_fill_manual(independent, 
                          values = qgraph:::colorblind(nLevel)) 
      
      
  } else {
    p <- ggplot(x$expected_values, 
                
                # With aesthetics:
                aes(
                  x = time_of_wave, # This is the variable that encodes the time of the wave.
                  y = est, # The estimate
                  ymin = ci.lower, # If we want CIs
                  ymax = ci.upper # If we want CIs
                ))
  }
  
  p <- p +
    
    # Seperate plot per variable:
    facet_grid(var ~ ., scales = "free_y") +
    
    # Add confidence region (uncomment to remove):
    geom_ribbon(alpha = 0.1, show.legend = FALSE, colour = NA) + 
    
    # Add lines:
    geom_line(lwd=1.2) + 
    
    # Add points:
    geom_point(cex = 3, pch = 21, fill = "white", stroke = 1.5) + 
    
    # Nicer theme:
    theme_bw() +
    
    # Set x axis scale to have waves instead of numbers:
    scale_x_continuous(breaks = sort(unique(x$expected_values$time_of_wave)),
                       labels = paste0(wave_label," ",sort(unique(x$expected_values$wave)))) + 
    
    # Nice labels:
    ylab(dependent) + 
    xlab("") +
    
    # Larger fonts:
    theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  
  return(p)
}


# Correlation plot:
corplot <- function(x, legend.position = "topright"){
  # library("dplyr")
  # library("qgraph")
  
  cors <- x$correlations
  cors <- cors %>% filter(lhs != rhs) %>% select(
    lhs = lhs, rhs = rhs, weight = std.all, pvalue = pvalue
  ) %>%
    mutate(
      lhs = gsub("\\_","\n",lhs),
      rhs = gsub("\\_","\n",rhs)
    ) %>% mutate(
      lhs = gsub("int\n","intercept\n",lhs),
      rhs = gsub("int\n","intercept\n",rhs)
    )
   
  cors$lty <- 3
  cors$lty[cors$pvalue < 0.05] <- 2
  cors$lty[cors$pvalue < 0.01] <- 1
  
  qgraph(cors[,1:3], directed = FALSE, theme = "colorblind", 
         edge.labels = TRUE, lty = cors$lty, vsize = 14)
  
  legend(legend.position, lty = c(1,2,3),legend = c("p < 0.01","p < 0.05","p > 0.05"), 
         bty = "n", cex = 1.5, lwd = 3)
  
}