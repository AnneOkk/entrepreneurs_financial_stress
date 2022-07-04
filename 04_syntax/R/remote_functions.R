# Installing uninstalled packages and loading packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


# Reading in .sav data
read_in <- function(df = files, directory, df_list) {
  for (fname in df) {
    df_list[[fname]] <- haven::read_sav(paste0(directory ,fname), encoding = NULL, user_na = FALSE, col_select = NULL,skip = 0, n_max = Inf, .name_repair = "unique")
  }
  names(df_list) <- paste0("", gsub(".sav","",names(df_list)))
  ff <- df_list
}

# Creating nicely formatted correlation table
corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower", "none"),
                    result=c("none", "html", "latex")){
  #Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method)
  numformat <- function(val) { sub("^(-?)0.", "\\1.", sprintf("%.2f", val)) }
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    ")))
  ## trunctuate the correlation matrix to two decimal
  R <- numformat(R)
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- " "
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  else if(removeTriangle[1]=="none"){
    Rnew <- as.matrix(Rnew)
    Rnew <- as.data.frame(Rnew)
  }
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex")
  }
}


corstars_no_stars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower", "none"),
                             result=c("none", "html", "latex")){
  #Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    ")))
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R,sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  else if(removeTriangle[1]=="none"){
    Rnew <- as.matrix(Rnew)
    Rnew <- as.data.frame(Rnew)
  }
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex")
  }
}


# Getting F and corresponding p value from lm object
Regressionp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)}


# Pasting together two character vectors; if corresponding element in either vector is NA, return only element from
# vector with !NA. If both elements are NA, return NA.
paste3 <- function(...,sep="") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}
val<- paste3(c("a","b", "c", NA), c("A","B", NA, NA))
val

# Saving indirect effect of mediation model
indirectsaved <- function(comp_df, random) {
  d <- comp_df[random, ] #rndomize by row
  apath <- lm(NA. ~ finstrain_center, data = d)
  bpath <- lm(satis ~ NA., data = d)
  indirect <- apath$coefficients[2]* bpath$coefficients[2]
  return(indirect)
}

# Creating time since business foundation in days
timediff <- function() {
  year_n <- as.numeric(names(attr(obs_df$T1year_1,"labels")[match(obs_df$T1year_1,attr(obs_df$T1year_1,"labels"))]))
  month_n <- as.numeric(obs_df$T1month_1 %>% remove_all_labels(.))
  time <- as.yearmon(paste(year_n, month_n), "%Y %m") %>% as_date(.)
  record <- as.yearmon(obs_df$T1RecordedDate, "%Y %m") %>% as_date(.)
  difftime <- difftime(record, time, UTC,
                       units = c("days"))
  return(difftime)
}

# Creating attention fail output T1. This may then be used to remove attention fails.
attention_fail_t1 <- function() {
  attention_f <- obs_df[(obs_df$T1chal_3 != 4) | (obs_df$T1threat_4 != 4) | (obs_df$T1satis_4 != 2),]
  ID_vals <- data.frame(table(attention_f$T1ID))
  Rows_fails <- attention_f$T1ID %in% ID_vals[ID_vals$Freq > 0,1]
  Att_fails <- attention_f[Rows_fails,]
  return(Att_fails)
}


# Creating attention fail output T2. This may then be used to remove attention fails.
attention_fail_t2 <- function() {
  attention_f <- t2[(t2$T2_lifesat_4 != 2) | (t2$T2_vita_5 != 3) | (t2$T2_entresat_4 != 6),]
  ID_vals <- data.frame(table(attention_f$T2_PID))
  Rows_fails <- attention_f$T2_PID %in% ID_vals[ID_vals$Freq > 0,1]
  Att_fails <- attention_f[Rows_fails,]
  return(Att_fails)
}


# Returning summary object for lavaan latent mediation model for models with variables in measurement part loading
# on three items each. With control variables = T1employable + T1employees + T1timebuiss + T1occ.

lavaan_sem <- function(IV, mediator, DV, controls = F) {
  if (controls == T) {
    controlinput <- paste('+ T1employable + T1employees + T1timebuiss + T1occ', '\n')
  } else {
    controlinput <- paste("", "\n")
  }

  measurement_c <- function() {
    measurement_part <- function(...) {
      params <- list(...)
      stopifnot(length(params)%%2==0)
      lefts = params[seq(1,length(params), by=2)]
      rights = params[seq(2,length(params), by=2)]
      rights <- Map(paste, rights, collapse="+")
      paste(paste0(lefts, " =~", rights), collapse="\n", "\n")
    }

    if(controls == F){
      meas <- measurement_part(paste0(IV), c(paste0(IV,"_1"), paste0(IV, "_2"), paste0(IV,"_3")),
                               paste0(mediator), c(paste0(mediator,"_1"), paste0(mediator, "_2"), paste0(mediator,"_3")),
                               paste0(DV), c(paste0(DV,"_1"), paste0(DV, "_2"), paste0(DV,"_3")))
    } else {
      meas <- measurement_part(paste0(IV), c(paste0(IV,"_1"), paste0(IV, "_2"), paste0(IV,"_3")),
                               paste0(mediator), c(paste0(mediator,"_1"), paste0(mediator, "_2"), paste0(mediator,"_3")),
                               paste0(DV), c(paste0(DV,"_1"), paste0(DV, "_2"), paste0(DV,"_3")),
                               "T1employable", c("T1employable_1", "T1employable_2", "T1employable_3"))

    }
    meas
  }
  model<- paste(measurement_c(),
                DV, "~ c*", IV, controlinput,
                mediator, "~ a*", IV, "\n",
                DV, "~ b*", mediator, "\n",
                "ab := a*b", '\n',
                "total := c + (a*b)")
  fit <- sem(model, data = obs_df)
  summary(fit, fit.measures=TRUE)
}


# center x and y variables in polynomial regression analysis on their grand mean
Center <- function(longer_comp) {

  # center x and y on grand mean across both variables
  grand.M <- mean( c(longer_comp$vitality, longer_comp$learn), na.rm = TRUE )
  longer_comp$x.c <- longer_comp$vitality-grand.M
  longer_comp$y.c <- longer_comp$learn-grand.M


  # compute higher order terms based on the centered predictors
  longer_comp$x2.c <- longer_comp$x.c^2
  longer_comp$xy.c <- longer_comp$x.c*longer_comp$y.c
  longer_comp$y2.c <- longer_comp$y.c^2

  # save means and compute their higher order terms
  longer_comp$x.mean  <- mean( longer_comp$vitality, na.rm = TRUE )
  longer_comp$y.mean  <- mean( longer_comp$learn, na.rm = TRUE )
  longer_comp$x2.mean <- longer_comp$x.mean^2
  longer_comp$xy.mean <- longer_comp$x.mean*longer_comp$y.mean
  longer_comp$y2.mean <- longer_comp$y.mean^2

  longer_comp

}

# Function to compute the intraclass correlation coefficient (ICC)
compute_icc <- function(lmer_object){
  var_dat <- lmer_object %>% VarCorr %>% as.data.frame
  icc <- var_dat$vcov[1]/(var_dat$vcov[1]+var_dat$vcov[2])
  return(icc)
}

# Create standardized coefficients
stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- round(fixef(object)*sdx/sdy, 2)
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc))
}


# create a lavaan mediation model for the multi studies paper
model_lavaan <- function(IV, mediator, DV, dat = obs_df, controls = F, fullfit = TRUE, T1control = F,
                         T1controlvar = NULL, T1controlmediator = F, T1controlmediatorvar = NULL, modificationindices_show = FALSE) {
  measurement_c <- function(...) {
    measurement_part <- function(...) {
      params <- list(...)
      stopifnot(length(params)%%2==0)
      lefts = params[seq(1,length(params), by=2)]
      rights = params[seq(2,length(params), by=2)]
      rights <- Map(paste, rights, collapse="+")
      paste(paste0(lefts, " =~", rights), collapse="\n", "\n")
    }

    if(controls == FALSE & T1control == FALSE & T1controlmediator == FALSE) {
      meas <- measurement_part(paste0(IV), grep(IV, names(dat), value = TRUE),
                               paste0(mediator), grep(mediator, names(dat), value = TRUE),
                               paste0(DV), grep(DV, names(dat), value = TRUE)
                               )
    } else if (controls == TRUE & T1control == FALSE & T1controlmediator == FALSE) {
      meas <- measurement_part(paste0(IV), grep(IV, names(dat), value = TRUE),
                               paste0(mediator), grep(mediator, names(dat), value = TRUE),
                               paste0(DV), grep(DV, names(dat), value = TRUE),
                               "T1employable", c("T1employable_1", "T1employable_2", "T1employable_3")
                               )
    } else if (controls == FALSE & T1control == TRUE & T1controlmediator == FALSE) {
      meas <- measurement_part(paste0(IV), grep(IV, names(dat), value = TRUE),
                               paste0(mediator), grep(mediator, names(dat), value = TRUE),
                               paste0(DV), grep(DV, names(dat), value = TRUE),
                               paste0(T1controlvar), grep(T1controlvar, names(dat), value = T)
                               )
    } else if (controls == TRUE & T1control == TRUE & T1controlmediator == FALSE) {
      meas <- measurement_part(paste0(IV), grep(IV, names(dat), value = TRUE),
                               paste0(mediator), grep(mediator, names(dat), value = TRUE),
                               paste0(DV), grep(DV, names(dat), value = TRUE),
                               paste0(T1controlvar), grep(T1controlvar, names(dat), value = T),
                               "T1employable", c("T1employable_1", "T1employable_2", "T1employable_3")
                               )
    } else if (controls == FALSE & T1control == FALSE & T1controlmediator == TRUE) {
      meas <- measurement_part(paste0(IV), grep(IV, names(dat), value = TRUE),
                               paste0(mediator), grep(mediator, names(dat), value = TRUE),
                               paste0(DV), grep(DV, names(dat), value = TRUE),
                               paste0(T1controlmediatorvar), grep(T1controlmediatorvar, names(dat), value = T)
                               )
    } else {
      meas <- measurement_part(paste0(IV), grep(IV, names(dat), value = TRUE),
                               paste0(mediator), grep(mediator, names(dat), value = TRUE),
                               paste0(DV), grep(DV, names(dat), value = TRUE),
                               paste0(T1controlmediatorvar), grep(T1controlmediatorvar, names(dat), value = T),
                               "T1employable", c("T1employable_1", "T1employable_2", "T1employable_3")
      )
    }
    meas
  }

  if (controls == FALSE & T1control == FALSE & T1controlmediator == FALSE) {
    controlinput <- paste("", "\n")
    } else if (controls == TRUE & T1control == FALSE & T1controlmediator == FALSE) {
     controlinput <- paste('+ T1employable + T1employees + T1timebuiss + T1occ', '\n')
    } else if (controls == FALSE & T1control == TRUE & T1controlmediator == FALSE) {
      controlinput <- paste('+', T1controlvar, '\n')
    } else if (controls == TRUE & T1control == TRUE & T1controlmediator == FALSE) {
      controlinput <- paste('+', T1controlvar, '+ T1employable + T1employees + T1timebuiss + T1occ', '\n')
    } else if (controls == FALSE & T1control == FALSE & T1controlmediator == TRUE) {
      controlinputmediator <- paste('+', T1controlmediatorvar, '\n')
      controlinput <- paste("", '\n')
    }else {
      controlinput <- paste('+ T1employable + T1employees + T1timebuiss + T1occ', '\n')
      controlinputmediator <- paste('+', T1controlmediatorvar, '\n')
    }
  if (T1controlmediator == FALSE) {
  model<- paste(measurement_c(IV, mediator, DV),
                DV, "~ c*", IV, controlinput,
                mediator, "~ a*", IV, "\n",
                DV, "~ b*", mediator, "\n",
                "ab := a*b", '\n',
                "total := c + (a*b)")
  } else {
    model<- paste(measurement_c(IV, mediator, DV),
                  DV, "~ c*", IV, controlinput,
                  mediator, "~ a*", IV, controlinputmediator, "\n",
                  DV, "~ b*", mediator, "\n",
                  "ab := a*b", '\n',
                  "total := c + (a*b)")
  }

  fit <- sem(model, data = obs_df)
  if (fullfit == TRUE) {
    summary(fit, fit.measures=TRUE)
  } else {
    subset(parameterEstimates(fit), rhs == "a*b") %>% select(-lhs, -op, -label, -z)
  }

  }


# create a lavaan mediation model for the multi studies paper
model_lavaan_multigroup <- function(IV, mediator, DV, dat = obs_df, controls = F, fullfit = TRUE, modificationindices_show = FALSE) {
  measurement_c <- function(...) {
    measurement_part <- function(...) {
      params <- list(...)
      stopifnot(length(params)%%2==0)
      lefts = params[seq(1,length(params), by=2)]
      rights = params[seq(2,length(params), by=2)]
      rights <- Map(paste, rights, collapse="+")
      paste(paste0(lefts, " =~", rights), collapse="\n", "\n")
    }

    if(controls == F){
      meas <- measurement_part(paste0(IV), grep(IV, names(dat), value = TRUE),
                               paste0(mediator), grep(mediator, names(dat), value = TRUE),
                               paste0(DV), grep(DV, names(dat), value = TRUE))
    } else {
      meas <- measurement_part(paste0(IV), grep(IV, names(dat), value = TRUE),
                               paste0(mediator), grep(mediator, names(dat), value = TRUE),
                               paste0(DV), grep(DV, names(dat), value = TRUE),
                               paste0(control_l), grep(control_l, names(dat), value = TRUE)
                               )

    }
    meas
  }

  if (controls == T) {
    controlinput <- paste('+ T1employable + T1employees + T1timebuiss + T1occ', '\n')
  } else {
    controlinput <- paste("", "\n")
  }
  model<- paste(measurement_c(IV, mediator, DV),
                DV, "~ c(c1, c11)*", IV, controlinput,
                mediator, "~ c(a1, a11)*", IV, "\n",
                DV, "~ c(b1, b11)*", mediator, "\n",
                "indirect1 := a1 * b1", "\n",
                "indirect2 := a11 * b11", "\n",
                "total1 := c1 + (a1*b1)", "\n",
                "total2 := c11 + (a11*b11)")
  fit <- sem(model, data = obs_df, group = "T1employees")
  if (fullfit == TRUE) {
    summary(fit, fit.measures=TRUE)
  } else {
    rbind(subset(parameterEstimates(fit), rhs == "a1*b1") %>% select(-lhs, -op, -label, -z),
    subset(parameterEstimates(fit), rhs == "a11*b11") %>% select(-lhs, -op, -label, -z))
  }
  if (modificationindices_show == TRUE) {
    lavaan::modindices(fit, sort = TRUE, maximum.number = 5)
  } else {
    NULL
  }
}


# create a lavaan mediation model for the multi studies paper
model_lavaan_multimed <- function(IV, mediator1, mediator2, DV, dat = obs_df, controls = F, fullfit = TRUE, modificationindices_show = FALSE) {
  measurement_c <- function(...) {
    measurement_part <- function(...) {
      params <- list(...)
      stopifnot(length(params)%%2==0)
      lefts = params[seq(1,length(params), by=2)]
      rights = params[seq(2,length(params), by=2)]
      rights <- Map(paste, rights, collapse="+")
      paste(paste0(lefts, " =~", rights), collapse="\n", "\n")
    }

    if(controls == F){
      meas <- measurement_part(paste0(IV), grep(IV, names(dat), value = TRUE),
                               paste0(mediator1), grep(mediator1, names(dat), value = TRUE),
                               paste0(mediator2), grep(mediator2, names(dat), value = TRUE),
                               paste0(DV), grep(DV, names(dat), value = TRUE))
    } else {
      meas <- measurement_part(paste0(IV), grep(IV, names(dat), value = TRUE),
                               paste0(mediator1), grep(mediator1, names(dat), value = TRUE),
                               paste0(mediator2), grep(mediator2, names(dat), value = TRUE),
                               paste0(DV), grep(DV, names(dat), value = TRUE),
                               "T1employable", c("T1employable_1", "T1employable_2", "T1employable_3")
                               )

    }
    meas
  }

  if (controls == T) {
    controlinput <- paste('+ T1employable + T1employees + T1timebuiss + T1occ + T1gender', '\n')
  } else {
    controlinput <- paste("", "\n")
  }
  model<- paste(measurement_c(IV, mediator1, nediator2, DV),
                DV, "~ c*", IV, controlinput,
                mediator1, "~ a1*", IV, "\n",
                mediator2, "~ a2*", IV, "\n",
                DV, "~ b1*", mediator1, "\n",
                DV, "~ b2*", mediator2, "\n",
                "indirect1 := a1 * b1", '\n',
                "indirect2 := a2 * b2", '\n',
                "contrast := indirect1 - indirect2", '\n',
                "total1 := c + (a1*b1)", '\n',
                "total2 := c + (a2*b2)")

  fit <- sem(model, data = obs_df)
  if (fullfit == TRUE) {
    summary(fit, fit.measures=TRUE)
  }
  else {
    rbind(subset(parameterEstimates(fit), label == "indirect1") %>% select(-lhs, -op, -label, -z),
          subset(parameterEstimates(fit), label == "indirect2") %>% select(-lhs, -op, -label, -z))
  }
  if (modificationindices_show == TRUE) {
    lavaan::modindices(fit, sort = TRUE, maximum.number = 5)
  } else {
    NULL
  }
}


medmodel <- function(step) {
  fit_table <- tidy(step, conf.int = T)
  output <- rbind(fit_table[fit_table$label == "c", ],
                  fit_table[fit_table$label == "a1", ],
                  fit_table[fit_table$label == "b1", ],
                  fit_table[fit_table$label == "a2", ],
                  fit_table[fit_table$label == "b2", ],
                  fit_table[fit_table$label %like% "control", ],
                  fit_table[fit_table$label == "a1b1", ],
                  fit_table[fit_table$label == "a2b2", ],
                  fit_table[fit_table$label == "total", ],
                  fit_table[fit_table$label == "total2", ]) %>%
    dplyr::select(dplyr::matches("term|estimate|p.value|conf.low|conf.high|std.all")) %>%
    mutate(p.value = p.value  %>% as.numeric(.)  %>% round(., digits = 3) %>% sprintf("%.3f",.)  %>%  as.character(.) %>% gsub("0\\.", "\\.", .)) %>%
    mutate(conf.low = formatC(conf.low, digits = 2, form = "f")) %>%
    mutate(conf.high = formatC(conf.high, digits = 2, form = "f")) %>%
    mutate(estimate = round(estimate, digits = 2)  %>% sprintf("%.2f",.)) %>%
    mutate(std.all = round(std.all, digits = 2) %>% sprintf("%.2f",.) %>%  as.character(.) %>% gsub("0\\.", "\\.", .)) %>%
    mutate(CI = paste0("[",conf.low, ", ", conf.high,"]")) %>%
    dplyr::select(-matches("conf"))
  out_df <- output[, c(1,2,5, 3, 4)]
  colnames(out_df) <- c("Predictor", "\\textit{B}", "95\\%CI", "\\textit{p}", "\\beta")
  return(out_df)
}





