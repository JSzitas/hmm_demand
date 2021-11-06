

library(magrittr)

intermittent_gen <- function(.length = 100, .prob = 0.5, .lambda = 150) {
  rpois(.length, lambda = .lambda) * rbinom(n = .length, size = 1, prob = .prob)
}

intermittent_example <- intermittent_gen(.prob = 0.3)
intermittent_short_example <- intermittent_gen(20, .prob = 0.4)

# idea: predict zero vs nonzero element of series as a regime switch via HMM:
hmm_zero_model <- function(series, lags = 1:5) {
  df <- tibble::tibble(target = series)  %>% 
    dplyr::mutate( target = target != 0 )
  
  df <- purrr::map(lags, ~ dplyr::lag(df[["target"]], n = .x)) %>%
    magrittr::set_names( paste0("lag_",lags) ) %>% 
    dplyr::bind_cols() %>%
    dplyr::bind_cols(df) %>% 
    na.omit() %>% 
    dplyr::mutate( target = as.factor(target) ) %>%
    as.data.frame()

  hmm <- depmixS4::depmix( target ~ ., data = df, nstates = 2, family = binomial(),
                           emcontrol = depmixS4::em.control(maxit = 1000, classification = "hard")) 

  fit_hmm <- depmixS4::fit(hmm)

  return(fit_hmm)
}

check_hmm_acc <- function( series, fitted_hmm ) {
  
  success_states_from_hmm <- as.logical(fitted_hmm@posterior[,1] -1)
  success_actual_series <- (series == 0)[
    (length(series) - length(success_states_from_hmm) + 1):length(series) ]
  
  tbl <- table( success_states_from_hmm, 
                success_actual_series )
  
  sum(diag(tbl))/sum(tbl)
}

forecast_hmm <- function( object, new_data = NULL, h = NULL, ... ) {

  transition_mat <- rbind(depmixS4::getpars(depmixS4::getmodel(object,"transition",1)),
                          depmixS4::getpars(depmixS4::getmodel(object,"transition",2))
                          )
  # extract the probability of the states at the final time point in the data (t=T)
  # this will act as a "prior" to compute the forecasted state distributions
  posterior_mat <- depmixS4::posterior(object, type = "viterbi")
  prior <- as.numeric( posterior_mat[nrow(posterior_mat),-1])

  # state-wise predictions for the observed variables
  coef <- rbind( 
    depmixS4::getpars(depmixS4::getmodel(object,"response",1)),
    depmixS4::getpars(depmixS4::getmodel(object,"response",2))
    )

  df <- object@response[[1]][[1]]@x
  df <- df[nrow(df),]
  project <- coef %*% df
  
  # return(prior %*% transition_mat)
  
  # return(1/(1+exp(project)))
  # return(1/(1+exp(-project)) )
  # return(1/(1+exp(-project)) )
  
  
  return( c(1/(1+exp(-project))) %*% c(prior %*% transition_mat)) 
  # return( pred_r_by_state * (prior %*% transition_mat)) 
  # return(df)
  
  
  
  
  
  # return(list(pred_r_by_state, prior %*% transition_mat))
  return(sum(df %*% pred_r_by_state * (prior %*% transition_mat)))
  # return(prior %*% transition_mat)
  # fcst <- rep(NA,h)#matrix( NA, ncol = 2, nrow = h )
  # 
  # mat_pow <- function( A, power = 1 ) {
  #   if( power == 1 ) return(A)
  # 
  #   B <- A
  #   for( i in seq_len(power) ) {
  #     B <- A %*% B
  #   }
  #   return(B)
  # }
  # 
  # for( step in seq_len(h) ) {
  #   current_transition <- mat_pow( transition_mat, step )
  #   fcst[step] <- log(sum(t(pred_r_by_state %*% t(prior %*% current_transition))))/2
  # }
  # return(fcst)
}

forecast_hmm2 <- function( object, h = 8 ) {
  transition_mat <- rbind(depmixS4::getpars(depmixS4::getmodel(object,"transition",1)),
                          depmixS4::getpars(depmixS4::getmodel(object,"transition",2))
  )
  posterior <- chk@posterior[["state"]]
  prior <- posterior[length(posterior)] - 1
  
  #define the Markov chain
  states <- c(0, 1)
  state_names <- as.character(states)
  mc <- new("markovchain", states = state_names,
            transitionMatrix = transition_mat, 
            # nrow = 2, 
            byrow = TRUE
            # dimnames = list(state_names, state_names)
  )
  as.numeric(markovchain::markovchainSequence(n = h, markovchain = mc,
                                   t0 = state_names[ prior == states ] ))
}

# data.frame( target = intermittent_example ) %>%
#   dplyr::mutate( lags = purrr::map( 1:5,  ) )
# hmm <- depmixS4::depmix( target ~ ., data = df_prepared$train, nstates = 2 )
#
# fit_hmm <- depmixS4::fit(hmm)
#
# posterior <- depmixS4::posterior(fit_hmm, type = "viterbi")
# depmixS4::logLik(fit_hmm)
