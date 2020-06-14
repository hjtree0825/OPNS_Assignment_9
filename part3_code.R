library(chebpol); library(np); library(EnvStats)

sim_data <- function(Q = 20,
                     mu = 1,
                     sigma = 1,
                     pr = 5,
                     pw = 1,
                     alpha = 0.3,
                     eta_q = 4,
                     beta = 0.95,
                     num_grid = 100,
                     threshold = 10^-5,
                     num_periods = 10^4) {
  
  inc_val_fun <- function(v1, v0){
    euler_const <- digamma(1)
    Vbar <- pmax(v1, v0) + log(1 + exp(pmin(v1, v0) - pmax(v1,v0))) + euler_const
  }
  
  Vbb_v <- rep(0,num_grid)
  
  
  delta <- 100
  while(delta > threshold){
    
    f <- function(i){
      integrate(g, i = 1, lower = 0, upper = i) %>%
        .$value
    }
    
    Vbb <- Vectorize(ipol(f, dims = num_grid, intervals = list(c(0,Q))))
    
    g <- function(i,s){
      v1 <- pr*s - pw*(Q-i+s) - alpha*(i-s) - eta_q + beta*Vbb(Q)
      v0 <- pr*s - alpha*(i-s) + beta*Vbb(i-s)
      inc_val_fun(v1, v0) * dlnormTrunc(s, mu, sigma, min = 0, max = i)
    }
    
    Vbb_n <- Vbb(seq(from = 0, to = Q, length.out = num_grid))
    
    delta <- max(abs(Vbb_v - Vbb_n))
    
    Vbb_v <- Vbb_n
  }
  
  CCP_calc <- function(df){
    df %>%
      mutate(
        v1 = pr*s - pw*(Q-i+s) - alpha*(i-s) - eta_q + beta * Vbb(Q),
        v0 = pr*s - alpha*(i-s) + beta * Vbb(i-s),
        f = v0 - v1,
        prob1 = pmax(.Machine$double.xmin, 1/(1+exp(f))),
        prob0 = pmax(.Machine$double.xmin, 1/(1+exp(-f))),
        action = ifelse(runif(1) > prob0, 1, 0)
      ) %>%
      select(demand, i, s, prob1, prob0, action)
  }
  
  data_within <- data.frame(
    demand = rlnorm(num_periods, mu, sigma),
    i = Q
  ) %>%
    mutate(
      s = min(i,demand),
      prob1 = 0.5,
      prob0 = 1 - prob1,
      action = 1
    )
  
  data_within[1,] <- CCP_calc(data_within[i,])

  for (ind in 2:num_periods){
    data_within$i[ind] <- ifelse(data_within$action[ind - 1] == 1, Q, i-s)
    data_within$s[ind] <- min(i,d)
    data_within[ind,] <- CCP_calc(data_within[i,])
  }

  data_within %>%
    saveRDS(paste0('simulated_data.RDS'))

}



estimation <- function(df,
                       Q = 10,
                       pr = 10, # start
                       pw = 1,
                       mu = 10, # start
                       sigma = 10, # start
                       alpha = 10, # start
                       eta_q = 10, # start
                       beta = 0.95){
  
  initial_val <- data.frame(
    mu = mu,
    sigma = sigma,
    pr = pr,
    alpha = alpha,
    eta_q = eta_q
  )
  
  dfn <- df %>%
    mutate(
      term_inv = i - s,
      CCP_hat = predict(np::npreg(action ~ term_inv))
    )
  
  likelihood <- function(param_vec){
    
    mu <- param_vec$mu
    sigma <- param_vec$sigma
    
    dfn %>%
      group_by(i,s) %>%
      mutate(vall = ifelse(i!=0, dlnormTrunc(x = s, mu, sigma, min = 0, max = i), NA)) %>%
      ungroup() %>%
      mutate(prob_hat = ifelse(action == 1, CCP_hat, 1 - CCP_hat)) %>%
      filter(!is.na(vall)) %>%
      summarise(sum(log(prob_hat * vall))) %>%
      prod(-1) %>%
      as.numeric
  }
  
  param_vec_n <- optim(initial_val[1:2], likelihood, method = "L-BFGS-B", lower = c(-Inf, 0), control = list(trace = TRUE))$par
  
  estimate <- function(param_vec){
    
    f_hat <- function(i,s, pr, alpha, eta_q){
      action <- ifelse(Q - i + s > 0, 1, 0)
      term_inven <- i - s
      eval <- (param_vec$pr*s - param_vec$pw*(param_vec$Q-i+s) - param_vec$alpha*(i-s) - action * param_vec$eta_q - log(predict(npreg(action ~ term_inven)))) * dlnormTrunc(s, param_vec_n[1], param_vec_n[2], min = 0, max = i)
    }
    
    f <- function(i){
      integrate(
        f_hat,
        i = i,
        pr = pr,
        alpha = alpha,
        eta_q = eta_q,
        lower = 0,
        upper = i
      ) %>%
        .$value
    }
    
    f_tilde <- Vectorize(ipol(f, dims = num_grid, intervals = list(c(0,Q))))
    
    v <- function(x){
      pw*(Q-x) + eta_q + beta*(f_tilde(x) - f_tilde(Q))
    }
    
    like_val <- dfn %>%
      mutate(
        gp0 = 1/(1 + exp(-v(term_inven))),
        gp1 = 1/(1 + exp(v(term_inven))),
        likelihood_val = ifelse(action == 1, gp1, gp0),
        loglk <- log(likelihood_val)
      ) %>%
      summarise(sum(loglk)) %>%
      prod(-1) %>%
      as.numeric()

  }
  
  param_vec_n[3:5] <- optim(initial_val[3:5], estimate, method = "L-BFGS-B", lower = c(0,0, -Inf), control = list(trace = TRUE))
  
  return(param_vec_n)
}
