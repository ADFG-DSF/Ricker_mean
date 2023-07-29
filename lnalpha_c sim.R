library(ggplot2)
library(tidyverse)
library(magrittr)


# Functions ---------------------------------------------------------------
Ricker <- function(S, lnalpha, beta) {
  S*exp(lnalpha - beta*S)
}

#Simulate a single stock
simulateSR_goal <- function(lnalpha, beta, sigW, N, goal, sigS = 0, phi = 0, power = 1, sigF = 0) {
  stopifnot(length(goal) %in% 1:2)
  if(length(goal) == 1){stopifnot(goal >= 0 & goal < 1)}
  if(length(goal) == 2){stopifnot(goal[1] > 1 & goal[2] < (lnalpha + 0.5 * sigW * sigW) / beta & goal[2] > goal[1])}
  stopifnot(power >= 0 & power <=1) # fishing power is 0 for no harvest and 1 when the fishery can fully control escapement
  
  # ----- initial values ----- #
  # initial value for S: Seq minus some harvest
  S <- 900
  
  # initial value for observed S
  Shat <- S*rlnorm(1, sdlog=sigS)
  
  # initializing all other values
  redresid <- 0
  E1R <- E2R <- whiteresid <- epsF <- U1 <- U2 <- U <- Ft <- H <- Rhat <- lnRhatShat <- fittedR <- NA
  R <- 1 #This values in never used but needed to evaluate U
  
  # recursive portion...
  for(i in 2:(N + 50 +1)) { #added 50 (will discard later)
    E1R[i] <- S[i-1]*exp(lnalpha - beta*S[i-1])
    E2R[i] <- E1R[i]*exp(phi*redresid[i-1])
    R[i] <- E2R[i]*rlnorm(1,0,sigW)
    redresid[i] <- log(R[i]/E1R[i])
    whiteresid[i] <- log(R[i]/E2R[i])
    epsF[i] <- rnorm(1,0,sigF)
    if(length(goal) == 2){
      U1[i] <- if(R[i] > goal[1]){(R[i] - goal[1]) / R[i]} else{0} # Find U to fish to lb of goal
      U2[i] <- if(R[i] > goal[2]){(R[i] - goal[2]) / R[i]} else{0} # Find U to fish to ub of goal
    }
    U[i] <- if(length(goal) == 1){goal} else{min(runif(1, U2[i], U1[i]), power)} # draw a harvest rate
    Ft[i] <- -log(1 - U[i])*exp(epsF[i])
    S[i] <- R[i]*exp(-Ft[i])
    Shat[i] <- S[i]*rlnorm(1, sdlog=sigS)
    H[i] <- R[i]-S[i]
    Rhat[i] <- Shat[i]+H[i]
    lnRhatShat[i] <- log(Rhat[i]/Shat[i])
  }
  
  return(list(S=Shat[1:N + 50],
              R=Rhat[2:(N+1) + 50],
              Strue=S[1:N + 50],
              Rtrue=R[2:(N+1) + 50],
              yield=R[2:(N+1) + 50] - S[1:N + 50],
              goal = goal,
              lnalpha = lnalpha, 
              beta = beta,
              sigW = sigW))
}

#Function to find U that archives S somewhere within the goal
H_goal <- function(R, goal, power = .85, ...){
  # stopifnot(length(goal) == 4)
  # stopifnot(goal[1, 1] > 1 & goal[1, 2] < (lnalpha + 0.5 * sigW * sigW) / beta & goal[1, 2] > goal[1, 1])
  # stopifnot(goal[2, 1] > 1 & goal[2, 2] < (lnalpha + 0.5 * sigW * sigW) / beta& goal[2, 2] > goal[2, 1])
  U_lb <- U_ub <- U <- c(NA, NA)
  
  for(j in 1:2){
    U_lb[j] <- if(R[j] > goal[j, 1]){(R[j] - goal[j, 1]) / R[j]} else{0} # Find U to fish to lb of goal
    U_ub[j] <- if(R[j] > goal[j, 2]){(R[j] - goal[j, 2]) / R[j]} else{0} # Find U to fish to ub of goal
    U[j] <- min(runif(1, U_ub[j], U_lb[j]), power)
  }
  return(list(U, goal = goal, power = power))
}

#Function to find U that archives Smsy for each goal
H_Smsy <- function(R, power = 0.85){
  lna <- with(parent.frame(), lnalpha)
  b <- with(parent.frame(), beta)
  sig <- with(parent.frame(), sigW)
  lnap <- lna + 0.5 * sig * sig
  Smsy <- lna / b * (0.5 - 0.07 * lna)
  Smsyp <- lnap / b * (0.5 - 0.07 * lnap)
  
  U <- c(if(R[1] > Smsy){min((R[1] - Smsy) / R[1], power)} else{0},
         if(R[2] > Smsyp){min((R[2] - Smsyp) / R[2], power)} else{0})
  
  return(list(U, power = power))
}

#Function to find U that archives S a some percentage of the goal range
H_target <- function(R, goal, target, power = .85, ...){
  S_target <- U <- c(NA, NA)
  S_target <- c(goal[1,1] + target * (goal[1,2] - goal[1,1]), 
                goal[2,1] + target * (goal[2,2] - goal[2,1]))
  
  U <- c(if(R[1] > S_target[1]){min((R[1] - S_target[1]) / R[1], power)} else{0},
         if(R[2] > S_target[2]){min((R[2] - S_target[2]) / R[2], power)} else{0})
  
  return(list(U, goal = goal, power = power))
}

#Simulate 2 stocks w same dynamics but fished to different goals.
simSR_goals <- function(lnalpha, beta, sigW, N, sigF = 0, Hfun, ...){
  # ----- initial values ----- #
  # initial value for S: Seq minus some harvest
  S <- matrix(c(900, rep(NA, N + 50), 900, rep(NA, N + 50)), N + 51, 2)
  
  # initializing all other values
  redresid <- 0
  E1R <- E2R <- redresid <- whiteresid <- epsF <- U <- Ft <- H <- Rhat <- lnRhatShat <- fittedR <-
    matrix(c(rep(NA, N + 51), rep(NA, N + 51)), N + 51, 2)
  R <-  matrix(c(1, rep(NA, N + 50), 1, rep(NA, N + 50)), N + 51, 2)
  
  # recursive portion...
  for(i in 2:(N + 50 + 1)) { #added 50 (will discard later)
    E1R[i, ] <- S[(i-1), ]*exp(lnalpha - beta*S[(i-1), ]) 
    R[i, ] <- E1R[i, ]*rlnorm(1,0,sigW)
    epsF[i] <- rnorm(1,0,sigF)
    temp <- Hfun(R = R[i, ], ...)
    U[i,] <- temp[[1]]
    Ft[i, ] <- -log(1 - U[i, ])*exp(epsF[i])
    S[i, ] <- R[i, ]*exp(-Ft[i, ])
    H[i, ] <- R[i, ]-S[i, ]
  }
  return(data.frame(sim = rep(1:N, times = 2),
                    goal = rep(c("mode", "mean"), each = N),
                    S = c(S[1:N + 50, 1], S[1:N + 50, 2]),
                    R = c(R[2:(N+1) + 50, 1], R[2:(N+1) + 50, 2]),
                    Y = c(R[2:(N+1) + 50, 1] - S[1:N + 50, 1], R[2:(N+1) + 50, 2] - S[1:N + 50, 2]),
                    U = c(U[2:(N+1) +50, 1], U[2:(N+1) +50, 2]),
                    power = temp[["power"]],
                    lnalpha = lnalpha,
                    beta = beta,
                    sigW = sigW))
}

#table of realized yield near Smsy... Not sure I need this
tab_Smsy <- function(list){
  lnalpha_c <- list$lnalpha + .5 * list$sigW * list$sigW
  data.frame(median = median(list$yield[list$S < list$lnalpha / list$beta * (0.5 - 0.07 * list$lnalpha) + 50 &
                                          list$S > list$lnalpha / list$beta * (0.5 - 0.07 * list$lnalpha) - 50]),
             mean = mean(list$yield[list$S < lnalpha_c / list$beta * (0.5 - 0.07 * lnalpha_c) + 10 & 
                                      list$S > lnalpha_c / list$beta * (0.5 - 0.07 * lnalpha_c) - 10]))
}



# Meaning of bias correction ----------------------------------------------
SRout <- simulateSR_goal(lnalpha = 1.5, beta = 0.001, sigW = .5, N = 100000, goal = c(405, 889), power = 0.9)
SRout$group <-  
  as.numeric(
    as.character(
      cut(SRout$S - 25, seq(0, max(SRout$S) + 50, 50), labels = seq(from = 50, by = 50, length.out= length(seq(0, max(SRout$S), 50))))))
ggplot(data = as.data.frame(SRout), aes(x = group, y = R)) +
  geom_boxplot(aes(group = group), color = "black") +
  stat_function(fun = Ricker, args = list(lnalpha = SRout$lnalpha, beta = SRout$beta), aes(color = "black")) +
  stat_function(fun = Ricker, args = list(lnalpha = SRout$lnalpha + 0.5 * SRout$sigW * SRout$sigW, 
                                          beta = SRout$beta), aes(color = "red")) +
  geom_point(data = aggregate(R ~ group, SRout, mean), aes(y = R, color = "red")) +
  geom_abline() +
  annotate("rect", xmin = SRout$goal[1], xmax = SRout$goal[2], ymin = 0, ymax = Inf, alpha = 0.1) +
  scale_x_continuous(name = "S") + 
  scale_y_continuous(name = "R", limits = c(0, quantile(SRout$R, 0.99))) +
  scale_color_manual(values = c("black", "red"), labels = c("ln(\u03B1)", "ln(\u03B1')"), name = "Central Tendency") +
  theme_bw(base_size = 18)

# Params w bias correction ------------------------------------------------
#change in Smsy and MSY
diff <-
  expand.grid(lnalpha = seq(0.5, 2.5, length.out = 10), beta = 0.001, sigW = c(.25, .5, .75, 1)) %>%
  dplyr::mutate(Smsy = lnalpha / beta *(0.5 - 0.07 * lnalpha),
                Rmsy = Ricker(Smsy, lnalpha, beta),
                lnalpha_c = lnalpha + .5 * sigW * sigW,
                Smsy_c = lnalpha_c / beta *(0.5 - 0.07 * lnalpha_c),
                Rmsy_c = Ricker(Smsy_c, lnalpha, beta),
                diff_Smsy = Smsy_c / Smsy,
                diff_Rmsy = Rmsy_c / Rmsy,
                diff_lnalpha = lnalpha_c / lnalpha) %>%
  tidyr::pivot_longer(dplyr::starts_with("diff")) %>%
  dplyr::mutate(name = factor(name, 
                              levels = c("diff_lnalpha", "diff_Smsy", "diff_Rmsy"), 
                              labels = c("ln(alpha^c)/ln(alpha)", "S[msy]^{c}/S[msy]", "R[msy]^{c}/R[msy]")))
ggplot(diff, aes(x = lnalpha, y = value, color = as.character(sigW))) + 
  geom_line() + 
  geom_point() +
  xlab(expression("ln"~alpha)) +
  ylab("Ratio") +
  #coord_cartesian(xlim = c(0, 3), ylim = c(1:3)) +
  scale_color_discrete(name = expression(sigma)) +
  theme_bw(base_size = 20) +
  facet_grid(. ~ name, labeller = label_parsed)

# SR params & goals ----------------------------------------------------
lnalpha <- rep(seq(.5, 2.5, length.out = 3), times = 3)
beta0 <- 0.001
beta <- rep(beta0, 9)
sigW <- rep(seq(0.5, 1, length.out = 3), each = 3)
N0 <- 10000
N <- rep(N0, 9)
lnalpha_c <- lnalpha + 0.5 * sigW * sigW
Smsy <- lnalpha / 0.001 * (0.5 - 0.07 * lnalpha)
Smsy_c <- lnalpha_c / 0.001 * (0.5 - 0.07 * lnalpha_c)
median_eg <- mapply(function(x, y) range(which(Ricker(1:(x/beta0), x, beta0) - 1:(x/beta0) > (Ricker(y, x, beta0) - y) * 0.9)), 
                    x = lnalpha, 
                    y = Smsy, 
                    SIMPLIFY = FALSE)
mean_eg <- mapply(function(x, y) range(which(Ricker(1:(x/beta0), x, beta0) -  1:(x/beta0) > (Ricker(y, x, beta0) - y) * 0.9)), 
                  x = lnalpha_c, 
                  y = Smsy_c, 
                  SIMPLIFY = FALSE)



# * Table/plot of goals --------------------------------------------------------
#Table of goals
data.frame(lnalpha = lnalpha, sigW = sigW, 
           median_goal = sapply(median_eg, function(x) paste0(x[1], "-", x[2])),
           mean_goal = sapply(mean_eg, function(x) paste0(x[1], "-", x[2])))

#Plot of goals
curves <- 
  expand.grid(S = seq(0, 1500, 10), lnalpha = lnalpha, sigW = sigW, beta = beta) %>%
  mutate(median = Ricker(S, lnalpha, beta) - S,
         mean = Ricker(S, lnalpha + .5*sigW*sigW, beta) - S) %>%
  pivot_longer(median:mean, names_to = "curve_type", values_to = "Y") %>%
  filter(Y > 0)

data.frame(lnalpha = lnalpha, sigW = sigW, 
           median_goal_lb = sapply(median_eg, function(x) x[1]),
           median_goal_ub = sapply(median_eg, function(x) x[2]),
           mean_goal_lb = sapply(mean_eg, function(x) x[1]),
           mean_goal_ub = sapply(mean_eg, function(x) x[2])) %>%
  pivot_longer(ends_with("b"), names_to = c("goal", "bound"), names_pattern = "(.*)_.*_(.*)", values_to = "val") %>%
  pivot_wider(names_from = "bound", values_from = "val") %>%
  ggplot(aes(xmin = lb, xmax = ub, ymin = 0, ymax = Inf, fill = goal)) +
    geom_rect(alpha = 0.2) +
    geom_line(aes(x = S, y = Y, color = curve_type), curves, inherit.aes = FALSE) +
    theme_bw() +
    scale_fill_manual(values = c("yellow", "blue")) +
    facet_grid(lnalpha ~ sigW, scales = "free", labeller = label_both)

#Make a list of goals of each type
goals_list <- mapply(function(x, y) matrix(c(unlist(x), unlist(y)), 2, 2, byrow = TRUE), median_eg, mean_eg, SIMPLIFY = FALSE)


# Simulate Smsy fishing  ---------------------------------------------------------------
sim_Smsy <- 
  mapply(simSR_goals, lnalpha = lnalpha, beta = beta,  sigW = sigW, MoreArgs = list(Hfun = H_Smsy, N = 500), SIMPLIFY = FALSE) %>%
  do.call(rbind, .) %>%
  mutate(close = ifelse(U == 0, 1, 0),
         cap = ifelse(U == power, 1, 0)) 

#Check harvest rates
sim_Smsy %>%
  ggplot(aes(x = S, y = lag(U), color = goal)) + 
  geom_point() +
  scale_x_continuous(limits = c(0, 5000)) +
  facet_grid(lnalpha~sigW, scales = "free_y", labeller = label_both)

#Mean R, S, and Y
hinge <- function(x){quantile(x, c(0.25, 0.5, 0.75)) %>% setNames(c("ymin", "y", "ymax"))}
sim_Smsy %>%
  pivot_longer(cols = c(S:Y), names_to = "stat", values_to = "value") %>%
  ggplot(aes(x = stat, y = value, color = goal)) + 
  stat_summary(fun = "mean", geom = "point", position = position_dodge(width = 0.75)) +
  stat_summary(fun.data = "hinge", geom = "crossbar", position = position_dodge(width = 0.75), linetype = 2, width = .25) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", position = position_dodge(width = 0.75)) +
  facet_grid(lnalpha~sigW, scales = "free_y", labeller = label_both)

#Mean % of generations where mode goal performs better than mean goal
#Not any advantage wiped away if the Ymode has to be %% grater than Ymean
sim_Smsy %>%
  pivot_wider(id_cols = c(sim, lnalpha, beta, sigW, power), names_from = "goal", values_from = c("R", "S", "Y")) %>%
  mutate(R_gt = ifelse(R_mode > 1.0*R_mean, 1, 0),
         S_gt = ifelse(S_mode > 1.0*S_mean, 1, 0),
         Y_gt = ifelse(Y_mode > 1.0*Y_mean, 1, 0)) %>%
  pivot_longer(cols = ends_with("gt"), names_to = "stat", values_to = "index") %>%
  ggplot(aes(x = lnalpha, y = index)) + 
  stat_summary(fun = "mean", geom = "point", position = position_dodge(width = 0.75)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", position = position_dodge(width = 0.75)) +
  geom_hline(aes(yintercept = 0.5)) +
  facet_grid(stat~sigW, labeller = label_both)

#% of time closed and fished at max power under each goal
sim_Smsy  %>%
  pivot_longer(cols = c(close, cap), names_to = "stat", values_to = "value") %>%
  ggplot(aes(x = stat, y = value, color = goal)) + 
  stat_summary(fun = "mean", geom = "point", position = position_dodge(width = 0.75)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", position = position_dodge(width = 0.75)) +
  facet_grid(lnalpha~sigW, scales = "free_y", labeller = label_both)




# Simulate target fishing  ---------------------------------------------------------------
#aim for bottom of goal
sim_target <- 
  mapply(simSR_goals, lnalpha = lnalpha, beta = beta,  sigW = sigW, goal = goals_list, 
         MoreArgs = list(Hfun = H_target, N = 500, target = .05, sigF = 0), 
         SIMPLIFY = FALSE) %>%
  do.call(rbind, .) %>%
  mutate(close = ifelse(U == 0, 1, 0),
         cap = ifelse(U == power, 1, 0)) 

#Check harvest rates
sim_target %>%
  ggplot(aes(x = S, y = lag(U), color = goal)) + 
  geom_point() +
  scale_x_continuous(limits = c(0, 5000)) +
  facet_grid(lnalpha~sigW, scales = "free_y", labeller = label_both)

#Mean R, S, and Y
sim_target %>%
  pivot_longer(cols = c(S:Y), names_to = "stat", values_to = "value") %>%
  ggplot(aes(x = stat, y = value, color = goal)) + 
  stat_summary(fun = "mean", geom = "point", position = position_dodge(width = 0.75)) +
  stat_summary(fun.data = "hinge", geom = "crossbar", position = position_dodge(width = 0.75), linetype = 2, width = .25) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", position = position_dodge(width = 0.75)) +
  facet_grid(lnalpha~sigW, scales = "free_y", labeller = label_both)

#% of time closed and fished at max power under each goal
sim_target  %>%
  pivot_longer(cols = c(close, cap), names_to = "stat", values_to = "value") %>%
  ggplot(aes(x = stat, y = value, color = goal)) + 
  stat_summary(fun = "mean", geom = "point", position = position_dodge(width = 0.75)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", position = position_dodge(width = 0.75)) +
  facet_grid(lnalpha~sigW, scales = "free_y", labeller = label_both)




# Simulate goal fishing  ---------------------------------------------------------------
sim_goal <- 
  mapply(simSR_goals, lnalpha = lnalpha, beta = beta,  sigW = sigW, goal = goals_list, 
         MoreArgs = list(Hfun = H_goal, N = 500), 
         SIMPLIFY = FALSE) %>%
  do.call(rbind, .) %>%
  mutate(close = ifelse(U == 0, 1, 0),
         cap = ifelse(U == power, 1, 0)) 

#Check harvest rates
sim_goal %>%
  ggplot(aes(x = S, y = lag(U), color = goal)) + 
  geom_point() +
  scale_x_continuous(limits = c(0, 5000)) +
  facet_grid(lnalpha~sigW, scales = "free_y", labeller = label_both)


#mean R, S, and Y under each goal
sim_goal %>%
  pivot_longer(cols = c(S:Y), names_to = "stat", values_to = "value") %>%
  ggplot(aes(x = stat, y = value, color = goal)) + 
  stat_summary(fun = "mean", geom = "point", position = position_dodge(width = 0.75)) +
  stat_summary(fun.data = "hinge", geom = "crossbar", position = position_dodge(width = 0.75), linetype = 2, width = .25) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", position = position_dodge(width = 0.75)) +
  facet_grid(lnalpha~sigW, scales = "free_y", labeller = label_both)

#% of time closed and fished at max power under each goal
sim_goal  %>%
  pivot_longer(cols = c(close, cap), names_to = "stat", values_to = "value") %>%
  ggplot(aes(x = stat, y = value, color = goal)) + 
  stat_summary(fun = "mean", geom = "point", position = position_dodge(width = 0.75)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", position = position_dodge(width = 0.75)) +
  facet_grid(lnalpha~sigW, scales = "free_y", labeller = label_both)



# Outlier example ---------------------------------------------------------
# A situation where the mode might be preferable
#dat_outlier <- sim_goal[sim_goal$lnalpha == 1.5 & sim_goal$sigW == 0.75, c("S", "R")][1:25, ]
#saveRDS(dat_outlier, "H:\\My Documents\\Collaboration\\Ricker mean\\dat_outlier.rds")
dat_outlier <- readRDS("H:\\My Documents\\Collaboration\\Ricker mean\\dat_outlier.rds")
dat_outlier
mod_outlier <- lm(log(dat_outlier$R/dat_outlier$S) ~ dat_outlier$S)
summary(mod_outlier)
mod_outlier2 <- lm(log(dat_outlier$R[-23]/dat_outlier$S[-23]) ~ dat_outlier$S[-23]) #estimation wo outlier
summary(mod_outlier2)

plot(50:1000, exp(1.5+.5*.75*.75)*50:1000*exp(-.001*50:1000), 
     type = "l", 
     ylim = c(0, 18000),
     xlab = "S",
     ylab = "R")
points(dat_outlier$S, dat_outlier$R)
lines(50:1000, exp(coef(mod_outlier)[1] + 0.5*sigma(mod_outlier)*sigma(mod_outlier))*50:1000*exp(coef(mod_outlier)[2]*50:1000), col = "red")
lines(50:1000, exp(1.6)*50:1000*exp(-.0007*50:1000), col = "green")
legend("topleft", 
       c("true Ricker mean",
         "estimated Ricker mean",
         "estiamted Ricker mode"),
       col = c("black", "red", "green"),
       lty = 1)
#lines(50:1000, exp(coef(mod_outlier2)[1] + 0.5*sigma(mod_outlier2)*sigma(mod_outlier2))*50:1000*exp(coef(mod_outlier2)[2]*50:1000), col = "blue")
