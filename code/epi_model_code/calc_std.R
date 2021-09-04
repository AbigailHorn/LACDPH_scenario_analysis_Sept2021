

##############################################################################
## CALCULATING STANDARD DEVIATION OF R*mu

# Inputs
R0_redux1 <- .21
R0_redux2 <- .33
R0_redux3 <- .55 #.575#.6 #65

prior.R0.redux1 <- c("unif", R0_redux1 - 0.05, R0_redux1 + 0.05) # .11
prior.R0.redux2 <- c("unif", R0_redux2 - 0.05, R0_redux2 + 0.05) # .11
prior.R0.redux3 <- c("unif", R0_redux3 - 0.05, R0_redux3 + 0.05) #.12

# Calculation
mu = R0_redux3*.65
b = mu + .05
a = mu - .05
var_mu = ((b-a)^2)/12

var_R <- .1^2

var_R_mu <- var_R*var_mu + var_R*(mu^2)+var_mu*(R0^2)

stdev_R_mu <- var_R_mu^.5 %>% round(digits=2)
stdev_R_mu

R0_y_print <- R0_y %>% round(digits=2)

## Output for paper table
R0_output <- paste0("$~N(",R0_y_print,",",stdev_R_mu,")$")
R0_output


## r(t)
prior.r1 <- c("unif",0.03, 0.35)
prior.r2 <- c("unif", 0.3, 0.65)

b = .65 #0.35
a = .35 #.03
mu.r <- (b+a)/2
mu.r
stdev.r = (((b-a)^2)/12)^.5
stdev.r
##############################################################################
## AKD

## Output for paper table
stdev=0.01
Alpha_output <- paste0("$~N(",Alpha_y,",",stdev,")$")
Kappa_output <- paste0("$~N(",Kappa_y,",",stdev,")$")
Delta_output <- paste0("$~N(",Delta_y,",",stdev,")$")

