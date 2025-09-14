human_SIS_ODE <- function(t, state, parameters) {
    H_S <- state[["H_S"]]
    H_I <- state[["H_I"]]
    N_H <- H_S + H_I
    N_V <- parameters[["N_V"]]
    I_V <- parameters[["I_V"]]
    
    # use inheritance cube to incorporate the
    # effect of genotype-specific transmission probability
    b <- parameters[["b"]]
    r <- parameters[["r"]]
    a <- parameters[["a"]]
    
    dH_S <- -a*H_S*(1/N_H)*(b%*%I_V) + r*H_I
    dH_I <- a*H_S*(1/N_H)*(b%*%I_V) - r*H_I 
    list(c(dH_S, dH_I))
}