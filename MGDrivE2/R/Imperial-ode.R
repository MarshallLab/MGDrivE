#' ODE describing the age-structured Imperial model
#' used in decoupled sampling, which will pass in values of I_V and
#' return the human states for usein the mosquito portion of the model
#' @name human_Imperial_ODE
#' @param t starting time of simulation
#' @param state distributon of disease states
#' @param parameters parameter set
#' @return matrix of disease states after integration
#' @export
human_Imperial_ODE <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

        # keep track of indices so magic numbers aren't flying around
        indices <- list()
        indices[["S"]] <- 1
        indices[["T"]] <- 2
        indices[["D"]] <- 3
        indices[["A"]] <- 4
        indices[["U"]] <- 5
        indices[["P"]] <- 6
        indices[["ICA"]] <- 7
        indices[["IB"]] <- 8
        indices[["ID"]] <- 9
        indices[["IVA"]] <- 10

        state <- matrix(state, ncol=length(indices))

        # human states
        S <- state[, indices[["S"]]]
        T <- state[, indices[["T"]]]
        D <- state[, indices[["D"]]]
        A <- state[, indices[["A"]]]
        U <- state[, indices[["U"]]]
        P <- state[, indices[["P"]]]
        H <- 1

        # immunity functions
        ICA <- state[, indices[["ICA"]]]
        IB <- state[, indices[["IB"]]]
        ID <- state[, indices[["ID"]]]
        IVA <- state[, indices[["IVA"]]]

        Y <- S + A + U

        # infected states
        infected <- A + U + T + D

        # human dynamics
        dS <- rep(0, na)
        dT <- rep(0, na)
        dD <- rep(0, na)
        dA <- rep(0, na)
        dU <- rep(0, na)
        dP <- rep(0, na)
        # immunity functions

        dIB <- rep(0, na)
        dID <- rep(0, na)
        FOI <- rep(0, na)
        FOIv <- 0
        
        I_Vtot <- sum(I_V)
        EIRd <- av0 * ((I_Vtot)/total_M) / omega
        EIR <- EIRd * foi_age * rel_foi

        # calculate the weighted average of onward infectivity
        # to humans based on circulating genotypes at the timestep
        # sort to ensure names match up

        I_V <- I_V[sort(names(I_V))]
        b0 <- b0[sort(names(b0))]
        b <- rep(0, na)

        if(I_Vtot != 0) {
            I_Vprop <- I_V/I_Vtot
            # dot product gives the weighted average - cast to vector to suppress R warning
            b <- as.vector(b0 %*% I_Vprop) * ((1 - b1) / (1 + (IB / IB0) ^ kB) + b1)
        }

        for (i in 1:(na)) {
            FOI[i] <- EIR[i] * (ifelse(IB[i] == 0, b0, b[i]))
        }

        rate_ibaq <- EIR / (EIR * uB + 1)
        rate_clinaq <- FOI / (FOI * uCA + 1)
        rate_detaq <- FOI / (FOI * uD + 1)
        rate_sevaq <- FOI / (FOI * uVA + 1)

        dICA <- rep(0, na)
        dIVA <- rep(0, na)

        IC <- ICM + ICA
        phi <- phi0 * ((1 - phi1) / (1 + (IC / IC0)^kC) + phi1)
        clin_inc <- phi * FOI * Y
        fv <- 1 - (1-fvS)/(1+(age_days/av)^gammaV)
        immunity_multiplier <- ((IVA+IVM)/iv0)^kv
        severe_disease <- theta0*(theta1 + ((1-theta1)/(1+fv*immunity_multiplier)))*FOI
        mort <- pctMort*severe_disease

        for (i in 1:(na)) {
            dS[i] <-
                -FOI[i] * S[i] + rP * P[i] + rU * U[i] + ifelse(i == 1, eta * H * het_wt, 0) - (eta + age_rate[i]) *
                    S[i] + ifelse(i == 1, 0, age_rate[i - 1] * S[i - 1])

            dT[i] <-
                ft * clin_inc[i] - rT * T[i] - (eta + age_rate[i]) * T[i] + ifelse(i == 1, 0, age_rate[i - 1] * T[i - 1])


            dD[i] <-
                (1 - ft) * clin_inc[i] - rD * D[i] - (eta + age_rate[i]) * D[i] + ifelse(i ==
                    1, 0, age_rate[i - 1] * D[i - 1])

            dA[i] <-
                (1 - phi[i]) * FOI[i] * Y[i] - FOI[i] * A[i] + rD * D[i] - recA *
                    A[i] - (eta + age_rate[i]) * A[i] + ifelse(i == 1, 0, age_rate[i - 1] * A[i - 1])

            dU[i] <-
                recA * A[i] - FOI[i] * U[i] - rU * U[i] - (eta + age_rate[i]) * U[i] + ifelse(i == 1, 0, age_rate[i - 1] * U[i - 1])
            dP[i] <-
                rT * T[i] - rP * P[i] - (eta + age_rate[i]) * P[i] + ifelse(i == 1, 0, age_rate[i - 1] * P[i - 1])

            dICA[i] <-
                rate_clinaq[i] - 1 / dc * ICA[i] + (ifelse(i == 1, -ICA[1] / x_I[1], -(ICA[i] - ICA[i - 1]) / x_I[i]))

            dIB[i] <-
                rate_ibaq[i] - IB[i] / db + (ifelse(i == 1, -IB[1] / x_I[1], -(IB[i] - IB[i -
                    1]) / x_I[i]))
            dID[i] <-
                rate_detaq[i] - ID[i] / dd + (ifelse(i == 1, -ID[1] / x_I[1], -(ID[i] - ID[i - 1]) / x_I[i]))

            dIVA[i] <- rate_sevaq[i] - IVA[i] / dv + (ifelse(i == 1, -IVA[1] / x_I[1], -(IVA[i] - IVA[i - 1]) / x_I[i]))
        }

        list(c(
            dS,
            dT,
            dD,
            dA,
            dU,
            dP,
            dICA,
            dIB,
            dID,
            dIVA
        ),
        clin_inc=clin_inc,
        mort=mort)
    })
}
