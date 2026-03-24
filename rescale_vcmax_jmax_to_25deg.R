# Standard Arrhenius correction
arrh_to_25 <- function(k25, Ea, Tk, RGAS) {
  k25 / exp((Ea * (Tk - 298.15)) / (298.15 * RGAS * Tk))
}

# Forward peaked Arrhenius
peaked_arrh <- function(k25, Ea, Tk, deltaS, Hd, RGAS) {
  arg1 <- k25 * exp((Ea * (Tk - 298.15)) / (298.15 * RGAS * Tk))
  arg2 <- 1 + exp((298.15 * deltaS - Hd) / (298.15 * RGAS))
  arg3 <- 1 + exp((Tk * deltaS - Hd) / (Tk * RGAS))
  return(arg1 * arg2 / arg3)
}

# Inversion to find k25
peaked_arrh_to_25 <- function(V_T, Ea, Tk, deltaS, Hd, RGAS) {
  f <- function(k25) peaked_arrh(k25, Ea, Tk, deltaS, Hd, RGAS) - V_T
  # initial guess = measured value
  result <- uniroot(f, c(0.5 * V_T, 2 * V_T))$root
  return(result)
}

# Constants
RGAS <- 8.314  # J mol-1 K-1

# Activation energies (J mol-1)
Eav <- 60000.0
Eaj <- 30000.0

# Deactivation energies (J mol-1)
Hdv <- 200000.0
Hdj <- 200000.0

# Entropy term (J mol-1 K-1)
Svv <- 650.0
Svj <- 650.0

# Leaf temperature
Tleaf_C <- 23.0
Tleaf_K <- Tleaf_C + 273.15
Tref_K <- 298.15  # 25 deg C

# Observed Vcmax and Jmax at leaf temperature
VcmaxT <- 40.0
JmaxT <- 80.0

# If your dataset is mostly around 15-30 deg C, standard Arrhenius is
# probably fine.
Vcmax25 <- arrh_to_25(VcmaxT, Eav, Tleaf_K, RGAS)
Jmax25 <- arrh_to_25(JmaxT, Eaj, Tleaf_K, RGAS)

# Peaked version accounts for the high-temperature decline in Vcmax,
# better for leaf temperatures >30 degC, but needs extra parameters
Vcmax25_peak <- peaked_arrh_to_25(VcmaxT, Eav, Tleaf_K, Svv, Hdv, RGAS)
Jmax25_peak <- peaked_arrh_to_25(JmaxT, Eaj, Tleaf_K, Svj, Hdj, RGAS)

cat("Measurement at leaf T\n")
cat("VcmaxT =", VcmaxT, "\n")
cat("JmaxT =", JmaxT, "\n\n")

cat("Corrected to 25 deg, standard Arrhenius:\n")
cat("Vcmax25 =", Vcmax25, "\n")
cat("Jmax25 =", Jmax25, "\n\n")

cat("\nCorrected to 25 deg, peaked Arrhenius:\n")
cat("Vcmax25_peak =", Vcmax25_peak, "\n")
cat("Jmax25_peak =", Jmax25_peak, "\n")
