
# Constants
RGAS <- 8.314   # J mol-1 K-1
Eav <- 65300.0  # J mol-1 (Bernacchi et al. 2004 PCE)

Tleaf_C <- 10.0
VcmaxT <- 22.0

# Convert to Kelvin
Tleaf_K <- Tleaf_C + 273.15
Tref_K <- 298.15  # 25 deg C

Vcmax25 <- VcmaxT * exp(Eav * (Tref_K - Tleaf_K) / (RGAS * Tleaf_K * Tref_K))


print(Vcmax25)
