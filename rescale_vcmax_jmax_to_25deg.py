import math
from scipy.optimize import fsolve

def arrh_to_25(k25, Ea, Tk, RGAS, Tref=298.15):
    """ Standard Arrhenius correction """
    return k25 / math.exp((Ea * (Tk - Tref)) / (Tref * RGAS * Tk))

# Peaked Arrhenius forward function
def peaked_arrh(k25, Ea, Tk, deltaS, Hd, RGAS):
    arg1 = k25 * math.exp((Ea * (Tk - 298.15)) / (298.15 * RGAS * Tk))
    arg2 = 1.0 + math.exp((298.15 * deltaS - Hd) / (298.15 * RGAS))
    arg3 = 1.0 + math.exp((Tk * deltaS - Hd) / (Tk * RGAS))
    return arg1 * arg2 / arg3

# Peaked Arrhenius inversion to find k25
def peaked_arrh_to_25(V_T, Ea, Tk, deltaS, Hd, RGAS):
    """Invert peaked Arrhenius to get k25 from measurement at Tk."""
    def func(k25):
        return peaked_arrh(k25, Ea, Tk, deltaS, Hd, RGAS) - V_T
    k25_guess = V_T  # initial guess
    k25_solution, = fsolve(func, k25_guess)
    return k25_solution

# Constants
RGAS = 8.314  # J mol-1 K-1

# Activation energies (J mol-1)
Eav = 60000.0
Eaj = 30000.0

# Deactivation energies (J mol-1)
Hdv = 200000.0
Hdj = 200000.0

# Entropy term (J mol-1 K-1)
Svv = 650.0
Svj = 650.0

# Leaf temperature
Tleaf_C = 23.0
Tleaf_K = Tleaf_C + 273.15
Tref_K = 298.15  # 25°C

# Observed Vcmax and Jmax at leaf temperature
VcmaxT = 40.0
JmaxT = 80.0

# If your dataset is mostly around 15-30 deg C, standard Arrhenius is
# probably fine.
Vcmax25 = arrh_to_25(VcmaxT, Eav, Tleaf_K, RGAS)
Jmax25 = arrh_to_25(JmaxT, Eaj, Tleaf_K, RGAS)

# Peaked version accounts for the high-temperature decline in Vcmax,
# better for leaf temperatures >30 degC, but needs extra parameters
Vcmax25_peak = peaked_arrh_to_25(VcmaxT, Eav, Tleaf_K, Svv, Hdv,  RGAS)
Jmax25_peak = peaked_arrh_to_25(JmaxT, Eaj, Tleaf_K, Svj, Hdj,  RGAS)

print("Measurement at leaf T")
print("VcmaxT =", VcmaxT)
print("JmaxT =", JmaxT)

print("\nCorrected to 25 deg, standard Arrhenius:")
print("Vcmax25 =", Vcmax25)
print("Jmax25 =", Jmax25)

print("\nCorrected to 25 deg, peaked Arrhenius:")
print("Vcmax25_peak =", Vcmax25_peak)
print("Jmax25_peak =", Jmax25_peak)
