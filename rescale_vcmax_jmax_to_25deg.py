import math

def arrh(k25, Ea, Tk, RGAS):
    """ Temperature dependence of kinetic parameters is described by an
    Arrhenius function.

    Parameters:
    ----------
    k25 : float
        rate parameter value at 25 degC or 298 K
    Ea : float
        activation energy for the parameter [J mol-1]
    Tk : float
        leaf temperature [deg K]

    Returns:
    -------
    kt : float
        temperature dependence on parameter

    References:
    -----------
    * Medlyn et al. 2002, PCE, 25, 1167-1179.
    """
    return k25 * math.exp((Ea * (Tk - 298.15)) / (298.15 * RGAS * Tk))


def peaked_arrh(k25, Ea, Tk, deltaS, Hd, RGAS):
    """ Temperature dependancy approximated by peaked Arrhenius eqn,
    accounting for the rate of inhibition at higher temperatures.

    Parameters:
    ----------
    k25 : float
        rate parameter value at 25 degC or 298 K
    Ea : float
        activation energy for the parameter [J mol-1]
    Tk : float
        leaf temperature [deg K]
    deltaS : float
        entropy factor [J mol-1 K-1)
    Hd : float
        describes rate of decrease about the optimum temp [J mol-1]

    Returns:
    -------
    kt : float
        temperature dependence on parameter

    References:
    -----------
    * Medlyn et al. 2002, PCE, 25, 1167-1179.

    """
    arg1 = arrh(k25, Ea, Tk, RGAS)
    arg2 = 1.0 + math.exp((298.15 * deltaS - Hd) / (298.15 * RGAS))
    arg3 = 1.0 + math.exp((Tk * deltaS - Hd) / (Tk * RGAS))

    return arg1 * arg2 / arg3


# Constants
RGAS = 8.314  # J mol-1 K-1

# Activation energies (J/mol)
Eav = 65300.0  # Vcmax
Eaj = 43790.0  # Jmax

# Deactivation energies (J/mol)
Hdv = 200000.0
Hdj = 200000.0

# Entropy term (J mol-1 K-1)
Svv = 650.0
Svj = 650.0

# Leaf temperature
Tleaf_C = 20.0
Tleaf_K = Tleaf_C + 273.15
Tref_K = 298.15  # 25°C

# Observed Vcmax and Jmax at leaf temperature
VcmaxT = 40.0
JmaxT = 80.0

# If your dataset is mostly around 15-30 deg C, standard Arrhenius is
# probably fine.
Vcmax25 = arrh(VcmaxT, Eav, Tleaf_K, RGAS)
Jmax25 = arrh(JmaxT, Eaj, Tleaf_K, RGAS)

# Peaked version accounts for the high-temperature decline in Vcmax,
# better for leaf temperatures >30 degC, but needs extra parameters 
Vcmax25_peak = peaked_arrh(VcmaxT, Eav, Tleaf_K, Svv, Hdv,  RGAS)
Jmax25_peak = peaked_arrh(JmaxT, Eaj, Tleaf_K, Svj, Hdj,  RGAS)

print("Measurement at leaf T")
print("VcmaxT =", VcmaxT)
print("JmaxT =", JmaxT)

print("\nCorrected to 25 deg, standard Arrhenius:")
print("Vcmax25 =", Vcmax25)
print("Jmax25 =", Jmax25)

print("\nCorrected to 25 deg, peaked Arrhenius:")
print("Vcmax25_peak =", Vcmax25_peak)
print("Jmax25_peak =", Jmax25_peak)
