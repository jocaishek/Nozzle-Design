import numpy as np
from scipy.optimize import fsolve

#initializers
pc = 3e6        # chamber pressure (Pa)
tc = 3500       # chamber temperature (K)
gamma = 1.4
r = 287         # J/kg-K; P=prt, r = P/pt
at = 0.01       # throat area (m^2)
ae = 0.05       # exit area (m^2)
pa = 101325     # sea level (Pa)
g0 = 9.81       
epsilon = ae/at     #how much nozzle expands after throat

#area
def area_mach(M):
    term1 = (2/(gamma+1)) * (1 + (gamma-1)/2 * M**2)
    return (1/M) * term1**((gamma+1)/(2*(gamma-1))) - epsilon
me_initial_guess = 2.0
me = fsolve(area_mach, me_initial_guess)[0]

#exit conditions
pe = pc * (1 + (gamma-1)/2 * me**2)**(-gamma/(gamma-1))
te = tc * (1 + (gamma-1)/2 * me**2)**(-1)
ve = me * np.sqrt(gamma * r * te)

#mass flow rate
mdot = (at * pc * np.sqrt(gamma/(r*tc)) * (2/(gamma+1))**((gamma+1)/(2*(gamma-1))))

#thrust
f = mdot * ve + (pe - pa) * ae

#impulse
imp = f / (mdot * g0)

#print statements
print("Exit Mach:", me)
print("Exit Velocity (m/s):", ve)
print("Mass Flow Rate (kg/s):", mdot)
print("Thrust (N):", f)
print("Specific Impulse (s):", imp)



