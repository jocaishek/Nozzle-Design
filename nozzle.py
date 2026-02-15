"""
Nozzle method to store values for different nozzles
"""

import numpy as np
from scipy.optimize import fsolve


class Nozzle():
    def __init__(self, pc = 3e6,        # chamber pressure (Pa),
                tc = 3500,       # chamber temperature (K)
                gamma = 1.4,
                r = 287,         # J/kg-K; P=prt, r = P/pt
                at = 0.01,       # throat area (m^2)
                ae = 0.05,       # exit area (m^2)
                pa = 101325):
        
        self.pc = pc
        self.tc = tc
        self.gamma = gamma
        self.r = r
        self.at = at
        self.ae = ae
        self.pa = pa
        self.g0 = 9.81
        self.epsilon = ae/at     #how much nozzle expands after throat

    def solve(self, me_initial_guess = 2.0):
        self.me = fsolve(self.area_mach, me_initial_guess)[0]
        #exit conditions
        self.pe = self.pc * (1 + (self.gamma-1)/2 * self.me**2)**(-self.gamma/(self.gamma-1))
        self.te = self.tc * (1 + (self.gamma-1)/2 * self.me**2)**(-1)
        self.ve = self.me * np.sqrt(self.gamma * self.r * self.te)
        #mass flow rate
        self.mdot = (self.at * self.pc * np.sqrt(self.gamma/(self.r*self.tc)) * (2/(self.gamma+1))**((self.gamma+1)/(2*(self.gamma-1))))
        #thrust
        self.f = self.mdot * self.ve + (self.pe - self.pa) * self.ae
        #impulse
        self.imp = self.f / (self.mdot * self.g0)
        return 

    def area_mach(self, M):
        term1 = (2/(self.gamma+1)) * (1 + (self.gamma-1)/2 * M**2)
        return (1/M) * term1**((self.gamma+1)/(2*(self.gamma-1))) - self.epsilon
    

if __name__ == "__main__":
    """
    Main function that is used to validate
    """
    import matplotlib.pyplot as plt

    nozzle1 = Nozzle()
    nozzle2 = Nozzle(at = .03)
    nozzle3 = Nozzle(at = .03, ae = .7)

    nozzle1.solve()
    nozzle2.solve()
    nozzle3.solve()
    
    print("Exit Mach:", nozzle1.me, "Exit Mach 2:", nozzle2.me, "Exit Mach 3:", nozzle3.me, )
    print("Thrust (N):", nozzle1.f, "Thrust (N) 2:", nozzle2.f, "Thrust (N) 3:", nozzle3.f)

    plt.bar([1,2,3],[nozzle1.me,nozzle2.me,nozzle3.me])
    plt.show()
    