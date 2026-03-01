"""
Nozzle method to store values for different nozzles
"""

import numpy as np
from scipy.optimize import fsolve
#from scipy.optimize import fsolve


class Nozzle():
    def __init__(self, pc = 3e6,        # chamber pressure (Pa),
                tc = 3500,       # chamber temperature (K)
                gamma = 1.4,
                r = 287,         # J/kg-K; P=prt, r = P/pt
                at = 0.01,       # throat area (m^2)
                ae = 0.05,       # exit area (m^2)
                pa = 101325,
                theta_deg = 15,
                L = 0.5):
        
        self.pc = pc
        self.tc = tc
        self.gamma = gamma
        self.r = r
        self.at = at
        self.ae = ae
        self.pa = pa
        self.g0 = 9.81
        self.epsilon = ae/at     #how much nozzle expands after throat
        
        self.theta = np.radians(theta_deg)
        self.L = L

        self.compute_geometry()


    def compute_geometry(self):
        self.rt = np.sqrt(self.at / np.pi)
        self.re = self.rt + self.L * np.tan(self.theta)
        self.ae = np.pi * self.re**2
        self.epsilon = self.ae / self.at
        self.eta_div = np.cos(self.theta)  # divergence efficiency


    def solve(self, me_initial_guess = 2.0):
        self.me = fsolve(self.area_mach, me_initial_guess)[0]
        #exit conditions
        self.pe = self.pc * (1 + (self.gamma-1)/2 * self.me**2)**(-self.gamma/(self.gamma-1))
        self.te = self.tc * (1 + (self.gamma-1)/2 * self.me**2)**(-1)
        self.ve = self.me * np.sqrt(self.gamma * self.r * self.te)
        #mass flow rate
        self.mdot = (self.at * self.pc * np.sqrt(self.gamma/(self.r*self.tc)) * (2/(self.gamma+1))**((self.gamma+1)/(2*(self.gamma-1))))
        #thrust
        '''self.f = self.eta_div * self.mdot * self.ve + (self.pe - self.pa) * self.ae'''
        #impulse
        self.imp = self.f / (self.mdot * self.g0)
        
    def area_mach(self, M):
       term1 = (2/(self.gamma+1)) * (1 + (self.gamma-1)/2 * M**2)
       return (1/M) * term1**((self.gamma+1)/(2*(self.gamma-1))) - self.epsilon
    

if __name__ == "__main__":
    """
    Main function that is used to validate
    """
    import matplotlib.pyplot as plt

    nozzle1 = Nozzle()
    nasars25 = Nozzle(pc = 2.1e7, tc = 3550, gamma = 1.22, r = 360, at = 0.01, ae=0.69)
    merlin1dsea = Nozzle(pc=9.7e6, tc=3400, gamma=1.22, r=310, at=0.01, ae=0.16, pa=101325) #sea level
    merlin1dvac = Nozzle(pc=9.7e6, tc=3400, gamma=1.22, r=310, at=0.01, ae=1.65, pa=0)
    saturnf1 = Nozzle(pc=7.0e6, tc=3300, gamma=1.23, r=300, at=0.01, ae=0.16, pa=101325)
    raptorvacuum = Nozzle(pc=3.0e7, tc=3500, gamma=1.20, r=370, at=0.01, ae=2.0, pa=0)

    nozzle1.solve()
    nasars25.solve()
    merlin1dsea.solve()
    merlin1dvac.solve()
    saturnf1.solve()
    raptorvacuum.solve()

    
    print("Exit Mach:", nozzle1.me, "Exit Mach nasars25:", nasars25.me, "Exit Mach merlin1dsea:", merlin1dsea.me, "Exit Mach marlin1dvac:", merlin1dvac.me, "Exit Mach saturnf1:", saturnf1.me, "Exit Mach raptorvacuum:", raptorvacuum.me)
    print("Thrust (N):", nozzle1.f, "Thrust (N) nasars25:", nasars25.f, "Thrust (N) merlin1dsea:", merlin1dsea.f, "Thrust (N) merlin1dvac:", merlin1dvac.f, "Thrust (N) saturnf1:", saturnf1.f, "Thrust (N) raptorvacuum:", raptorvacuum.f)

    plt.bar([1,2,3],[nozzle1.me,nasars25.me,merlin1dsea.me])
    plt.show()

def optimality_score(nozzle):
    if nozzle.pa == 0:
        delta = nozzle.pe
    else:
        delta = abs(nozzle.pe - nozzle.pa) / nozzle.pa
    return nozzle.f / (1 + delta)

nozzles = [nozzle1, nasars25, merlin1dsea, merlin1dvac, saturnf1, raptorvacuum]

names = ["Nozzle1", "NASA RS-25", "Merlin 1D Sea-Level", "Merlin 1D Vacuum", "Saturn F-1", "Raptor Vacuum"]

def optimality_score(nozzle):
    # dimensionless mismatch: Pe/Pa for sea-level, or Pe/(Pe + 1e5) for vacuum
    if nozzle.pa == 0:
        # vacuum environment, normalize by a reference pressure ~1e5 Pa
        delta = nozzle.pe / 1e5
    else:
        # sea level
        delta = abs(nozzle.pe - nozzle.pa) / nozzle.pa

    # normalize thrust too
    f_norm = nozzle.f / 1e6  # scale F to millions of N for numerical stability
    score = f_norm / (1 + delta)
    return score

#scores
scores = [optimality_score(n) for n in nozzles]

#Make scores to 0-1
max_score = max(scores)
normalized_scores = [s / max_score for s in scores]

for name, score in zip(names, normalized_scores):
    print(f"{name}: {score:.2f}")
'''
def optimize_nozzle(base_nozzle):

    def objective(x):
        theta_deg, L = x

        if theta_deg <= 1 or theta_deg >= 40:
            return 1e9
        if L <= 0.1 or L >= 5:
            return 1e9

        test_nozzle = Nozzle(pc=base_nozzle.pc,
                             tc=base_nozzle.tc,
                             gamma=base_nozzle.gamma,
                             r=base_nozzle.r,
                             at=base_nozzle.at,
                             theta_deg=theta_deg,
                             L=L,
                             pa=base_nozzle.pa)

        test_nozzle.solve()

        return -test_nozzle.f   # negative for maximization

    result = minimize(objective, x0=[15, 0.5], method='Nelder-Mead')

    return result
'''