"""
Full-featured realistic nozzle optimization
"""

import numpy as np
from scipy.optimize import fsolve, minimize
import matplotlib.pyplot as plt

# -------------------- NOZZLE CLASS -------------------- #
class Nozzle:
    def __init__(self, pc=3e6, tc=3500, gamma=1.4, r=287,
                 at=0.01, r_profile=None, pa=101325, L=0.5):
        self.pc = pc
        self.tc = tc
        self.gamma = gamma
        self.r = r
        self.pa = pa
        self.L = L
        self.f = 0

        if r_profile is None:
            # default: linear expansion from throat
            self.r_profile = np.linspace(np.sqrt(at/np.pi), np.sqrt(at/np.pi)+0.05, 10)
        else:
            self.r_profile = np.array(r_profile)

        self.compute_geometry()

    def compute_geometry(self):
        self.x = np.linspace(0, self.L, len(self.r_profile))
        self.A = np.pi*self.r_profile**2
        self.at = self.A[0]
        self.ae = self.A[-1]
        self.epsilon = self.ae/self.at
        self.eta_div = np.cos(np.arctan((self.r_profile[-1]-self.r_profile[0])/self.L))

    def area_mach(self, M):
        term1 = (2/(self.gamma+1))*(1+(self.gamma-1)/2*M**2)
        return (1/M)*term1**((self.gamma+1)/(2*(self.gamma-1)))-self.epsilon

    def solve(self, me_initial_guess=2.0):
        self.compute_geometry()
        self.me = fsolve(self.area_mach, me_initial_guess)[0]
        self.pe = self.pc*(1+(self.gamma-1)/2*self.me**2)**(-self.gamma/(self.gamma-1))
        self.te = self.tc*(1+(self.gamma-1)/2*self.me**2)**-1
        self.ve = self.me*np.sqrt(self.gamma*self.r*self.te)
        self.mdot = (self.at*self.pc*np.sqrt(self.gamma/(self.r*self.tc))*
                     (2/(self.gamma+1))**((self.gamma+1)/(2*(self.gamma-1))))
        self.f = self.eta_div*self.mdot*self.ve + (self.pe - self.pa)*self.ae
        self.imp = self.f/(self.mdot*9.81)


# -------------------- SCORING FUNCTION -------------------- #
def optimality_score(nozzle):
    if nozzle.pa == 0:
        delta = nozzle.pe / 1e5
    else:
        delta = abs(nozzle.pe - nozzle.pa)/nozzle.pa
    f_norm = nozzle.f/1e6
    return f_norm/(1+delta)


# -------------------- FULL-GEOMETRY OPTIMIZATION WITH REALISTIC CONSTRAINTS -------------------- #
def optimize_nozzle_full(base_nozzle):
    N = len(base_nozzle.r_profile)
    r_initial = base_nozzle.r_profile.copy()
    
    # realistic constraints
    r_min = np.sqrt(base_nozzle.at/np.pi)   # throat radius
    r_max = r_min * 50  # maximum realistic expansion ratio ~50

    def objective(r_profile):
        # 1. Throat fixed
        if abs(r_profile[0]-r_min) > 1e-6:
            return 1e9
        # 2. Monotonic expansion
        if np.any(np.diff(r_profile) < 0):
            return 1e9
        # 3. Min/Max radius
        if np.any(r_profile < r_min) or np.any(r_profile > r_max):
            return 1e9

        test_nozzle = Nozzle(pc=base_nozzle.pc,
                             tc=base_nozzle.tc,
                             gamma=base_nozzle.gamma,
                             r=base_nozzle.r,
                             at=base_nozzle.at,
                             r_profile=r_profile,
                             pa=base_nozzle.pa,
                             L=base_nozzle.L)
        test_nozzle.solve()
        return -test_nozzle.f  # maximize thrust

    result = minimize(objective, r_initial, method='Nelder-Mead', options={'maxiter':500})
    optimized_nozzle = Nozzle(pc=base_nozzle.pc,
                              tc=base_nozzle.tc,
                              gamma=base_nozzle.gamma,
                              r=base_nozzle.r,
                              at=base_nozzle.at,
                              r_profile=result.x,
                              pa=base_nozzle.pa,
                              L=base_nozzle.L)
    optimized_nozzle.solve()
    return optimized_nozzle


# -------------------- MAIN SCRIPT -------------------- #
if __name__ == "__main__":
    # Define multiple nozzles
    nozzle1 = Nozzle()
    nasars25 = Nozzle(pc=2.1e7, tc=3550, gamma=1.22, r=360, at=0.01, r_profile=np.linspace(0.056,0.15,10))
    merlin1dsea = Nozzle(pc=9.7e6, tc=3400, gamma=1.22, r=310, at=0.01, r_profile=np.linspace(0.056,0.08,10), pa=101325)
    merlin1dvac = Nozzle(pc=9.7e6, tc=3400, gamma=1.22, r=310, at=0.01, r_profile=np.linspace(0.056,0.2,10), pa=0)
    saturnf1 = Nozzle(pc=7e6, tc=3300, gamma=1.23, r=300, at=0.01, r_profile=np.linspace(0.056,0.08,10), pa=101325)
    raptorvacuum = Nozzle(pc=3e7, tc=3500, gamma=1.2, r=370, at=0.01, r_profile=np.linspace(0.056,0.22,10), pa=0)

    nozzles = [nozzle1, nasars25, merlin1dsea, merlin1dvac, saturnf1, raptorvacuum]
    names = ["Nozzle1","NASA RS-25","Merlin 1D Sea-Level","Merlin 1D Vacuum","Saturn F-1","Raptor Vacuum"]

    # Solve all nozzles
    for n in nozzles:
        n.solve()

    # Compute normalized scores
    scores = [optimality_score(n) for n in nozzles]
    max_score = max(scores)
    normalized_scores = [s/max_score for s in scores]

    print("Normalized Scores:")
    for name, score in zip(names, normalized_scores):
        print(f"{name}: {score:.2f}")

    # -------------------- CHOOSE NOZZLE TO OPTIMIZE -------------------- #
    selected_nozzle = nasars25   # <-- CHANGE THIS to optimize a different nozzle
    optimal_nozzle = optimize_nozzle_full(selected_nozzle)

    # -------------------- PLOTS -------------------- #
    # Contour comparison
    plt.figure(figsize=(8,4))
    plt.plot(selected_nozzle.x, selected_nozzle.r_profile, 'b-o', label='Original')
    plt.plot(optimal_nozzle.x, optimal_nozzle.r_profile, 'r-o', label='Optimized')
    plt.xlabel("Axial position (m)")
    plt.ylabel("Radius (m)")
    plt.title("Nozzle Contour: Original vs Optimized")
    plt.legend()
    plt.show()

    # Thrust comparison
    plt.figure(figsize=(6,4))
    plt.bar(["Original","Optimized"], [selected_nozzle.f, optimal_nozzle.f])
    plt.ylabel("Thrust (N)")
    plt.title("Thrust Comparison")
    plt.show()

    # Exit Mach comparison
    plt.figure(figsize=(6,4))
    plt.bar(["Original","Optimized"], [selected_nozzle.me, optimal_nozzle.me], color=['blue','orange'])
    plt.ylabel("Exit Mach")
    plt.title("Exit Mach Comparison")
    plt.show()

    # -------------------- RESULTS -------------------- #
    print("Original Thrust: {:.2f} N, Exit Mach: {:.2f}".format(selected_nozzle.f, selected_nozzle.me))
    print("Optimized Thrust: {:.2f} N, Exit Mach: {:.2f}".format(optimal_nozzle.f, optimal_nozzle.me))
    print("Thrust Improvement: {:.2f}%".format((optimal_nozzle.f-selected_nozzle.f)/selected_nozzle.f*100))