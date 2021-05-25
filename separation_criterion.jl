# Final Project - CHEME 7770 
# Expansion of the Bell et al. (1984) Thermodynamic Model for Cell-Cell Adhesion
# Part 2: Evaluation of a Criterion for Self-Organization

using Roots
using NLsolve
using Plots
gr(show = true)

# Model Parameters
kb = 1.38e-16                           # dyn-cm
T = 310.15                              # Kelvin
A1 = 2e-6                               # cm^2
A2 = 2e-6                               # cm^2
Amax = 1e-6                             # cm^2
N1t = 1e5                               # tot num receptors (cell 1)
N2t = 1e5                               # tot num receptors (cell 2)
KL = 1e-8                               # cm^2
L = 2e-6                                # cm
sp_const = 0.1                          # dyn/cm
Nk = 20                                 # num Kuhn segments in a glycocalyx polymer
Lk = 15e-7                              # Kuhn length, cm
Qg = 100                                # num charged residues per polymer
Cm = 150                                # ionic strength of media, mM
d = 15e-7                               # diameter of glycocalyx molecule, cm
vK = Lk^2*d                             # excluded vol, cm^3



# -------------- Compute Adhesion Free Energies for Each Cell-Cell Interaction --------------
# Calculate standard-state chemical potentials using Eq (13). 
# Here, I approximate μ1₀ = μ2₀ = μb₀, and call this constant mu.
mu = kb*T*(log(KL) - 1)

Cg1 = 500
Cg2 = 1000
Cg = Cg1 + Cg2
# Define function for finding root of Eqn (9) to solve for the equilibrium value of S
# x: S
function f(x)
    t1 = (-3*Cg*sp_const)/(8*Nk*Lk^2)
    t2 = (3*Cg*sp_const*L)/(8*Nk*Lk^2)
    t3 = -sp_const - (3*kb*T)/(4*Nk*Lk^2)
    t4 = L - 2*vK*Cg*Nk^2 - (Qg^2*Cg)/(2*Cm)
    t5 = 2*vK*Nk^2 + Qg^2/(2*Cm)
    t6 = 2*kb*T*vK*Cg^2*Nk^2 + (kb*T*Qg^2*Cg^2)/(2*Cm)
    return x^5*t1 + x^4*t2 + x^3*Cg*t3 + x^2*Cg*sp_const*t4 + x*Cg^2*sp_const*L*t5 + t6
end

S = find_zero(f, (0,1))


# Define function for calculating equilibrium values of Nb and Ac 
function eq_params(Cg, Qg=Qg)
    repulsion = kb*T*(Cg + ((3*S^2*Cg)/(8*Nk*Lk^2)) + ((2*vK*Cg^2*Nk^2)/S) + ((Qg^2*Cg^2)/(2*Cm*S)))                                                                                         
    K = KL*exp(-0.5*sp_const*(S-L)*(S-L)/(kb*T))
    zeta = A1*A2*repulsion/(kb*T*K)
    Nb = 0.5*(N1t + N2t - sqrt((N1t - N2t)^2 + 4*zeta))     #Eq (10)
    Ac = ((kb*T)/(2*repulsion))*(N1t + N2t - sqrt((N1t-N2t)^2 + 4*zeta))     #Eq (11)
    return Nb, Ac
end

# Calculate the free energies associated with each type of adhesion.
function free_energy(Cg, Nb, Ac, Qg=Qg)
    mu1 = mu + kb*T*log((N1t-Nb)/A1)
    mu2 = mu + kb*T*log((N2t-Nb)/A2)
    mub = mu + 0.5*sp_const*(S-L)^2 + kb*T*log(Nb/Ac)
    repulsion = kb*T*(Cg + ((3*S^2*Cg)/(8*Nk*Lk^2)) + ((2*vK*Cg^2*Nk^2)/S) + ((Qg^2*Cg^2)/(2*Cm*S)))
    return (N1t-Nb)*mu1 + (N2t-Nb)*mu2 + Nb*mub + Ac*kb*T*repulsion
end


# Type 1-1 adhesion
G_11 = free_energy(2*Cg1, eq_params(2*Cg1)[1], eq_params(2*Cg1)[2])
println("G_11 = ", G_11)

# Type 2-2 adhesion
G_22 = free_energy(2*Cg2, eq_params(2*Cg2)[1], eq_params(2*Cg2)[2])
println("G_22 = ", G_22)

# Average of the homotypic cell Interactions
G_avg = (G_11 + G_22)/2
println("G_avg = ", G_avg)

# Type 1-2 adhesion
G_12 = free_energy(Cg1+Cg2, eq_params(Cg1+Cg2)[1], eq_params(Cg1+Cg2)[2])
println("G_12 = ", G_12)