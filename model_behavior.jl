# Final Project - CHEME 7770 
# Expansion of the Bell et al. (1984) Thermodynamic Model for Cell-Cell Adhesion 
# Part 1: Examining the behavior of the extended thermodynamic model
# Bell's model was extended by replacing the simplified phenomenological equation with Equation (10) from Gandhi et al. (2019).

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
Cg1 = 500                               # spatially avgd glycocalyx number density on Cell 1, num/cm^2
Cg2 = 500                               # spatially avgd glycocalyx number density on Cell 2, num/cm^2
Qg = 100                                # num charged residues per polymer
Cm = 150                                # ionic strength of media, mM
d = 15e-7                               # diameter of glycocalyx molecule, cm
vK = Lk^2*d                             # excluded vol, cm^3
Cg = Cg1 + Cg2


# -------------------- Calculate Equilibrium Values of (Cg1+Cg2), Nb, and Ac --------------------
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
println("S_eq = ", S, " cm")


# Define function for calculating equilibrium values of Nb and Ac 
function eq_params(Cg, Qg=Qg)
    repulsion = kb*T*(Cg + ((3*S^2*Cg)/(8*Nk*Lk^2)) + ((2*vK*Cg^2*Nk^2)/S) + ((Qg^2*Cg^2)/(2*Cm*S)))                                                                                         
    K = KL*exp(-0.5*sp_const*(S-L)*(S-L)/(kb*T))
    zeta = A1*A2*repulsion/(kb*T*K)
    Nb = 0.5*(N1t + N2t - sqrt((N1t - N2t)^2 + 4*zeta))    #Eq. (10)
    Ac = ((kb*T)/(2*repulsion))*(N1t + N2t - sqrt((N1t-N2t)^2 + 4*zeta))   #Eq. (11)
    return Nb, Ac
end

# Calculate equilibrium values of Nb and Ac
Nb_eq = eq_params(Cg)[1]
Ac_eq = eq_params(Cg)[2]
println("Nb_eq = ", Nb_eq)
println("Ac_eq = ", Ac_eq)




# -------------------- Generate Plots for Varying Cg --------------------

# Define function for finding root of Eqn () to solve for a critical value of Cg at which Ac = Amax
# x: Cg
function crit_Cg(x)
    repulsion = kb*T*(x + ((3*S^2*x)/(8*Nk*Lk^2)) + ((2*vK*x^2*Nk^2)/S) + ((Qg^2*x^2)/(2*Cm*S)))
    K = KL*exp(-0.5*sp_const*(S-L)*(S-L)/(kb*T))
    zeta = A1*A2*repulsion/(kb*T*K)
    return ((kb*T)/(2*repulsion))*(N1t + N2t - sqrt((N1t-N2t)^2 + 4*zeta)) - Amax
end

Cg_crit = find_zero(crit_Cg, (0, 1e9))
println("Cg_crit = ", Cg_crit)

# Find value of Cg at which Ac = 0
function Cg_zero(x)
    repulsion = kb*T*(x + ((3*S^2*x)/(8*Nk*Lk^2)) + ((2*vK*x^2*Nk^2)/S) + ((Qg^2*x^2)/(2*Cm*S)))
    K = KL*exp(-0.5*sp_const*(S-L)*(S-L)/(kb*T))
    zeta = A1*A2*repulsion/(kb*T*K)
    return ((kb*T)/(2*repulsion))*(N1t + N2t - sqrt((N1t-N2t)^2 + 4*zeta))
end
Cg0 = find_zero(Cg_zero, (0, 1e9))
println("Cg0 = ", Cg0)
# Calculate Nb and S at Cg0
Nb0 = eq_params(Cg0)[1]
Ac0 = eq_params(Cg0)[2]


# Use NLsolve.jl to solve Eqns (7a), (7b), and (7c) to get Nb, S, and Ac
# First, define a function for the system of nonlinear equations
# x[1]: Nb
# x[2]: S
# x[3]: Ac
# U[1]: Eq (7a)
# U[2]: Eq (7b)
# U[3]: Eq (7c)
function u!(U, x)
    K = KL*exp(-0.5*sp_const*(x[2]-L)*(x[2]-L)/(kb*T))
    U[1] = x[1]/x[3] - ((N1t-x[1])/A1)*((N2t-x[1])/A2)*K
    U[2] = x[1]/x[3] - Cg - (3*x[2]^2*Cg)/(8*Nk*Lk^2) - (2*vK*Cg^2*Nk^2)/x[2] - (Qg^2*Cg^2)/(2*Cm*x[2])
    U[3] = x[1]/x[3] - (kb*T)/(sp_const*(L-x[2]))*
        ((3*x[2]*Cg)/(4*Nk*Lk^2) - (2*vK*Cg^2*Nk^2)/(x[2]^2) - (Qg^2*Cg^2)/(2*Cm*x[2]^2))
end
# Also need to define a function for the Jacobian of the system
# x[1]: Nb
# x[2]: S
# x[3]: Ac
# J[1,1]:dU[1]/dNb; J[1,2]:dU[1]/dS; J[1,3]:dU[1]/dAc
# J[2,1]:dU[2]/dNb; J[2,2]:dU[2]/dS; J[2,3]:dU[2]/dAc
# J[3,1]:dU[3]/dNb; J[3,2]:dU[3]/dS; J[3,3]:dU[3]/dAc
function j!(J, x)
    K = KL*exp(-0.5*sp_const*(x[2]-L)*(x[2]-L)/(kb*T))
    n1 = (N1t - x[1])/A1
    n2 = (N2t - x[1])/A2
    g1 = Cg/(Nk*Lk^2)
    g2 = vK*Cg^2*Nk^2
    g3 = (Qg^2*Cg^2)/Cm
    J[1,1] = 1/x[3] - K*((-N1t-N2t+2*x[1])/(A1*A2))
    J[1,2] = -n1*n2*KL/(kb*T)*sp_const*(L-x[2]) *
        exp((-sp_const*(x[2]-L)^2)/(2*kb*T))
    J[1,3] = -x[1]/(x[3]^2)
    J[2,1] = 1/x[3]
    J[2,2] = (-3*x[2]*Cg)/(4*Nk*Lk^2) + (2*vK*Cg^2*Nk^2)/(x[2]^2) + (Qg^2*Cg^2)/(2*Cm*(x[2]^2))
    J[2,3] = -x[1]/(x[3]^2)
    J[3,1] = 1/x[3]
    J[3,2] = -(kb*T)/(sp_const*(L-x[2])^2)*(0.75*x[2]*g1 - 2*g2/(x[2]^2) - g3/(2*(x[2]^2))) - 
        (kb*T)/(sp_const*(L-x[2]))*(0.75*g1 + 4*g2/(x[2]^3) + g3/(x[2]^3))
    J[3,3] = -x[1]/(x[3]^2)
end

# Intialize storage vectors
Cg_vals = Vector{Float64}()
Nb_vals = Vector{Float64}() 
S_vals = Vector{Float64}() 
Ac_vals = Vector{Float64}() 

# Use zero values for Cg, Nb and S as a starting point. Loop through values of Cg in the range (Cg_crit, Cg0)
# and calculate Nb, S, and Ac for each Cg value.
Cg = Cg0
xₒ=[Nb0; S; Ac0]
while Cg > Cg_crit
    global sol = nlsolve(u!, j!, xₒ)
    append!(Cg_vals, Cg)
    append!(Nb_vals, sol.zero[1]/N1t)
    append!(S_vals, (sol.zero[2]-L)/L)
    append!(Ac_vals, sol.zero[3]/Amax)
    global xₒ = sol.zero
    global Cg = Cg - 1
end

# Plot the results
plt1 = plot(Cg_vals, [Nb_vals, S_vals, Ac_vals], xaxis="Cg1 + Cg2 (#/cm^2)", label = ["Nb/N1t" "(S-L)/L" "Ac/Amax"])
savefig(plt1, "./plot1.png")
display(plt1)

# -------------------- Effect of Qg on Dimensionless Variables --------------------
# Here, keep Cg1 and Cg2 constant
Cg1 = 500
Cg2 = 500
Cg = Cg1 + Cg2

# Define function for finding root of Eqn (9) to solve for a critical value of Qg at which Ac = Amax
# x: Qg
function crit_Qg(x)
    repulsion = kb*T*(Cg + ((3*S^2*Cg)/(8*Nk*Lk^2)) + ((2*vK*Cg^2*Nk^2)/S) + ((x^2*Cg^2)/(2*Cm*S)))
    K = KL*exp(-0.5*sp_const*(S-L)*(S-L)/(kb*T))
    zeta = A1*A2*repulsion/(kb*T*K)
    return ((kb*T)/(2*repulsion))*(N1t + N2t - sqrt((N1t-N2t)^2 + 4*zeta)) - Amax
end

Qg_crit = find_zero(crit_Qg, (0, 1e9))
println("Qg_crit = ", Qg_crit)

# Find value of Qg at which Ac = 0
function Qg_zero(x)
    repulsion = kb*T*(Cg + ((3*S^2*Cg)/(8*Nk*Lk^2)) + ((2*vK*Cg^2*Nk^2)/S) + ((x^2*Cg^2)/(2*Cm*S)))
    K = KL*exp(-0.5*sp_const*(S-L)*(S-L)/(kb*T))
    zeta = A1*A2*repulsion/(kb*T*K)
    return ((kb*T)/(2*repulsion))*(N1t + N2t - sqrt((N1t-N2t)^2 + 4*zeta))
end
Qg0 = find_zero(Qg_zero, (0, 1e9))
println("Qg0 = ", Qg0)

# Use NLsolve.jl to solve Eqns (10a), (10b), and (10c) to get Nb and Ac
# Intialize storage vectors
Qg_vals = Vector{Float64}()
Nb_vals = Vector{Float64}() 
S_vals = Vector{Float64}() 
Ac_vals = Vector{Float64}() 

# Use zero values for Qg, Nb and S as a starting point. Loop through values of Qg in the 
# range (Qg_crit, Qg0) and calculate Nb, S, and Ac for each Cg value.
Qg = Qg0
xₒ=[Nb0; S; Ac0]
while Qg > Qg_crit
    global sol = nlsolve(u!, j!, xₒ)
    append!(Qg_vals, Qg)
    append!(Nb_vals, sol.zero[1]/N1t)
    append!(S_vals, (sol.zero[2]-L)/L)
    append!(Ac_vals, sol.zero[3]/Amax)
    global xₒ = sol.zero
    global Qg = Qg - 0.5
end

# Plot the results
plt2 = plot(Qg_vals, [Nb_vals, S_vals, Ac_vals], xaxis="Qg (#/polymer)", label = ["Nb/N1t" "(S-L)/L" "Ac/Amax"])
savefig(plt2, "./plot2.png")
display(plt2)




# -------------------- Effect of N1t, N2t on Range of Allowable Cg values --------------------

# Define function for finding root of Eqn (9) to solve for a critical value of N1t, N2t at which Ac = Amax
# x: N1t, N2t
function crit_N1t(x)
    repulsion = kb*T*(Cg + ((3*S^2*Cg)/(8*Nk*Lk^2)) + ((2*vK*Cg^2*Nk^2)/S) + ((Qg^2*Cg^2)/(2*Cm*S)))
    K = KL*exp(-0.5*sp_const*(S-L)*(S-L)/(kb*T))
    zeta = A1*A2*repulsion/(kb*T*K)
    return ((kb*T)/(2*repulsion))*(x + x - sqrt(4*zeta)) - Amax
end

N1t_crit = find_zero(crit_N1t, (0, 1e9))

# Define function for finding root of Eqn (9) to solve for a critical value of N1t, N2t at which Ac = 0
# x: N1t, N2t
function zero_N1t(x)
    repulsion = kb*T*(Cg + ((3*S^2*Cg)/(8*Nk*Lk^2)) + ((2*vK*Cg^2*Nk^2)/S) + ((Qg^2*Cg^2)/(2*Cm*S)))
    K = KL*exp(-0.5*sp_const*(S-L)*(S-L)/(kb*T))
    zeta = A1*A2*repulsion/(kb*T*K)
    return ((kb*T)/(2*repulsion))*(x + x - sqrt(4*zeta))
end

N1t0 = find_zero(zero_N1t, (0, 1e9))

# Initialize storage vectors
N1t_vals = Vector{Float64}()
Cgcrit_vals = Vector{Float64}() 
Cg0_vals = Vector{Float64}() 

N1t = N1t_crit
N2t = N1t_crit
while N1t > N1t0
    crit = find_zero(crit_Cg, (0, 1e9))
    zero = find_zero(Cg_zero, (0, 1e9))
    append!(Cgcrit_vals, crit)
    append!(Cg0_vals, zero)
    append!(N1t_vals, N1t)
    global N1t = N1t - 1e4
    global N2t = N2t - 1e4
end

plt3 = plot(N1t_vals, Cgcrit_vals, xaxis="N1t = N2t", yaxis = "Cg_crit (#/cm^2)", label = false)
savefig(plt3, "./plot3.png")
display(plt3)

plt4 = plot(N1t_vals, Cg0_vals, xaxis="N1t = N2t", yaxis = "Cg0 (#/cm^2)", label = false)
savefig(plt4, "./plot4.png")
display(plt4)

plt5 = plot(N1t_vals, [Cgcrit_vals, Cg0_vals], xaxis="N1t = N2t (#/cell)", yaxis = "Cg (#/cm^2)", label=["Cg_crit" "Cg0"])
savefig(plt5, "./plot5.png")
display(plt5)