using Plots
include("numerov_cooley.jl")

function V_harmonic(x)
    return 0.5*x*x
end

# Parameters
x_min = -20.0
x_max = 20.0
x_match = 2.0
E_min = 90
E_max = 95
E_con, x_grid, ψ = find_eigenstate(x_min, x_max, x_match, V_harmonic, E_min, E_max)

# Plots

plot(x_grid, abs2.(ψ), label="E = $(round(E_con;digits=2))")
title!("Eigenfunction for harmonic oscillator potential.")
xlabel!("x")
ylabel!("|ψ(x)|²")
savefig("results/harmonic_psi_2_E_$(round(E_con;digits=2)).png")