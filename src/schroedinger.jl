module Schroedinger
export find_eigenstate

using Numerov

# ψ is a array to store the values inplace
function match_wavefn!(ψ, E_n, V_x, Nx, Δx, n_match; m=1.0, ħ=1.0, σ=1e-10)
    P(i) = -2m / ħ^2 * (V_x[i] - E_n)
    R(i) = 0
    ψ[2] = ψ[end - 1] = σ
    for i =  2:n_match - 1
        ψ[i+1] = numerov_iteration(Δx, ψ[i], ψ[i-1], P(i+1), P(i), P(i-1), R(i+1), R(i), R(i-1))
    end

    for i = Nx - 1:-1:n_match + 2
        ψ[i-1] = numerov_iteration(Δx, ψ[i], ψ[i+1], P(i-1), P(i), P(i+1), R(i-1), R(i), R(i+1))
    end

    # Scale the right part equal to left part
    ψ_left = ψ[n_match]
    ψ_right = 2ψ[n_match+1] * ((m * Δx^2) * (V_x[n_match+1] - E_n) / ħ^2 + 1) - ψ[n_match+2]
    ψ[n_match+1:end] = (ψ_left / ψ_right) * ψ[n_match+1:end]
    normalize!(ψ)
end

function find_eigenstate(xgrid, V_x, E_min, E_max, x_match; m=1.0, ħ=1.0, σ=1e-10, err=1e-10)
    Nx = length(xgrid)
    Δx = abs(xgrid[2] - xgrid[1])
    n_match = 0
    for i = 1:Nx
        (xgrid[i] <= x_match) ? continue : ((n_match = i); break)
    end

    ψ = zeros(Float64, length(xgrid))
    # TODO: Cooley Correction
    # TODO: Iteration

    # match_wavefn!(ψ, E_converged, V_x, Nx, Δx, n_match; m=m, ħ=ħ, σ=σ)
    # return E_converged, ψ
end

end