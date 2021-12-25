using LinearAlgebra

"""
    numerov_iteration(x::Real, Δ::Real, yj::Real, yi::Real, E::Real; direction::Int64=1)

Calculate the next value of a function y using numerov method.
"""
function numerov_iteration(x::Real, Δ::Real, yj::Number, yi::Number, E::Real, P::Function; direction::Int64=1)
    yk = 2yj*(1 - (5Δ^2 / 12) * P(x, E))
    yk -= yi*(1 + (Δ^2 / 12) * P(x- direction*Δ, E))
    yk /= (1 + (Δ^2 / 12) * P(x + direction*Δ, E))
    return yk
end

"""
    cooley_correction!(x::AbstractArray{Real}, Δ::Real, ψ::AbstractArray{Number}, n_match::Int64, p::Function, v::Function, E; ħ=1, m=1)

Determine the correction in energy eigenvalue using Cooley method.
"""
function cooley_correction!(x::Vector{Float64}, Δ::Real, ψ::Vector{Float64}, n_match::Int64, p::Function, v::Function, E; ħ=1, m=1)
    ΔE  = (1 + (Δ^2 / 12) * p(x[n_match] + Δ, E)) * ψ[n_match + 1]  # Y_{m+1}
    ΔE += (1 + (Δ^2 / 12) * p(x[n_match] - Δ, E)) * ψ[n_match - 1]  # Y_{m-1}
    ΔE -= 2*(1 + (Δ^2 / 12) * p(x[n_match], E)) * ψ[n_match]  # 2*Y_{m}
    ΔE *= (-ħ^2 / (2*m*Δ^2))
    ΔE += ( v(x[n_match]) - E ) * ψ[n_match]
    ΔE *= conj(ψ[n_match]) / sum(abs2.(ψ))
    return ΔE
end

"""
    match_wavefn!(ψ::AbstractArray{Number}, E_n::Real, x::AbstractArray{Real}, Δ::Real, n_match::Int64, P::Function;s=1e-10)

Find the wavefunction for given energy using the Matching method.
"""
function match_wavefn!(
    ψ::Vector{Float64}, 
    E::Real, x::Vector{Float64}, 
    Δ::Real, 
    n_match::Int64,
    P::Function;
    s=1e-10
    )
    Nx = length(x)
    ψ[2] = ψ[end - 1] = s
    for i =  2:n_match - 1
        ψ[i+1] = numerov_iteration(x[i], Δ, ψ[i], ψ[i-1], E, P)
    end

    for i = Nx - 1:-1:n_match + 2
        ψ[i-1] = numerov_iteration(x[i], Δ, ψ[i], ψ[i+1], E, P; direction=-1)
    end

    # Scale the right part equal to left part
    ψ_left = ψ[n_match]
    ψ_right = numerov_iteration(
        x[n_match+1], Δ, 
        ψ[n_match+1], ψ[n_match+2], E, P; 
        direction=-1
    )
    ψ[n_match+1:end] = (ψ_left / ψ_right) * ψ[n_match+1:end]
    return normalize!(ψ)
end

"""
    find_eigenstate(x_min::Real, x_max::Real, x_match::Real, V::Function, E_g::Real; Δ=0.001, m=1.0, ħ=1.0, s=1e-10, err=1e-6)

Find the eigenstate for given potential `V` using the initial guess `E_g`.
"""
function find_eigenstate(
    x_min::Real, 
    x_max::Real, 
    x_match::Real, 
    V::Function,  
    E_min::Real,
    E_max::Real;
    Δ=0.001, m=1.0, ħ=1.0,
    s=1e-10, err=1e-6
    )
    x = collect(x_min:Δ:x_max)
    Nx = length(x)
    n_match = 0
    for i = 1:Nx
        (x[i] <= x_match) ? continue : ((n_match = i); break)
    end

    # Define P(x) for numerov method
    P(x, E) = -2m / ħ^2 * (V(x) - E)

    # Array to store wavefunction
    ψ = zeros(Float64, length(x))

    # Initial Guess for Energy
    E = 0.5*(E_min + E_max)
    match_wavefn!(ψ, E, x, Δ, n_match, P)
    ΔE = cooley_correction!(x, Δ, ψ, n_match, P, V, E; m=m, ħ=ħ)

    while(ΔE >= err)
        E += ΔE
        match_wavefn!(ψ, E, x, Δ, n_match, P)
        ΔE = cooley_correction!(x, Δ, ψ, n_match, P, V, E; m=m, ħ=ħ)
    end

    E_converged = E + ΔE
    match_wavefn!(ψ, E, x, Δ, n_match, P)
    return E_converged, x, ψ
end