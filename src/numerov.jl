module Numerov

export  numerov_iteration

function numerov_iteration(Δx, yj, yi, Pk, Pj, Pi, Rk, Rj, Ri)
    yk = 2yj*(1 - 5(Δx)^2 / 12 * Pj)
    yk -= yi*(1 + (Δx)^2 / 12 * Pi)
    yk += (Δx)^2 / 12 * (Rk + 10Rj + Ri)
    yk /= (1 + (Δx)^2 / 12 Pk)
    return yk
end

end