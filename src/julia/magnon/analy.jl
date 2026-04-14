# vacuum energy
function vacuum_energy(n_sites, Jz)
    return -Jz * (n_sites - 1)
end

# single-magnon excitation
function ek(k, n_sites, alpha, Jz)
    e = sum([-4 * cos(k * r) / r^alpha for r in 1:n_sites]) + 4 * Jz
    return e
end
