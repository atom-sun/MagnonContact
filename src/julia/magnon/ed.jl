"""Exact Diagonalization of long-range coupling spin model.

"""

using LinearAlgebra
using SparseArrays
include("util.jl")


function hamiltonian_one(n_sites, alpha, Jz; boundary_compensate=false)
    # model Hamiltonian
    hop = [l == m ? 0.0 : -2.0 / abs(l - m)^alpha
           for l in 1:n_sites, m in 1:n_sites]
    ezz0 = -(n_sites - 1) * Jz
    hzz = Vector{Float64}(undef, n_sites)
    hzz[1] = ezz0 + 2 * Jz
    hzz[end] = ezz0 + 2 * Jz
    hzz[2:end-1] .= ezz0 + 4 * Jz
    hzz = Diagonal(hzz)
    ham = hop + hzz

    # boundary magnon excitation energy compensation
    if boundary_compensate
        comp = zeros(n_sites)
        comp[1] = 2 * Jz
        comp[end] = 2 * Jz
        ham += Diagonal(comp)
    end
    return ham
end


function operator_z_one(n_sites)
    # op_z(j) = Diagonal([1 - 2 * (l == j) for l in 1:n_sites])
    # op_z(j) = Diagonal([l == j ? -1.0 : 1.0 for l in 1:n_sites])
    # pre-allocated memory the fastest
    function op_z(j)
        # Zj operator
        v = ones(n_sites)
        v[j] = -1.0
        return Diagonal(v)
    end
    return op_z
end


function operator_zz_one(n_sites)
    function op_zz(i, j)
        # ZiZj operator
        v = ones(n_sites)
        if i != j
            v[i] = -1.0
            v[j] = -1.0
        end
        return Diagonal(v)
    end
    return op_zz
end


function operator_pm_one(n_sites)
    function op_pm(i, j)
        # S+S- operator
        if i == j
            v = ones(n_sites)
            v[i] = 0
            return Diagonal(v)
        else
            return sparse([j], [i], [1.0], n_sites, n_sites)
        end
    end
    return op_pm
end


function ed_one_magnon(n_sites, alpha, Jz; boundary_compensate=false)
    """Exact diagonalization in one-magnon subspace."""

    # get model Hamiltonian
    ham = hamiltonian_one(n_sites, alpha, Jz;
        boundary_compensate=boundary_compensate)

    # solve ground state
    ev = eigen(ham)
    e0 = ev.values[1]
    v0 = ev.vectors[:, 1]
    @info "ground state energy:  $e0"

    # observables/correlation functions
    op_z = operator_z_one(n_sites)
    op_zz = operator_zz_one(n_sites)
    op_pm = operator_pm_one(n_sites)

    # expectation values
    z = [v0' * op_z(j) * v0 for j in 1:n_sites]
    @info "total sz:  $(sum(z))"
    zz = [v0' * op_zz(i, j) * v0 for i in 1:n_sites, j in 1:n_sites]
    pm = [v0' * op_pm(i, j) * v0 for i in 1:n_sites, j in 1:n_sites]

    return e0, v0, z, zz, pm

end


function time_evolution_one(
    n_sites, alpha, Jz, t;
    boundary_compensate=false,
    time_step=nothing,
    nsteps=nothing,
)

    # set evolution intervals
    dt, nsteps = time_step_and_nsteps(t, time_step, nsteps)

    # set initial state: symmetric one-magnon state
    psi0 = ones(n_sites)
    psi0 /= norm(psi0)  # norm(psi0, 2) = sqrt(psi0' * psi0)

    # set evolve Hamiltonian
    ham = hamiltonian_one(n_sites, alpha, Jz;
        boundary_compensate=boundary_compensate)

    # pre-allocate observables storage.
    erg_t = zeros(ComplexF64, nsteps)
    sz_t = zeros(ComplexF64, nsteps, n_sites)
    pm_t = zeros(ComplexF64, nsteps, n_sites, n_sites)
    zz_t = zeros(ComplexF64, nsteps, n_sites, n_sites)

    # define operators
    # unitary evolution: exp(-i H dt)
    ut = exp(-im * dt * ham)
    # observables
    op_z = operator_z_one(n_sites)  # Zj: op_z(j)
    op_zz = operator_zz_one(n_sites)  # ZiZj: op_zz(i, j)
    op_pm = operator_pm_one(n_sites)  # S+iS-j: op_pm(i, j)

    # time evolution
    @info "Start unitary evolution: one-magnon N=$n_sites dt=$dt nsteps=$nsteps"
    psi_t = deepcopy(psi0)
    t = 0.0
    for i in 1:nsteps
        psi_t = ut * psi_t
        e = psi_t' * ham * psi_t
        z = [psi_t' * op_z(j) * psi_t for j in 1:n_sites]
        pm = [psi_t' * op_pm(i, j) * psi_t for i in 1:n_sites, j in 1:n_sites]
        zz = [psi_t' * op_zz(i, j) * psi_t for i in 1:n_sites, j in 1:n_sites]
        # save observables
        erg_t[i] = e
        sz_t[i, :] = z
        pm_t[i, :, :] = pm
        zz_t[i, :, :] = zz
        t += dt
        @info "$i step done.  t=$t, e=$e, total sz:$(sum(z))."
    end

    @info "Unitary time evolution done."
    return psi_t, erg_t, sz_t, pm_t, zz_t
end


function two_magnon_basis(n_sites)
    # define Hilbert space basis.
    # using in operator construction in two-magnon subspace.
    basis = [(i, j) for i in 1:n_sites-1 for j in i+1:n_sites]
    return basis
end


function hamiltonian_two(n_sites, alpha, Jz; boundary_compensate=false)

    # set Hilbirt space with basis
    basis = two_magnon_basis(n_sites)
    # n_states = Int(n_sites * (n_sites - 1) / 2)
    n_states = length(basis)

    # model Hamiltonian
    # hopping matrix
    hop = zeros(n_states, n_states)
    for (ix1, (k, l)) in enumerate(basis)
        for (ix2, (m, n)) in enumerate(basis)
            if l == m
                hop[ix1, ix2] = -2.0 / (n - k)^alpha
            elseif n == k
                hop[ix1, ix2] = -2.0 / (l - m)^alpha
            elseif (k == m) && (l != n)
                hop[ix1, ix2] = -2.0 / abs(l - n)^alpha
            elseif (k != m) && (l == n)
                hop[ix1, ix2] = -2.0 / abs(k - m)^alpha
            end
        end
    end
    # interaction term
    ezz0 = -(n_sites - 1) * Jz * ones(n_states)
    ez = 2 * Jz * [(l > 1) + 2 * (m - l > 1) + (m < n_sites)
                   for (l, m) in basis]
    hzz = Diagonal(ezz0 + ez)
    # Hamiltonian := hopping + interaction
    ham = hop + hzz

    # boundary magnon excitation energy compensation
    if boundary_compensate
        ec = 2 * Jz * [(l == 1) + (m == n_sites) for (l, m) in basis]
        ham += Diagonal(ec)
    end
    return ham
end


function operator_z_two(n_sites)
    # defined Hilbirt space with basis
    basis = two_magnon_basis(n_sites)

    # Zj operator
    op_z(j) = Diagonal([1.0 - 2.0 * (j == l || j == m) for (l, m) in basis])
    return op_z
end


function operator_zz_two(n_sites)
    # defined Hilbirt space with basis
    basis = two_magnon_basis(n_sites)

    # ZiZj operator
    op_zz(i, j) = Diagonal([
        1.0 - 2.0 * (((l == i) + (l == j) + (m == i) + (m == j)) % 2)
        for (l, m) in basis
    ])
    return op_zz
end


function operator_pm_two(n_sites)
    # defined Hilbirt space with basis
    basis = two_magnon_basis(n_sites)

    # SiSj operator
    function op_pm(i, j)
        if i == j
            return Diagonal([1 - (l == j || m == j) for (l, m) in basis])
        else
            v = [(k == j && l != i && m == i && n != j && l == n) ||
                 (k == j && l != i && n == i && m != j && l == m) ||
                 (l == j && k != i && m == i && n != j && k == n) ||
                 (l == j && k != i && n == i && m != j && k == m)
                 for (k, l) in basis, (m, n) in basis]
            return v
        end
    end
end


function ed_two_magnon(n_sites, alpha, Jz; boundary_compensate=false)
    """Exact diagonalization in two-magnon subspace."""

    # get model Hamiltonian and Hilbert space basis
    ham = hamiltonian_two(n_sites, alpha, Jz;
        boundary_compensate=boundary_compensate)

    # solve ground state
    ev = eigen(ham)
    e0 = ev.values[1]
    v0 = ev.vectors[:, 1]
    @info "ground state energy:  $e0"

    # observables/correlation functions
    op_z = operator_z_two(n_sites)
    op_zz = operator_zz_two(n_sites)
    op_pm = operator_pm_two(n_sites)

    # expectation values
    z = [v0' * op_z(j) * v0 for j in 1:n_sites]
    @info "total sz:  $(sum(z))"
    zz = [v0' * op_zz(i, j) * v0 for i in 1:n_sites, j in 1:n_sites]
    pm = [v0' * op_pm(i, j) * v0 for i in 1:n_sites, j in 1:n_sites]

    return e0, v0, z, zz, pm

end


function time_evolution_two(
    n_sites, alpha, Jz, t;
    boundary_compensate=false,
    time_step=nothing,
    nsteps=nothing,
)

    # set evolution intervals
    dt, nsteps = time_step_and_nsteps(t, time_step, nsteps)

    # set initial state: symmetric two-magnon state
    n_states = Int(n_sites * (n_sites - 1) / 2)
    psi0 = ones(n_states)
    psi0 /= norm(psi0)  # norm(psi0, 2) = sqrt(psi0' * psi0)

    # set evolve Hamiltonian
    ham = hamiltonian_two(n_sites, alpha, Jz;
        boundary_compensate=boundary_compensate)

    # pre-allocate observables storage.
    erg_t = zeros(ComplexF64, nsteps)
    sz_t = zeros(ComplexF64, nsteps, n_sites)
    pm_t = zeros(ComplexF64, nsteps, n_sites, n_sites)
    zz_t = zeros(ComplexF64, nsteps, n_sites, n_sites)

    # define operators
    # unitary evolution: exp(-i H dt)
    ut = exp(-im * dt * ham)
    # observables
    op_z = operator_z_two(n_sites)  # Zj: op_z(j)
    op_zz = operator_zz_two(n_sites)  # ZiZj: op_zz(i, j)
    op_pm = operator_pm_two(n_sites)  # S+iS-j: op_pm(i, j)

    # time evolution
    @info "Start unitary evolution: two-magnon N=$n_sites dt=$dt nsteps=$nsteps"
    psi_t = deepcopy(psi0)
    t = 0.0
    for i in 1:nsteps
        psi_t = ut * psi_t
        e = psi_t' * ham * psi_t
        z = [psi_t' * op_z(j) * psi_t for j in 1:n_sites]
        pm = [psi_t' * op_pm(i, j) * psi_t for i in 1:n_sites, j in 1:n_sites]
        zz = [psi_t' * op_zz(i, j) * psi_t for i in 1:n_sites, j in 1:n_sites]
        # save observables
        erg_t[i] = e
        sz_t[i, :] = z
        pm_t[i, :, :] = pm
        zz_t[i, :, :] = zz
        t += dt
        @info "$i step done.  t=$t, e=$e, total sz:$(sum(z))."
    end

    @info "Unitary time evolution done."
    return psi_t, erg_t, sz_t, pm_t, zz_t
end
