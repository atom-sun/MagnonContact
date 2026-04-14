using ITensors, ITensorMPS
using Random
include("util.jl")


# Hamiltonian for the long-range quantum spin chain
function longrange_spinchain_mpo(n_sites, alpha, Jz; boundary_compensate=false)
    """Long-range quantum spin chain.

    H = - ∑ₙₘ J/|m-n|ᵅ(σₙˣσₘˣ + σₙʸσₘʸ) - ∑ₙ Jz σₙᶻσₙ₊₁ᶻ

    In practice, we set J = 1. And Jz denotes the value of relative ratio of 
    the interaction strength over the hopping strength Jz/J.

    This model has the rotation symmetry of Sz, hence conserve the quantum 
    number of magnons. (That is, subspaces with different number of magnons 
    decouples, in this case.)

    To help implementing charge conservation, the Hamiltonian can be expressed 
    in terms of ladder operators σ⁺, σ⁻, instead of the original Pauli matrices.

    H = - ∑ₙₘ 2J/|m-n|ᵅ(σₙ⁻σₘ⁺ + σₙ⁺σₘ⁻) - ∑ₙ Jz σₙᶻσₙ₊₁ᶻ

    Here the ladder operators σ⁺ σ⁻ are introduced as,

        σˣ = σ⁺ + σ⁻,  σʸ = -i(σ⁺ - σ⁻).

    or,
        σ⁻ = (σˣ + iσʸ) / 2,  σ⁺ = (σˣ - iσʸ) / 2.

    In ITensor pkg, σ⁺ and σ⁻ are denoted as S+ and S-, respectively.

    """

    # assume open boundary conditions
    # the model conserves quantum number of magnons
    sites = siteinds("S=1/2", n_sites; conserve_qns=true, conserve_sz=true)

    ampo = AutoMPO()

    # hopping term
    for n in 1:n_sites, m in (n+1):n_sites
        hopp = 1 / (m - n)^alpha
        # To help ITensor identify the conservation, we express the model in 
        # terms of S+,S- operators, instead of the original Pauli matrices.
        ampo .+= -2 * hopp, "S-", n, "S+", m
        ampo .+= -2 * hopp, "S+", n, "S-", m
    end

    # interaction term
    for n in 1:(n_sites-1)
        ampo .+= -Jz, "Z", n, "Z", n + 1
    end

    # TO facilitate dmrg in QN conserved subspaces:
    # compensation for the boundary magnon excitation energy
    # add this term could largely speed up dmrg convergency
    if boundary_compensate
        ampo .+= -Jz, "Z", 1
        ampo .+= -Jz, "Z", n_sites
        ampo .+= Jz, "Id", 1
        ampo .+= Jz, "Id", n_sites
    end

    H = MPO(ampo, sites)
    return H, sites
end


function solve_magnon(
    n_sites, alpha, Jz, num_magnons;
    boundary_compensate=false,
    nsweeps=10,
    bonddims=[10, 20, 100, 200, 400, 800],
    ctf=1E-13,
)

    # Set model.
    H, sites = longrange_spinchain_mpo(
        n_sites, alpha, Jz; boundary_compensate=boundary_compensate)

    # Prepare initial state MPS.
    # Since the model conserves the quantum number of magnons' particle number,
    # we set initial state with specified number of magnons.
    init = ["Up" for _ in 1:n_sites]
    # Here, rand is not available for allowoing duplicated indices.
    # dn = sort(rand(1:n_sites, num_magnons))
    # Use shuffle to avoid duplicated indices.
    dn = shuffle(1:n_sites)[1:num_magnons]
    for n in dn
        init[n] = "Dn"
    end
    psi_i = MPS(sites, init)

    # dmrg parameters
    sweeps = Sweeps(nsweeps)
    setmaxdim!(sweeps, bonddims...)
    setcutoff!(sweeps, ctf)

    # solve dmrg
    energy, psi = dmrg(H, psi_i, sweeps)
    return sites, psi, energy
end


function magnon_creation_symmetric(sites)
    ampo = AutoMPO()
    for n in 1:length(sites)
        ampo += 1.0, "S-", n
    end
    O = MPO(ampo, sites)
    return O
end


function solve_tdvp(
    n_sites, alpha, Jz, num_magnons, t;
    boundary_compensate=false,
    time_step=nothing,
    nsteps=nothing,
    ctf=1E-13,
    maxdim=800,
)
    # get model.
    H, sites = longrange_spinchain_mpo(
        n_sites, alpha, Jz; boundary_compensate=boundary_compensate)

    # initial state: symmetric n-magnon state
    psi0 = productMPS(sites, ["Up" for _ in 1:n_sites])  # vacuum
    createO = magnon_creation_symmetric(sites)
    # create symmetric n-magnon state.
    for _ in 1:num_magnons
        psi0 = apply(createO, psi0; cutoff=ctf, maxdim=maxdim)
    end
    psi0 = normalize!(psi0)
    @info "Initial state prepared: $(flux(psi0))/$(n_sites)sites"

    # define time intervals
    time_step, nsteps = time_step_and_nsteps(t, time_step, nsteps)
    t = time_step * nsteps

    # time evolution
    @info "Start TDVP: N=$n_sites $(num_magnons)magnons t=$t: " *
          " time_step=$time_step nsteps=$nsteps"
    psi_t = tdvp(H, -im * t, psi0;
        time_step=-im * time_step, nsteps=nsteps,
        cutoff=ctf, maxdim=maxdim, normalize=true)
    @info "TDVP evolution done."

    return psi_t, sites, H
end
