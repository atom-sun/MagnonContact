using Test


src_path = joinpath(@__DIR__, "..", "src/julia/magnon")
include(joinpath(src_path, "mps.jl"))
include(joinpath(src_path, "ed.jl"))


@testset "GroundState" begin

    using ITensors
    X = array([[0 1]; [1 0]])
    Y = array([[0 -im]; [im 0]])
    Z = array([[1 0]; [0 -1]])
    Sx = X / 2
    Sy = Y / 2
    Sz = Z / 2
    Sp = array([[0 1]; [0 0]])
    Sm = array([[0 0]; [1 0]])

    @testset "spin_basics" begin
        @test Sp == (X + im * Y) / 2
        @test Sm == (X - im * Y) / 2
        @test Sx == (Sp + Sm) / 2
        @test Sy == -im * (Sp - Sm) / 2
    end

    @testset "ITensor_consistency" begin
        site = siteinds("S=1/2", 1)[1]

        @test array(op("X", site)) == X
        @test array(op("Y", site)) == Y
        @test array(op("Z", site)) == Z
        @test array(op("Sx", site)) == Sx
        @test array(op("Sy", site)) == Sy
        @test array(op("Sz", site)) == Sz
        @test array(op("S+", site)) == Sp
        @test array(op("S-", site)) == Sm
    end

    n_sites = 20
    alpha = 1.75
    Jz = 2.5

    @testset "polarized_state" begin
        num_magnons = 0
        sites, psi, energy = solve_magnon(n_sites, alpha, Jz, num_magnons)
        @test energy == -Jz * (n_sites - 1)

        sz_theo = n_sites - 2 * num_magnons
        @test flux(psi) == QN("Sz", sz_theo)

        sz = expect(psi, "Z")
        totsz = sum(sz)
        @test abs(totsz - sz_theo) < 1E-10

    end

    @testset "single_magnon" begin

        # MPS solving
        num_magnons = 1
        sites, psi, energy = solve_magnon(n_sites, alpha, Jz, num_magnons)
        corr_zz = correlation_matrix(psi, "Z", "Z")
        corr_pm = correlation_matrix(psi, "S+", "S-")

        # test quantum number conservation
        sz_theo = n_sites - 2 * num_magnons
        @test flux(psi) == QN("Sz", sz_theo)

        # test total Sz conservation
        sz = expect(psi, "Z")
        totsz = sum(sz)
        @test abs(totsz - sz_theo) < 1E-10

        # ED solving
        e0, v0, z, zz, pm = ed_one_magnon(n_sites, alpha, Jz)
        @test abs(sum(z) - sz_theo) < 1E-10

        # compare MPS with ED
        @test abs(e0 - energy) < 1E-8
        @test findmax(abs.(zz[1:5, 1:5] - corr_zz[1:5, 1:5]))[1] < 1E-8
        @test findmax(abs.(pm[1:5, 1:5] - corr_pm[1:5, 1:5]))[1] < 1E-8

    end

    @testset "two_magnon" begin
        # compensate boundary magnon excitation energy
        # could largely speed up MPS convergency

        # MPS solving
        num_magnons = 2
        sites, psi, energy = solve_magnon(
            n_sites, alpha, Jz, num_magnons; boundary_compensate=true)
        corr_zz = correlation_matrix(psi, "Z", "Z")
        corr_pm = correlation_matrix(psi, "S+", "S-")

        sz_theo = n_sites - 2 * num_magnons
        @test flux(psi) == QN("Sz", sz_theo)

        sz = expect(psi, "Z")
        totsz = sum(sz)
        @test abs(totsz - sz_theo) < 1E-10

        # ED solving
        e0, v0, z, zz, pm = ed_two_magnon(
            n_sites, alpha, Jz; boundary_compensate=true)
        @test abs(sum(z) - sz_theo) < 1E-10

        # compare MPS with ED
        @test abs(e0 - energy) < 1E-10
        @test findmax(abs.(zz[1:5, 1:5] - corr_zz[1:5, 1:5]))[1] < 1E-9
        @test findmax(abs.(pm[1:5, 1:5] - corr_pm[1:5, 1:5]))[1] < 1E-9

    end

    # TODO: test energy: two-body > single-body > vacuum

end


@testset "TimeEvolution" begin

    # model parameters and time-evolution steps
    n_sites = 20
    alpha = 1.9
    Jz = 1.9
    num_magnons = 2
    t = 1.0
    time_step = 0.05
    nsteps = 20
    t == time_step * nsteps

    @testset "tdvp" begin
        @time psi_t, sites, H =
            solve_tdvp(n_sites, alpha, Jz, num_magnons, t;
                boundary_compensate=true,
                time_step=time_step, nsteps=nsteps,
                ctf=1E-13, maxdim=800)

        # test energy real
        erg = inner(psi_t', H, psi_t)
        @test abs(imag(erg)) < 1E-15
        @show real(erg)

        # test num_magnons conserved
        sz_theo = n_sites - 2 * num_magnons
        @test flux(psi_t) == QN("Sz", sz_theo)
        sz = expect(psi_t, "Z")
        totsz = sum(sz)
        @test abs(totsz - sz_theo) < 1E-14
        sz = expect(psi_t, "Z")
        @test isreal(sz)

        # test correlation functions
        pm = correlation_matrix(psi_t, "S+", "S-")
        zz = correlation_matrix(psi_t, "Z", "Z")

        # pm is hermittian, but not necessarily to be real.
        # * test hermittian:
        # @test ishermitian(pm)  # not work
        @test sum(abs.(pm' - pm)) < 1E-15
        # * pm should not have to be symmetric.
        @test findmax(abs.(transpose(pm) - pm))[1] > 0.001
        @show findmax(abs.(imag.(pm)))

        # zz is hermittian and symmetric, so it's real.
        # physically, zz is equivalent to particle number correlations.
        # * test hermittian:
        # @test ishermitian(zz)  # not work
        @test sum(abs.(zz' - zz)) < 1E-15
        # * test symmetric:
        # @test issymmetric(zz)  # not work
        @test sum(abs.(transpose(zz) - zz)) < 1E-14
        @test all(abs.(imag.(zz)) .< 1E-15)

        @show pm[1:3, 1:3]
        @show zz[1:3, 1:3]
    end

    @testset "ED_unitary_time_evolution" begin
        for (tefunc, num_magnons) in
            zip((time_evolution_one, time_evolution_two), (1, 2))
            @time psi_t, erg_t, sz_t, pm_t, zz_t =
                tefunc(n_sites, alpha, Jz, t; boundary_compensate=true)

            # test energy real
            @test abs(imag(erg_t[end])) < 1E-15

            # test Sz real
            @test isreal(sz_t[end, :])

            # test magnon number conservation
            sz_theo = n_sites - 2 * num_magnons
            totsz = sum(sz_t[end, :])
            @test abs(totsz - sz_theo) < 1E-13

            # test correlation functions
            pm = pm_t[end, :, :]
            zz = zz_t[end, :, :]

            # pm is hermittian, but not necessarily to be real.
            # * test hermittian:
            # @test ishermitian(pm)  # not work
            # * pm should not have to be symmetric.
            @test findmax(abs.(transpose(pm) - pm))[1] > 0.001
            @show findmax(abs.(imag.(pm)))

            # zz is hermittian and symmetric, so it's real.
            # physically, zz is equivalent to particle number correlations.
            # * test hermittian:
            @test ishermitian(zz)  # zz' == zz
            # * test symmetric:
            @test issymmetric(zz)  # transpose(zz) == zz
            # since here is unitary evolution, zz is rigorously real.
            @test isreal(zz)

            @show pm[1:3, 1:3]
            @show zz[1:3, 1:3]
        end
    end
end
