using HDF5

src_path = joinpath(@__DIR__, "..", "src/julia/magnon")
include(joinpath(src_path, "ed.jl"))


# demo: one-magnon ED and two-magnon ED with small system.
n_sites = 20
alpha = 1.9
Jz = 1.9

# run one-magnon ED
@time e0, v0, z, zz, pm =
    ed_one_magnon(n_sites, alpha, Jz; boundary_compensate=true);
println("ground state energy: ", e0)
println("total sz: ", sum(z))
@show zz[1:5, 1:5]
@show pm[1:5, 1:5]

# run two-magnon ED
@time e0, v0, z, zz, pm =
    ed_two_magnon(n_sites, alpha, Jz; boundary_compensate=true);
println("ground state energy: ", e0)
println("total sz: ", sum(z))
@show zz[1:5, 1:5]
@show pm[1:5, 1:5]


# computing large system
# at two-body resonance
n_sites = 100
alpha = 1.9
Jz = 1.9

@time e0, v0, z, zz, pm =
    ed_two_magnon(n_sites, alpha, Jz; boundary_compensate=true);
println("ground state energy: ", e0)
println("total sz: ", sum(z))
@show zz[1:5, 1:5]
@show pm[1:5, 1:5]
h5open("run/data/twomagnon-n$n_sites-a$alpha-Jz$Jz-bc-ED.h5", "w") do f
    write(f, "zz", zz)
    write(f, "pm", pm)
end


# time evolution computing
# one-step unitary evolution
# small system benchmark
n_sites = 20
alpha = 1.9
Jz = 1.9
t = 1.0
# run one-magnon
@time psi_t, erg_t, sz_t, pm_t, zz_t =
    time_evolution_one(n_sites, alpha, Jz, t; boundary_compensate=true);
erg = erg_t[end];
sz = sz_t[end, :];
pm = pm_t[end, :, :];
zz = zz_t[end, :, :];
@show erg
@show isreal(sz)
@show isreal(zz)
@show sum(abs.(pm' - pm))
@show findmax(abs.(imag(pm)))
# run two-magnon
@time psi_t, erg_t, sz_t, pm_t, zz_t =
    time_evolution_two(n_sites, alpha, Jz, t; boundary_compensate=true);
erg = erg_t[end];
sz = sz_t[end, :];
pm = pm_t[end, :, :];
zz = zz_t[end, :, :];
@show erg
@show isreal(sz)
@show isreal(zz)
@show sum(abs.(pm' - pm))
@show findmax(abs.(imag(pm)))


# compute large systems
# parallel run
using Distributed
addprocs(5)
@show nprocs()
@everywhere begin
    using HDF5
    using ProgressMeter
    src_path = joinpath(@__DIR__, "..", "src/julia/magnon")
    include(joinpath(src_path, "ed.jl"))

    # cost should scale as (n^2)^2 (or other n's power law larger than that).
    # 100sites around 450s~640s.
    # 150sites around 7200s.
    # test scale as: 7200/630 ~ 1.5^6, n^6 scaling, i.e. (n^3)^2.
    # 200sites run out of memory.
    n_sites = 150
    alpha = 1.9
    Jz = 1.9
end

ts = (1.0, 2.0, 3.0, 4.0, 5.0)
@sync begin
    @showprogress @distributed for t in ts
        @info "Unitary evolution t=$t task starts."
        @time psi_t, erg_t, sz_t, pm_t, zz_t =
            time_evolution_two(n_sites, alpha, Jz, t; boundary_compensate=true)
        erg = erg_t[end]
        sz = sz_t[end, :]
        pm = pm_t[end, :, :]
        zz = zz_t[end, :, :]

        # save full observables data
        cfg_str = "n$n_sites-a$alpha-Jz$Jz-magn2-bc-t$t"
        fname = "run/data/tunitary-$cfg_str-full.h5"
        h5open(fname, "w") do f
            write(f, "erg", erg)
            write(f, "sz", sz)
            write(f, "pm", pm)
            write(f, "zz", zz)
        end

        # * erg, sz, zz are real. only save real part.
        # * pm has both real part and imaginary part.
        #   save pm_real, pm_imag separately.
        erg = real(erg)
        sz = real(sz)
        pm_real = real(pm)
        pm_imag = imag(pm)
        zz = real(zz)

        fname = "run/data/tunitary-$cfg_str.h5"
        h5open(fname, "w") do f
            write(f, "erg", erg)
            write(f, "sz", sz)
            write(f, "pm_real", pm_real)
            write(f, "pm_imag", pm_imag)
            write(f, "zz", zz)
        end
        @info "Unitary evolution t=$t task done."
    end
end
