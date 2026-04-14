using HDF5

src_path = joinpath(@__DIR__, "..", "src/julia/magnon")
include(joinpath(src_path, "mps.jl"))


# demo: small system
n_sites = 20
alpha = 1.75
Jz = 2.5

# run magnon MPS
for num_magnons in (0, 1, 2)

    @time sites, psi, energy =
        solve_magnon(n_sites, alpha, Jz, num_magnons; boundary_compensate=true)
    @show flux(psi)
    println(energy)

    sz = expect(psi, "Z")
    totsz = sum(sz)
    println(totsz)

    corr_zz = correlation_matrix(psi, "Z", "Z")
    corr_pm = correlation_matrix(psi, "S+", "S-")
    @show corr_zz[1:5, 1:5]
    @show corr_pm[1:5, 1:5]
    # cannot compute XX or YY correlation, not conservation.
    # correlation_matrix(psi, "X", "X")
    # correlation_matrix(psi, "Y", "Y")

end


# computing large system
n_sites = 100
alpha = 1.75
Jz = 2.45
num_magnons = 2
@time sites, psi, energy =
    solve_magnon(n_sites, alpha, Jz, num_magnons; boundary_compensate=true);
@show flux(psi)
println(energy)

sz = expect(psi, "Z");
totsz = sum(sz);
println(totsz)

corr_zz = correlation_matrix(psi, "Z", "Z");
corr_pm = correlation_matrix(psi, "S+", "S-");

h5open("run/data/twomagnon-n$n_sites-a$alpha-Jz$Jz-bc-MPS.h5", "w") do f
    write(f, "zz", corr_zz)
    write(f, "pm", corr_pm)
    write(f, "sz", sz)
end


# few-body to many-body systems
# near two-body resonance
# universal behavior
n_sites = 200
alpha = 1.9
Jz = 1.9
@time for num_magnons in (2, 3, 4)
    @time sites, psi, energy =
        solve_magnon(n_sites, alpha, Jz, num_magnons; boundary_compensate=true)
    @show flux(psi)
    println(energy)
    sz = expect(psi, "Z")
    totsz = sum(sz)
    println(totsz)
    corr_zz = correlation_matrix(psi, "Z", "Z")
    corr_pm = correlation_matrix(psi, "S+", "S-")
    cfg_str = "num_magnons$num_magnons-n$n_sites-a$alpha-Jz$Jz-bc"
    fname = "run/data/magnon-$cfg_str-MPS.h5"
    h5open(fname, "w") do f
        write(f, "zz", corr_zz)
        write(f, "pm", corr_pm)
        write(f, "sz", sz)
    end
end


# time evolution MPS using TDVP method.
# small system benchmark
# experiment results: 20-step much more correct than 1-step compared with ED!!!
n_sites = 20
alpha = 1.9
Jz = 1.9
# test num_magnons 1 or 2 to compare with ED unitary evolution
num_magnons = 2
# 1 step evolution
t = 1.0
# nsteps = 1
@time psi_t, sites, H =
    solve_tdvp(n_sites, alpha, Jz, num_magnons, t;
        boundary_compensate=true,
        # time_step=time_step, nsteps=nsteps,
        ctf=1E-13, maxdim=800);
@show flux(psi_t)
# energy
erg = inner(psi_t', H, psi_t)
@show erg
@show abs(imag(erg))

# other observables
sz = expect(psi_t, "Z");
@show isreal(sz)
@show sum(sz)

pm = correlation_matrix(psi_t, "S+", "S-");
@show sum(abs.(pm' - pm))
@show findmax(abs.(transpose(pm) - pm))
@show findmax(abs.(imag(pm)))
@show pm[1:3, 1:3]

zz = correlation_matrix(psi_t, "Z", "Z");
@show sum(abs.(zz' - zz))
@show sum(abs.(transpose(zz) - zz))
@show findmax(abs.(imag(zz)))
@show zz[1:3, 1:3]


# 20-steps evolution
t = 1.0
time_step = 0.05
nsteps = 20
t == time_step * nsteps
@time psi_t, sites, H =
    solve_tdvp(n_sites, alpha, Jz, num_magnons, t;
        boundary_compensate=true,
        time_step=time_step, nsteps=nsteps,
        ctf=1E-13, maxdim=800);
@show flux(psi_t)
# energy
erg = inner(psi_t', H, psi_t);
@show erg
@show abs(imag(erg))

# other observables
sz = expect(psi_t, "Z");
@show isreal(sz)
@show sum(sz)

pm = correlation_matrix(psi_t, "S+", "S-");
@show sum(abs.(pm' - pm))
@show findmax(abs.(transpose(pm) - pm))
@show findmax(abs.(imag(pm)))
@show pm[1:3, 1:3]

zz = correlation_matrix(psi_t, "Z", "Z");
@show sum(abs.(zz' - zz))
@show sum(abs.(transpose(zz) - zz))
@show findmax(abs.(imag(zz)))
@show zz[1:3, 1:3]


# compute large system long time evolution
n_sites = 200
alpha = 1.9
Jz = 1.9
t = 1.0
time_step = 0.05
nsteps = 20
t == time_step * nsteps
# time cost:
# * t = 1.0
#   ** 2-magn: ~310s
#   ** 3-magn: ~420s
#   ** 4-magn: ~930s
@time for num_magnons in (2, 3, 4,)
    @time psi_t, sites, H =
        solve_tdvp(n_sites, alpha, Jz, num_magnons, t;
            boundary_compensate=true,
            time_step=time_step, nsteps=nsteps,
            ctf=1E-13, maxdim=800)
    @show flux(psi_t)

    # energy
    erg = inner(psi_t', H, psi_t)
    @show erg
    @show abs(imag(erg))

    # other observables
    sz = expect(psi_t, "Z")
    @show isreal(sz)

    pm = correlation_matrix(psi_t, "S+", "S-")
    @show sum(abs.(pm' - pm))
    @show findmax(abs.(transpose(pm) - pm))
    @show findmax(abs.(imag(pm)))

    zz = correlation_matrix(psi_t, "Z", "Z")
    @show sum(abs.(zz' - zz))
    @show sum(abs.(transpose(zz) - zz))
    @show findmax(abs.(imag(zz)))

    # save full observables data
    cfg_str = "n$n_sites-a$alpha-Jz$Jz-magn$num_magnons-bc-t$t" *
              "-time_step$time_step-nsteps$nsteps"
    fname = "run/data/tdvp-$cfg_str-full.h5"
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

    fname = "run/data/tdvp-$cfg_str.h5"
    h5open(fname, "w") do f
        write(f, "erg", erg)
        write(f, "sz", sz)
        write(f, "pm_real", pm_real)
        write(f, "pm_imag", pm_imag)
        write(f, "zz", zz)
    end
end
