using HDF5

src_path = joinpath(@__DIR__, "..", "src/julia/magnon")
include(joinpath(src_path, "mps.jl"))
include(joinpath(src_path, "data.jl"))


data_path = joinpath(@__DIR__, "..", "data")
filename = joinpath(data_path, "twobody_resonance.txt")
params = process_resonance_parameters(filename)


# two-magnon systems universal behavior
# near two-body resonance
n_sites = 200
num_magnons = 2
@time for (alpha, Jz) in params
    sites, psi, energy =
        solve_magnon(n_sites, alpha, Jz, num_magnons; boundary_compensate=true)
    sz = expect(psi, "Z")
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
