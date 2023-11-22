push!(LOAD_PATH, Base.source_dir())

using ITPM_SimX
using NearestNeighbors, UnicodePlots, LinearAlgebra, SuiteSparse, SparseArrays, Statistics
include("libs.jl")
include("funs_2D.jl")
include("funs_sat.jl")
include("linAlgLib.jl")


@info "Тестирование прямого расчёта симулятора"
grd, gdm_prop, prp, x, nt = make_gdm(kp_init = 0.5, he_init = 5, nt_init = 480)

grd, gdm_prop, prp, x, nt = make_gdm(kp_init = 10,
                                     he_init = 0.5,
                                     nt_init = 60,
                                     nx_init = 100,
                                     ny_init = 100,
                                     Lx_init = 1000,
                                     Ly_init = 1000)

gdm_sat = make_gdm_prop_sat(mu_o = 3f0)

wxy9 = collect(Iterators.product(x,x))[:]
well = make_well(wxy9,grd)
nw = length(well)

satf = calc_sat_step(prp, grd, gdm_prop, gdm_sat, well, nt)
    sim_calc = make_sim2f(grd, gdm_prop, well, prp, nt, satf)

qw = rand(0.:2.:64., nw, nt);
#qw[[1,3,7,9],:] .= (-abs.(qw[[1,3,7,9],:]).-10);
#qw[[2,4,5,6,8],:] .= abs.(qw[[2,4,5,6,8],:].+10);
#qw[5,:] .= 72
qw .= mean(qw,dims=2)

rsl = sim_calc(qw = qw)
wtc = calc_wtc(rsl.SW, gdm_sat.fkrp, well);
wtc[rsl.qw .< 0.0] .= 0.0;
qo = rsl.qw.*(1 .- wtc);  qo[rsl.qw .< 0.0] .= 0.0;
kin = cumsum(sum(qo, dims=1)[:])/(sum(prp.Vp.*(1.0 .- gdm_sat.Sw0))).*gdm_prop.dt
dpo = (cumsum(sum(ifelse.(qw.>0,qw,0.), dims=1)[:]).*gdm_prop.dt)./sum(prp.Vp)
ppo = (cumsum(sum(abs, ifelse.(qw.<0,qw,0.), dims=1)[:]).*gdm_prop.dt)./sum(prp.Vp)

lineplot(rsl.ppl[5,:])|>println
lineplot(wtc[5,:])|>println
lineplot(dpo, kin, xlabel= "ППО", ylabel="КИН")|>println
heatmap(reshape(rsl.PM[:,25], grd.nx, grd.ny))|>println
heatmap(reshape(rsl.SW[:,25], grd.nx, grd.ny))|>println

#@profiler sim_calc(qw = qw)
