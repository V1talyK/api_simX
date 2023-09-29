push!(LOAD_PATH, Base.source_dir())

using ITPM_SimX, UnicodePlots
using NearestNeighbors, UnicodePlots, LinearAlgebra, SuiteSparse, SparseArrays, Statistics

@info "Тестирование прямого расчёта симулятора"
grd, gdm_prop, prp, x, nt = make_gdm(kp_init = 0.5, nt_init = 36000)
gdm_sat = make_gdm_prop_sat()

wxy9 = collect(Iterators.product(x,x))[:]
well = make_well(wxy9,grd)
nw = length(well)

satf = calc_sat_step(prp, grd, gdm_prop, gdm_sat, well, nt)
    sim_calc = make_sim2f(grd,gdm_prop, well, prp, nt, satf)

qw = rand(-1:0.1:1, nw, nt);
qw[[1,3,7,9],:] .= (-abs.(qw[[1,3,7,9],:]).-1)*0;
qw[[2,4,5,6,8],:] .= abs.(qw[[2,4,5,6,8],:].+1)*0
qw[5,:] .= 0.5
qw .= mean(qw,dims=2)

rsl = sim_calc(qw = qw)
wtc = calc_wtc(rsl.SW, gdm_sat.fkrp, well);
wtc[rsl.qw .< 0.0] .= 0.0;
qo = rsl.qw.*(1 .- wtc);  qo[rsl.qw .< 0.0] .= 0.0;
kin = cumsum(sum(qo, dims=1)[:])/(sum(prp.Vp.*(1.0 .- gdm_sat.Sw0))).*gdm_prop.dt
ppo = (cumsum(sum(qw, dims=1)[:]).*gdm_prop.dt)./sum(prp.Vp)

lineplot(rsl.ppl[1,:])|>println
lineplot(wtc[2,:])|>println
lineplot(ppo, kin)|>println
heatmap(reshape(rsl.PM[:,1], grd.nx, grd.ny))|>println
heatmap(reshape(rsl.SW[:,10000], grd.nx, grd.ny))|>println
