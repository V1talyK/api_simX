push!(LOAD_PATH, Base.source_dir())

using ITPM_SimX, UnicodePlots
using Revise

@info "Тестирование прямого расчёта симулятора"
grd, gdm_prop, prp, x, nt = make_gdm(kp_init = 0.5)

wxy9 = collect(Iterators.product(x,x))[:]
well = make_well(wxy9,grd)
nw = length(well)

sim_calc = make_sim(grd,gdm_prop, well, prp, nt)

qw = rand(-1:0.1:1, nw, nt);
qw[[1,3,7,9],:] .= -abs.(qw[[1,3,7,9],:]).-1;
qw[[2,4,5,6,8],:] .= abs.(qw[[2,4,5,6,8],:].+1)

rsl = sim_calc(qw = qw)

lineplot(rsl.ppl[1,:])|>println
heatmap(reshape(rsl.PM[:,1], grd.nx, grd.ny))|>println

uf = falses(nw, nt);
uf[[2,4,5,6,8],:] .= true;
pw = zeros(nw, nt);
#pw[[2,4,5,6,8],:] .= 5f0;
pw[[2,4,5,6,8],:] .= rsl.pw[[2,4,5,6,8],:];

rsl = sim_calc(qw = qw, uf = uf, pw = pw)
