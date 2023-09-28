push!(LOAD_PATH, Base.source_dir())

using ITPM_SimX, UnicodePlots

@info "Тестирование прямого расчёта симулятора"
grd, gdm_prop, prp, x, nt = make_gdm()

wxy9 = collect(Iterators.product(x,x))[:]
well = make_well(wxy9,grd)
nw = length(well)

sim_calc = make_sim(grd,gdm_prop, well, prp, nt)

qw = rand(-1:0.1:1, nw, nt);
qw[[1,3,7,9],:] .= -abs.(qw[[1,3,7,9],:]).-1;
qw[[2,4,5,6,7],:] .= abs.(qw[[2,4,5,6,7],:].+1)

rsl = sim_calc(qw = qw)

lineplot(rsl.ppl[1,:])|>println
heatmap(reshape(rsl.PM[:,1], grd.nx, grd.ny))|>println
