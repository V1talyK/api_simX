@info "Тестирование прямого расчёта симулятора, смена управлений на скважинах"
grd, gdm_prop, prp, x, nt = make_gdm(kp_init = 10)

wxy9 = collect(Iterators.product(x,x))[:]
well = make_well(wxy9,grd)
nw = length(well)

sim_calc, cIWC = make_sim(grd,gdm_prop, well, prp, nt)

qw = rand(-1:0.1:1, nw, nt);
qw[[1,3,7,9],:] .= -abs.(qw[[1,3,7,9],:]).-1;
qw[[2,4,5,6,8],:] .= abs.(qw[[2,4,5,6,8],:].+1)

rsl = sim_calc(qw = qw)
AIW, _ = cIWC()

uf = falses(nw, nt);
uf[[2,4,5,6,8],:] .= true;
pw = zeros(nw, nt);
pw[[2,4,5,6,8],:] .= rsl1.pw[[2,4,5,6,8],:];

rsl2 = sim_calc(qw = qw, uf = uf, pw = pw)
@test (sum(abs,rsl1.pw.-rsl2.pw)<1e-3) & (sum(abs,rsl1.qw.-rsl2.qw)<1e-3)

plt = plot(rsl.ppl[2,:])
    plot!(plt, rsl.ppla[2,:])
    plot!(plt, rsl.ppls[2,:])
