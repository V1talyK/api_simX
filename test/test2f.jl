@info "Тестирование прямого расчёта 2х фазного симулятора"
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
    sim_calc, _ = make_sim2f(grd, gdm_prop, well, prp, nt, satf)

qw = rand(0.:2.:64., nw, nt);
qw[[1,3,7,9],:] .= (-abs.(qw[[1,3,7,9],:]).-10);
qw[[2,4,5,6,8],:] .= abs.(qw[[2,4,5,6,8],:].+10);
qw .= sum(qw,dims=2)./nt

rsl = sim_calc(qw = qw)
wtc = calc_wtc(rsl.SW, gdm_sat.fkrp, well);
wtc[rsl.qw .< 0.0] .= 0.0;
qo = rsl.qw.*(1 .- wtc);  qo[rsl.qw .< 0.0] .= 0.0;
kin = cumsum(sum(qo, dims=1)[:])/(sum(prp.Vp.*(1.0 .- gdm_sat.Sw0))).*gdm_prop.dt
println(kin)
println("wtc=", wtc)
@test (kin[end] <= 1) & all(wtc.<=1)
