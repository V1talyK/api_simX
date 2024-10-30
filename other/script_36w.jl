using Plots, SimScriptTool
grd, gdm_prop, prp, x, nt = make_gdm(;he_init = 40.,
                                     kp_init = 100,
                                     mp_init = 0.2,
                                     nt_init = 120,
                                     nx_init = 51,
                                     ny_init = 51,
                                     Lx_init = 2500,
                                     Ly_init = 2500,
                                     bet = 1e-4,
                                     Paq = 8,
                                     λb = 2.0)

x = make_well_grid(grd, 0.2, 6)
wxy36 = collect(zip(1:36,collect(Iterators.product(x,x))[:]))
well = make_well(wxy36,grd)
nw = length(unique(getindex.(well,2)))

sim_calc = make_sim(grd,gdm_prop, well, prp, nt)
qw = rand(-1:0.1:1, nw, nt);
qw, _uf = SimScriptTool.gen_real_rand_qw(nw, nt; mult = 1.0)
qw .= abs.(qw)
qw[iw_inj,:] .*= -1
qw = SimScriptTool.scaleRateByGeoProp(qw, prp.Vp, prp.kp, prp.he, gdm_prop.dt, gdm_prop.P0, grd.ds)

pw = 2*ones(nw, nt);
uf =  falses(nw, nt)
qw[[9,14,15,23,27],:] .= 0.0
g1 = [1,2,3,7,9,13,14,15]
g2 = [4,5,6,7,10,12,16,17,18]
g3 = [19,20,21,25,27,31,32,33]
g4 = [22,23,24,28,30,34,35,36]
iw_inj = [8,11,26,29]
uf[iw_inj,:] .= true
pw[iw_inj,:] .= 16
qw[8,1:24] .= 0.0; uf[8,1:24] .= false
qw[11,1:36] .= 0.0; uf[11,1:36] .= false
qw[26,1:48] .= 0.0; uf[26,1:48] .= false
qw[29,1:60] .= 0.0; uf[29,1:60] .= false

pw[8,37:end] .=15
pw[8,49:end] .=14
pw[8,61:end] .=13

pw[11,37:end] .=15
pw[11,49:end] .=14
pw[11,61:end] .=10

pw[26,37:end] .=15
pw[26,49:end] .=14
pw[26,61:end] .=13

pw[29,37:end] .=15
pw[29,49:end] .=14
pw[29,61:end] .=16

pw[8,72:84] .= 13 .+ (1:13)./5
pw[11,84:96] .= 9
pw[26,72:107] .= 13 .+sqrt.(1:36)./2

qw[1,72:end] .= 0.0
qw[g2,1:16] .= 0.0
qw[g3,1:36] .= 0.0

rsl = sim_calc(qw = qw, uf = uf, pw = pw)
sum(rsl.qw[rsl.qw.>0])*30.4./(2500*2500*40*0.2)
iw = 4
  plt = plot(rsl.ppl[iw,:])
  plot!(plt, rsl.pw[iw,:])

plot(rsl.qw[iw,:])
qp = sum.(filter.(x->x>0, eachcol(rsl.qw)))
qi = sum.(filter.(x->x<0, eachcol(rsl.qw)))
plt = plot(qp)
  plot(plt, -qi)
using StatsBase
plot(mean(rsl.ppl, dims =1)[:])

plt = plot()
  plt2 = plot()
for (k,v) in enumerate([g1, g2, g3, g4])
  qp = sum.(filter.(x->x>0, eachcol(rsl.qw[v,:])))
  qi = sum.(filter.(x->x<0, eachcol(rsl.qw[[iw_inj[k]],:])))
  pplG = mean.(filter.(x->x>0, eachcol(rsl.ppl[v,:])))
  plot!(plt2, pplG)
  plot!(plt, qp)
  plot!(plt, -qi)
end
display(plt)
display(plt2)

    #
using Dates

pth = "/home/lik/rgm_v1"
tlb = []
vd = range(Date("2000-01-01"), step = Dates.Month(1), length=nt)
for iw = 1:nw
  for t = 1:nt
    prod = rsl.qw[iw,t]>0 ? rsl.qw[iw,t] : 0.0
    inj = rsl.qw[iw,t]<0 ? -rsl.qw[iw,t] : 0.0
    push!(tlb,[iw, vd[t], rsl.pw[iw,t]*10,prod, inj])
  end
end
write_to_csv(pth, tlb, ["well", "date", "pw, атм.", "prod, м3/сут.", "inj, м3/сут."])
