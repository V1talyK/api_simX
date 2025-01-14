using Plots, SimScriptTool, UnicodePlots
tag = "13w"
grd, gdm_prop, prp, x, nt = make_gdm(;he_init = 10.,
                                     kp_init = 100,
                                     mp_init = 0.2,
                                     nt_init = 240,
                                     nx_init = 51,
                                     ny_init = 51,
                                     Lx_init = 2500,
                                     Ly_init = 2500,
                                     bet = 1e-4,
                                     Paq = 10,
                                     λb = 0.25)

x = ITPM_SimX.make_well_grid(grd, 0.25, 5)
wxy13 = collect(zip(1:13,collect(Iterators.product(x,x))[1:2:25]))
well = make_well(wxy13,grd)
nw = length(unique(getindex.(well,2)))
iw_inj = [4,5,9,10]
iw_prod = setdiff(1:13,iw_inj)


prp.kp[.&(abs.((grd.X[well[2][1]]+grd.X[well[5][1]])/2 .- grd.X).<20,
        (grd.Y[well[5][1]]+grd.Y[well[6][1]])/2 .<=grd.Y,
        (grd.Y[well[12][1]]+grd.Y[well[13][1]])/2 + 600 .>=grd.Y)] .*= 0.01

prp.kp[.&((grd.X[well[2][1]]+grd.X[well[5][1]])/2 .< grd.X,
        abs.((grd.Y[well[5][1]]+grd.Y[well[6][1]])/2 .- grd.Y).<20)] .*= 0.01


prp.kp[.&(3*grd.X.+grd.Y.>3100,
          3*grd.X.+grd.Y.<4100,
          grd.X.>150.0,
          grd.X.<2250.0,
          grd.Y.>150.0,
          grd.Y.<2250.0)] .*= 4

#UnicodePlots.heatmap(reshape(prp.kp, grd.nx, grd.ny)')|>println
wxy = getindex.(wxy13,2)
plt = plot_map_and_well(range(0, 2500, grd.nx),
                  range(0, 2500, grd.ny),
                  reshape(prp.kp, grd.nx, grd.ny)', wxy, [iw_prod, iw_inj], ["Доб.", "Наг."],
                  [:circle, :dtriangle])
Plots.annotate!(plt, getindex.(wxy,1).+120, getindex.(wxy,2).-50, text.(string.("№", 1:nw),12))

pth = joinpath(Base.source_dir(),"map_13w")
Plots.savefig(pth)
Plots.svg(pth)

sim_calc, cIWC = make_sim(grd,gdm_prop, well, prp, nt)
qw = rand(-1:0.1:1, nw, nt);
qw, _uf = SimScriptTool.gen_real_rand_qw(nw, nt; mult = 1.0)
qw .= abs.(qw)
qw[iw_inj,:] .*= -1
qw = SimScriptTool.scaleRateByGeoProp(qw, prp.Vp, prp.kp, prp.he, gdm_prop.dt, gdm_prop.P0, grd.ds)

pw = 2*ones(nw, nt);
uf =  falses(nw, nt)
uf[iw_inj,:] .= true
pw[iw_inj,:] .= 14
uf[iw_inj,1:36] .= false
qw[iw_inj,1:12] .= 0.0
qw[iw_inj[1],13:36] .= -qw[iw_inj[1],13:36]./4
qw[iw_inj[2],19:36] .= -qw[iw_inj[2],19:36]./4; qw[iw_inj[2],13:18] .=0;
qw[iw_inj[3],25:36] .= -qw[iw_inj[3],25:36]./4; qw[iw_inj[3],13:24] .=0;
qw[iw_inj[4],31:36] .= -qw[iw_inj[4],31:36]./4; qw[iw_inj[4],13:30] .=0;

pw[iw_inj,:] .+= rand(-0.2:0.1:0.2,length(iw_inj),nt)
pw[iw_inj[1],25:36] .-=1; pw[iw_inj[1],160:180] .-=4
pw[iw_inj[2],48:60] .+=2; pw[iw_inj[2],100:120] .+=3; pw[iw_inj[2],140:180] .+=3
pw[iw_inj[3],54:80] .-=(54:80)./80; pw[iw_inj[3],140:200] .-=(0:60)./30
pw[iw_inj[4],36:nt] .+=sin.((36:nt)./10)-cos.((36:nt)./12)

qw[iw_prod[4],:] .= 2000 .+qw[iw_prod[4],:]./(maximum(qw[iw_prod[4],:])-minimum(qw[iw_prod[4],:]))*500
qw[iw_prod[4],72:end] .= 0.0
#qw[iw_prod[4],qw[iw_prod[4],:].<1800] .= 1800
#uf[iw_prod[4],:] .= true
#pw[iw_prod[4],:]

qw[iw_prod[7],:] .= 1800 .+qw[iw_prod[7],:]./(maximum(qw[iw_prod[7],:])-minimum(qw[iw_prod[7],:]))*300
#uf[iw_prod[7],:] .= true
#pw[iw_prod[7],:]

flag = true;
k = 0
while flag
  k+=1
  rsl = sim_calc(qw = qw, uf = uf, pw = pw)
  ia = findall(rsl.pw.<0.05)
  flag = (&)(k<10, length(ia)>0)
  if flag
    qw[ia] .*=0.7
  end
  println(length(ia))
end
rsl = sim_calc(qw = qw, uf = uf, pw = pw)

sum(rsl.qw[rsl.qw.>0])*30.4./(2500*2500*10*0.2)
iw =1
  plt = plot(rsl.pw[iw,:], lw = 2, label = "заб.")
  plot!(plt, rsl.ppl[iw,:], lw = 2, label = "яч.")
  plot!(plt, rsl.ppla[iw,:], lw = 2, label = "зона.")
  plot!(plt, rsl.ppls[iw,:], lw = 2, label = "ост.")
  Plots.ylabel!(plt, "Давление, МПа")


plot(rsl.qw[iw,:], st = :step)
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

plt = plot_map_and_well(range(0, 2500, grd.nx),
                  range(0, 2500, grd.ny),
                  reshape(rsl.PM[:,end], grd.nx, grd.ny)', wxy, [iw_prod, iw_inj], ["Доб.", "Наг."],
                  [:circle, :dtriangle])
Plots.annotate!(plt, getindex.(wxy,1).+120, getindex.(wxy,2).-50, text.(string.("№", 1:nw),12))
#
#     #
# using Dates
#
# pth = "/home/lik/rgm_m13_v0"
# tlb = []
# vd = range(Date("2000-01-01"), step = Dates.Month(1), length=nt)
# for iw = 1:nw
#   for t = 1:nt
#     prod = rsl.qw[iw,t]>0 ? rsl.qw[iw,t] : 0.0
#     inj = rsl.qw[iw,t]<0 ? -rsl.qw[iw,t] : 0.0
#     push!(tlb,[iw, vd[t], rsl.pw[iw,t]*10,prod, inj])
#   end
# end
# write_to_csv(pth, tlb, ["well", "date", "pw, атм.", "prod, м3/сут.", "inj, м3/сут."])

tlb = Vector(undef,0)
for iw = 1:nw
    for t = 1:nt
        qp = ifelse(rsl.qw[iw,t] >= 0.0, rsl.qw[iw,t], 0.0)
        qi = ifelse(rsl.qw[iw,t] <= 0.0, rsl.qw[iw,t], 0.0)
        opra_i = qi > 0.0
        opra_p = qp > 0.0
        tmp = Vector{Any}(undef, 0)
        append!(tmp, iw, t, gdm_prop.dt*t, qp, 0.0, qp, qi, rsl.pw[iw,t], rsl.ppl[iw,t], opra_p, opra_i)
        push!(tlb, tmp)
    end
end

head = ["id","numb","date","liquid","oil","water","injection","bhp","cellp_m","opra_p","opra_i"]
pth = joinpath(Base.source_dir(),"script_$(tag)_2500.csv")
    write_to_csv(pth, tlb, head)

prm = Dict("bnd"=>[(0.0, 0.0),
                   (grd.nx*grd.dx, 0.0),
                   (grd.nx*grd.dx, grd.ny*grd.dx),
                   (0.0, grd.ny*grd.dx) ],
           "wxy"=>wxy13,
           "wkp"=>prp.kp[getindex.(well,1)],
           "whe"=>prp.he[getindex.(well,1)],
           "kp"=>round.(prp.mp, sigdigits=3),
           "he"=>round.(prp.he, sigdigits=3),
           "mp"=>round.(prp.mp, sigdigits=3),
           "prop"=>gdm_prop)

pth = joinpath(Base.source_dir() ,"script_$(tag)_2500.json")
    write_to_json(pth, prm)

gdm_sat = make_gdm_prop_sat(mu_o = 3.0f0, n_o = 2, n_w = 2)
satf = calc_sat_step(prp, grd, gdm_prop, gdm_sat, well, nt)
    sim_calc, cIWf = make_sim2f(grd, gdm_prop, well, prp, nt, satf)

qw_base = copy(qw)
qw = copy(qw_base)
flag = true;
k = 0
while flag
  k+=1
  rsl = sim_calc(qw = qw, uf = uf, pw = pw)
  satf.Sw0i.=gdm_sat.Sw0
  ia = findall(rsl.pw.<0.05)
  flag = (&)(k<10, length(ia)>0)
  if flag
    qw[ia] .*=0.7
  end
  println(length(ia))
end
satf.Sw0i.=gdm_sat.Sw0
rsl = sim_calc(qw = qw, uf = uf, pw = pw)
w2w, _ = cIWf(qw=qw)
wtc = calc_wtc(rsl.SW, gdm_sat.fkrp, well);
wtc[rsl.qw .< 0.0] .= 0.0;
qo = rsl.qw.*(1 .- wtc);  qo[rsl.qw .< 0.0] .= 0.0;
iw = 8
  plt = plot(qw[iw,:], label = "жид.")
  plt_twin = twinx(plt);
  plot!(plt, qo[iw,:], label = "нефть.")
  plot!(plt_twin, wtc[iw,:], label = "обв.", lc = :black)
kin = cumsum(sum(qo, dims=1)[:])/(sum(prp.Vp.*(1.0 .- gdm_sat.Sw0))).*gdm_prop.dt
plot(kin)

plot(getindex.(w2w,1,5))

plt = plot_map_and_well(range(0, 2500, grd.nx),
                  range(0, 2500, grd.ny),
                  reshape(rsl.SW[:,240], grd.nx, grd.ny)', wxy, [iw_prod, iw_inj], ["Доб.", "Наг."],
                  [:circle, :dtriangle])
Plots.annotate!(plt, getindex.(wxy,1).+120, getindex.(wxy,2).-50, text.(string.("№", 1:nw),12))


tlb = Vector(undef,0)
for iw = 1:nw
    for t = 1:nt
        qp = ifelse(rsl.qw[iw,t] >= 0.0, rsl.qw[iw,t], 0.0)
        qi = ifelse(rsl.qw[iw,t] <= 0.0, rsl.qw[iw,t], 0.0)
        #qo = ifelse(rsl.qw[iw,t] <= 0.0, rsl.qw[iw,t], 0.0)
        opra_i = qi > 0.0
        opra_p = qp > 0.0
        tmp = Vector{Any}(undef, 0)
        append!(tmp, iw, t, gdm_prop.dt*t, qp, qo[iw, t], qp.-qo[iw, t], qi, rsl.pw[iw,t], rsl.ppl[iw,t], opra_p, opra_i)
        push!(tlb, tmp)
    end
end

head = ["id","numb","date","liquid","oil","water","injection","bhp","cellp_m","opra_p","opra_i"]
pth = joinpath(Base.source_dir(),"script_$(tag)_2500_2f.csv")
    write_to_csv(pth, tlb, head)

prm = Dict("bnd"=>[(0.0, 0.0),
                   (grd.nx*grd.dx, 0.0),
                   (grd.nx*grd.dx, grd.ny*grd.dx),
                   (0.0, grd.ny*grd.dx) ],
           "wxy"=>wxy13,
           "wkp"=>prp.kp[getindex.(well,1)],
           "whe"=>prp.he[getindex.(well,1)],
           "kp"=>round.(prp.mp, sigdigits=3),
           "he"=>round.(prp.he, sigdigits=3),
           "mp"=>round.(prp.mp, sigdigits=3),
           "prop"=>gdm_prop,
           "2f"=>Dict("mu_w"=>1, "mu_o"=>3, "n_w"=>2, "n_o"=>2))

pth = joinpath(Base.source_dir() ,"script_$(tag)_2500_2f.json")
    write_to_json(pth, prm)
