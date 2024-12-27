function make_sim(grd, gdm_prop, well, prp, nt)
    nc = grd.nc
    nw = length(unique(getindex.(well,2)))
    nwc = length(well)
    w1 = getindex.(well,1)
    w2 = getindex.(well,2)

    rc = grd.rc
    λbi = grd.λbi
    he = copy(prp.he)

    Paq = gdm_prop.Paq
    P0 = gdm_prop.P0
    dt = gdm_prop.dt

    λbc = zeros(Float32, nc);
    PM = zeros(Float32, nc+nw,nt)
    PM0 = zeros(Float32, nc+nw)
    pwc = zeros(Float32, nw,nt)
    pplc = zeros(Float32, nw,nt) #Давление в ячейке
    ppla = zeros(Float32, nw,nt) #Давление осреднено по площади
    ppls = zeros(Float32, nw,nt) #Давление в остановленной скважине
    qwc = zeros(Float32, nw,nt)
    qcl = zeros(Float32, nwc)

    GM = ones(Float32, nc);
    bb = zeros(Float32,nc+nw)
    bs = zeros(Float32, nc);
    PMs = zeros(Float32, nc);

    λbc[λbi] .= gdm_prop.λb;
    M2Wa = make_M2W(grd, well)

    dw = inv.(accumarray(w2,ones(length(w2))))[w2]
    tM = (M2M = sparse(w2,w1,dw,nw,nc),
          M2Mw = sparse(w1,w2,1,nc+nw,nw),
          W2Ma = M2Wa')
    actW = trues(nw,nt)
    ufl = falses(nw,nt)

    WI = 2*pi./fill(log(0.14*sqrt(2)*grd.dx/0.05),nwc)
    A, W1 = makeA(view(rc,:,1),view(rc,:,2),nc,nw,w1,w2)
    AS = sparse(view(rc,:,1),view(rc,:,2),1.0,nc,nc)
    makeAG = make_fun_AG(grd.nc,grd.rc,grd.dx,grd.ds);

    qw0 = zeros(Float32, nw, nt)
    pw0 = ones(Float32, nw, nt)
    wc0 = ones(Float32, nw, nt)
    uf0 = falses(nw, nt)
    kp0 = prp.kp
    he0 = prp.he
    flag_stop_well0 = false; #Флаг на моделирование остановки скважин

    function msim(; qw = qw0, pw = pw0, kp = kp0, he = he0, uf = uf0, wc = wc0, fsw = flag_stop_well0)
        GM.=kp.*he.*10. *8.64*1e-3;
        AG, T = makeAG(kp.*10. *8.64*1e-3,he)
        wct = view(wc, w2, 1)
        uft = view(uf, :, 1)
        updA!(A,W1,AG,view(rc,:,1),view(rc,:,2),nc,nw,T,λbc,w1,w2,GM,WI,wct,uft,prp.eVp)
        ACL = cholesky(-A)
        CL = make_CL_in_julia(ACL, Threads.nthreads())
        updateCL!(CL, ACL)

        updAS!(AS,AG,view(rc,:,1),view(rc,:,2),nc,T,λbc,prp.eVp)
        ACLS = cholesky(-AS)
        CLS = make_CL_in_julia(ACLS, Threads.nthreads())
        updateCL!(CLS, ACLS)

        PM0 .= P0
        AM = transpose(convert(Array{Float32,2}, ACL\tM.M2Mw))
        uuf = !all(allunique.(eachrow(uf)))

        for t=1:nt
            if uuf
                uft = view(uf, :, t)
                updA!(A,W1,AG,view(rc,:,1),view(rc,:,2),nc,nw,T,λbc,w1,w2,GM,WI,wct,uft,prp.eVp)
                cholesky!(ACL, -A)
                updateCL!(CL, ACL)
            end

            wct = view(wc, w2, t)
            bs.=0f0;
            bs .= bs .- T.*λbc*Paq .- prp.eVp.*view(PM0, 1:nc)

            PM[:,t], pwc[:,t], pplc[:,t], qwc[:,t] = sim_step!(PM0, qcl, CL, bb,
                            nc,nw,Paq,T,well,
                            view(uf,:,t),view(qw,:,t), view(pw,:,t),
                            λbc, WI, wct, prp.eVp, tM, w1, w2)
            if fsw
                bs[w1] .= view(bs, w1) .+ qcl
                for iw = 1:nw
                    fl = w2.==iw
                    bs[w1[fl]] -= qcl[fl]
                    #PMs .= ACLS\bs
                    back_slash_slvr!(PMs, CLS, bs)
                    ppls[iw,t] = - sum(view(PMs, w1[fl]))/count(fl)
                    bs[w1[fl]]  += qcl[fl]
                end
            end
            ppla[:,t] .= tM.W2Ma*view(PM0,1:nc)
        end
        pwc[uf].=pw[uf]
        rsl = (ppl = pplc, qw = qwc, pw = pwc,  PM = PM[1:nc,:], ppla = ppla,
               ppls = ppls)
        return rsl
    end

    function cIWC(; qw = qw0, pw = pw0, kp = kp0, he = he0, uf = uf0, wc = wc0)
        #Расчёт матрицы связи скважин
        GM.=kp.*he.*10. *8.64*1e-3;
        AG, T = makeAG(kp.*10. *8.64*1e-3,he)
        wct = view(wc, w2, 1)
        uft = view(uf, :, 1)
        updA!(A,W1,AG,view(rc,:,1),view(rc,:,2),nc,nw,T,λbc,w1,w2,GM,WI,wct,uft,prp.eVp)
        ACL = cholesky(-A)
        # CL = make_CL_in_julia(ACL, Threads.nthreads())
        # updateCL!(CL, ACL)
        PM0 .= P0
        AM = transpose(convert(Array{Float32,2}, ACL\tM.M2Mw))
        rAdf, rBdf = make_reduce_ma3x_dims(A, w1, w2, nc)
        AIW = rAdf(A, T, λbc)

        # for t = 1:1
        #     wct = view(wc, w2, t)
        #     PM[:,t], pwc[:,t], pplc[:,t], qwc[:,t] = sim_step!(PM0, ACL, bb,
        #                     nc,nw,Paq,T,well,
        #                     view(uf,:,t),view(qw,:,t), view(pw,:,t),
        #                     λbc, WI, wct, prp.eVp, tM, w1, w2)
        # end

        BIW = 0.0#rBdf(bb)

        return AIW, BIW
    end

    return msim, cIWC
end

function sim_step!(PM0, qcl, ACL, bb, nc,nw,Paq,T,well,uft,qwt,pwt,λbc,
    WI, wct, eVp, tM, w1, w2)

    makeB!(bb, nc,nw,Paq,T,well,uft,qwt,pwt,λbc, view(PM0,1:nc), WI, wct, eVp);

    #PM0 .= ACL\bb;
    back_slash_slvr!(PM0, ACL, bb)
    PM0 .= .-PM0
    pwc = view(PM0,nc+1:nc+nw)
    pplc = tM.M2M*view(PM0,1:nc);
    qcl .= WI.*view(T,w1).*(view(PM0,w1).-view(pwt,w2)).*wct
    qwc = accumarray(w2, qcl, nw)
    nuft = .!uft
    qwc[nuft] .= qwt[nuft]
    pwc[uft] .= pwt[uft]
    qcl .= WI.*view(T,w1).*(view(PM0,w1).-view(pwc,w2)).*wct
    #println(sum(abs,pplcBt.-temp_ppl))
    return PM0, pwc, pplc, qwc
end

function sim_step!(PM0, ACL, bb, nc,nw,Paq,T,well,uft,qwt,pwt,λbc,
    WI, wct, eVp, tM, w1, w2)
    @info "Старый вызов"
    qcl = zeros(Float32, length(w1))
    return sim_step!(PM0, qcl, ACL, bb, nc,nw,Paq,T,well,uft,qwt,pwt,λbc,
        WI, wct, eVp, tM, w1, w2)
end

function make_sim2f(grd, gdm_prop, well, prp, nt, satc)
    nc = grd.nc
    nw = length(unique(getindex.(well,2)))
    nwc = length(well)
    w1 = getindex.(well,1)
    w2 = getindex.(well,2)

    rc = grd.rc
    bnd_ind = grd.λbi
    he = copy(prp.he)

    Paq = gdm_prop.Paq
    P0 = gdm_prop.P0
    dt = gdm_prop.dt

    λbc = zeros(Float32, nc);
    PM = zeros(Float32, nc+nw,nt)
    PM0 = zeros(Float32, nc+nw)
    SW = zeros(Float32, nc,nt)
    pwc = zeros(nw,nt)
    pplc = zeros(nw,nt)
    qwc = zeros(nw,nt)
    qcl = zeros(Float32, nwc)

    GM = ones(Float32, nc);
    Mbt = ones(Float32, nc); Mbt .= satc.fkrp.w.(satc.Sw0i).+satc.fkrp.o.(1f0 .- satc.Sw0i)
    Tp = zeros(Float32, size(rc,1))
    bb = zeros(Float32,nc+nw)

    λbc[bnd_ind] .= gdm_prop.λb;

    dw = inv.(accumarray(w2,ones(length(w2))))[w2]
    tM = (M2M = sparse(w2,w1,dw,nw,nc),
          M2Mw = sparse(w1,w2,1,nc+nw,nw))
    actW = trues(nw,nt)
    ufl = falses(nw,nt)

    WI = 2*pi./fill(log(0.14*sqrt(2)*grd.dx/0.05),nwc)
    A, W1 = makeA(view(rc,:,1),view(rc,:,2),nc,nw,w1,w2)
    makeAG = make_fun_AG(grd.nc,grd.rc,grd.dx,grd.ds);

    qw0 = zeros(Float32, nw, nt)
    pw0 = ones(Float32, nw, nt)
    wc0 = ones(Float32, nw, nt)
    uf0 = falses(nw, nt)
    kp0 = prp.kp
    he0 = prp.he
    stream_flag = trues(length(Tp))
    Tpa = zeros(Float32, length(bnd_ind))
    Tpa1 = zeros(Float32, nc)
    bnd_stream_flag = trues(length(Tpa))

    WTp = ones(Float32, nw)
    AIW = Vector(undef, nt)

    function msim(; qw = qw0, pw = pw0, kp = kp0, he = he0, uf = uf0, wc = wc0)
        GM.=kp.*he.*10. *8.64*1e-3;
        AG, T = makeAG(kp.*10. *8.64*1e-3,he)
        PM0 .= P0
        wct = view(wc, w2, 1)
        uft = view(uf, :, 1)
        updA!(A,W1,AG,view(rc,:,1),view(rc,:,2),nc,nw,T,λbc,w1,w2,GM,WI,wct, uft, prp.eVp)
        ACL = cholesky(-A)
        #AM = transpose(convert(Array{Float32,2}, ACL\tM.M2Mw))
        for t=1:nt
            uft = view(uf, :, t)
            stream_flag .= view(PM0,view(rc,:,1)) .> view(PM0,view(rc,:,2))
            Tp[stream_flag] = Mbt[view(rc,stream_flag,1)]
            stream_flag .= view(PM0,view(rc,:,1)) .< view(PM0,view(rc,:,2))
            Tp[stream_flag] = Mbt[view(rc,stream_flag,2)]
            stream_flag .= view(PM0,view(rc,:,1)) .== view(PM0,view(rc,:,2))
            Tp[stream_flag] = 2*Mbt[view(rc,stream_flag,1)].*Mbt[view(rc,stream_flag,2)]./(Mbt[view(rc,stream_flag,1)].+Mbt[view(rc,stream_flag,2)])

            WTp.=1f0;
            WTp[qw[:,t].>0.] .= Mbt[w1[qw[:,t].>0.]]

            bnd_stream_flag .= view(PM0,bnd_ind) .> Paq
            Tpa[bnd_stream_flag] = view(Mbt[bnd_ind], bnd_stream_flag)
            # println(satc.fkrp.w.(ones(length(count(.!bnd_stream_flag)))))
            # println(Tpa[.!bnd_stream_flag])
            Tpa[.!bnd_stream_flag] = satc.fkrp.w.(ones(Float32, count(.!bnd_stream_flag)))

            wct = view(wc, w2, t)
            updA!(A,W1,AG.*Tp,view(rc,:,1),view(rc,:,2),nc,nw,T,λbc,w1,w2,GM,WI, WTp.*wct,uft,prp.eVp)
            cholesky!(ACL, -A)

            PM[:,t], pwc[:,t], pplc[:,t], qwc[:,t] = sim_step!(PM0, qcl, ACL, bb,
                            nc,nw,Paq,T,well,
                            uft,view(qw,:,t), view(pw,:,t),
                            λbc, WI.*WTp, wct, prp.eVp, tM, w1, w2)

            SW[:,t], Mbt[:] = satc(PM0, view(qwc,:,t), GM, AG, false)
        end
        pwc[uf].=view(pw,uf)
        rsl = (ppl = pplc, qw = qwc, pw = pwc,  PM = PM[1:nc,:], SW)
        return rsl
    end

    function cIWC(; qw = qw0, pw = pw0, kp = kp0, he = he0, uf = uf0, wc = wc0)
        #Расчёт матрицы связи скважин
        GM.=kp.*he.*10. *8.64*1e-3;
        AG, T = makeAG(kp.*10. *8.64*1e-3,he)
        PM0 .= P0
        wct = view(wc, w2, 1)
        uft = view(uf, :, 1)
        updA!(A,W1,AG,view(rc,:,1),view(rc,:,2),nc,nw,T,λbc,w1,w2,GM,WI,wct, uft, prp.eVp)
        rAdf, rBdf = make_reduce_ma3x_dims(A, w1, w2, nc)
        #AM = transpose(convert(Array{Float32,2}, ACL\tM.M2Mw))
        for t=1:nt
            uft = view(uf, :, t)
            stream_flag .= view(PM0,view(rc,:,1)) .> view(PM0,view(rc,:,2))
            Tp[stream_flag] = Mbt[view(rc,stream_flag,1)]
            stream_flag .= view(PM0,view(rc,:,1)) .< view(PM0,view(rc,:,2))
            Tp[stream_flag] = Mbt[view(rc,stream_flag,2)]
            stream_flag .= view(PM0,view(rc,:,1)) .== view(PM0,view(rc,:,2))
            Tp[stream_flag] = 2*Mbt[view(rc,stream_flag,1)].*Mbt[view(rc,stream_flag,2)]./(Mbt[view(rc,stream_flag,1)].+Mbt[view(rc,stream_flag,2)])

            WTp.=1f0;
            WTp[qw[:,t].>0.] .= Mbt[w1[qw[:,t].>0.]]

            bnd_stream_flag .= view(PM0,bnd_ind) .> Paq
            Tpa[bnd_stream_flag] = view(Mbt[bnd_ind], bnd_stream_flag)
            Tpa[.!bnd_stream_flag] = satc.fkrp.w.(ones(Float32, count(.!bnd_stream_flag)))

            wct = view(wc, w2, t)
            updA!(A,W1,AG.*Tp,view(rc,:,1),view(rc,:,2),nc,nw,T,λbc,w1,w2,GM,WI, WTp.*wct,uft,prp.eVp)
            ACL = cholesky(-A)

            Tpa1[bnd_ind] .= T[bnd_ind].*Tpa;
            AIW[t] = rAdf(-A, Tpa1, λbc)

            PM[:,t], pwc[:,t], pplc[:,t], qwc[:,t] = sim_step!(PM0, qcl, ACL, bb,
                            nc,nw,Paq,T,well,
                            uft,view(qw,:,t), view(pw,:,t),
                            λbc, WI.*WTp, wct, prp.eVp, tM, w1, w2)

            SW[:,t], Mbt[:] = satc(PM0, view(qwc,:,t), GM, AG, false)
        end
        BIW = 0.0#rBdf(bb)
        return AIW, BIW
    end

    return msim, cIWC
end

function make_gdm_prop(;bet0 = 1e-4,
                        Paq0 = 20,
                        λb0 = 2.0)
    bet = bet0;
    Paq = Paq0;
    P0 = Paq0;
    dt = 30.5;
    λb = λb0;
    return (bet = bet, Paq = Paq, P0 = P0, dt = dt, λb = λb)
end

function make_grid(nx,ny,Lx,Ly)
     nc = nx*ny;
     xr = range(0, Lx, length=nx+1)
     x = xr[1:end-1].+step(xr)/2

     yr = range(0, Ly, length=ny+1)
     y = yr[1:end-1].+step(yr)/2

     dx = step(xr)
     ds = step(xr)
     #x = range(dx/2; length=nx, step = dx)
     #y = range(dx/2; length=ny, step = dx)

     XY = collect(Iterators.product(x,y))[:];
     X = getindex.(XY,1)
     Y = getindex.(XY,2)
     rc = make_rc(nx,ny);

     irc = findall(.|(rc[:,1].<=nx,rc[:,1].>nx*(ny-1),mod.(rc[:,1],nx).==1,mod.(rc[:,1],ny).==0))
     λbi = sort(unique(rc[irc,1]))

     return (nc = nc, nx = nx, ny = ny, dx=dx, ds=ds, X = X, Y=Y, rc = rc, λbi=λbi)
end

function make_rc(nx,ny)
    ic = collect(1:nx*ny)
    rc = vcat(hcat(ic,ic.+1),hcat(ic,ic.-1),hcat(ic,ic.+nx),hcat(ic,ic.-nx))
    rc = rc[.&(rc[:,2].>0,rc[:,2].<=nx*ny),:];
    rc = rc[.!.&(mod.(rc[:,1],nx).==0,mod.(rc[:,2],ny).==1),:];
    rc = rc[.!.&(mod.(rc[:,1],nx).==1,mod.(rc[:,2],ny).==0),:];
    return rc
end

function make_gdm(;he_init = 1.,
                   kp_init = 0.2,
                   mp_init = 0.2,
                   nt_init = 360,
                   nx_init = 21,
                   ny_init = 21,
                   Lx_init = 1000,
                   Ly_init = 1000,
                   bet = 1e-4,
                   Paq = 20,
                   λb = 2.0)
    #Создаём всё что надо
    nx, ny, nt = nx_init, ny_init, nt_init;
    Lx, Ly = Lx_init, Ly_init;
    grd = make_grid(nx,ny,Lx,Ly); #кол-во ячеек x, кол-во ячеек y, размер X, размер Y
    gdm_p = make_gdm_prop(bet0 = bet, Paq0 = Paq, λb0 = λb)

    kp = kp_init*ones(grd.nx,grd.ny);   kp = kp[:];
    he = he_init*ones(grd.nx,grd.ny);    he = he[:];
    mp = mp_init*ones(grd.nx,grd.ny);    mp = mp[:];

    #Эффективный поровый объём ячеек (упругоёмкость)
    eVp = gdm_p.bet.*he.*grd.ds.*grd.ds/gdm_p.dt;
    #Поровый объём ячеек
    Vp = he.*grd.ds.*grd.ds.*mp

    prp = (kp = kp, he = he, mp = mp, eVp = eVp, Vp = Vp)

    x = make_well_grid(grd, 0.2, 3)
    return grd, gdm_p, prp, x, nt
end

function make_well_grid(grd, a = 0.25, ncol = 3)
    #a - отступ от края
    #ncol - кол-во рядов
    stp = (1 - a*2)/(ncol-1)
    x = range(grd.dx*grd.nx*a; length=ncol, step = grd.dx*grd.nx*stp)
    return x
end
function make_well(wxy,grd)
    XY = hcat(grd.X,grd.Y)
    well = make_Won(wxy,XY)
    return [well[i,:] for i in 1:length(wxy)]
end

function make_Won(wxy,XY)
    kdtree = KDTree(XY')
    point = hcat(getindex.(wxy,1), getindex.(wxy,2))'
    idx, d = nn(kdtree, point)
    return hcat(getindex.(idx,1),1:length(wxy))
end

function make_Won(wixy::Vector{Tuple{Int64, Tuple{Float64, Float64}}}, XY)
    kdtree = KDTree(XY')
    wi = getindex.(wixy,1)
    wxy = getindex.(wixy,2)
    point = hcat(getindex.(wxy,1), getindex.(wxy,2))'
    idx, d = nn(kdtree, point)
    return hcat(getindex.(idx,1),wi)
end

function make_M2M(grd, well)
    M2M = zeros(grd.nc,9)
    for (km,v) in enumerate(getindex.(well,1))
        vi = grd.rc[v.==grd.rc[:,2],1]
        vj = zeros(Int64,4,4)
        for (k,v1) in enumerate(vi)
            vj[:,k] = grd.rc[v1.==grd.rc[:,2],1]
        end
        gh = zeros(length(unique(vj)))
        for (k,v2) in enumerate(unique(vj))
            gh[k] = sum(v2.==vj)
        end

        M2M[vcat(unique(vj)[gh.>1],vi),km] .= 1/9
    end
    return M2M
end

function makeA(r,c,nc,nw,w1,w2)
    A2 = zeros(nc)
    A = sparse(r,c,1.0,nc+nw,nc+nw);
    accumarray!(A2,r,1)
    updatesp!(A,1:nc,1:nc,.-A2.-1)
    updatesp!(A,nc+1:nc+nw,nc+1:nc+nw,-1.0)
    updatesp!(A,w1,nc.+w2,1.0)
    updatesp!(A,nc.+w2,w1,1.0)

    W1 = sparse(w1,w2,1.,nc,nw)
    return A, W1
end


function makeB!(b, nx,nw,Pk,T,well,uf,qw,pw,λb,p0,WI, wct, eV=0)
    b .= zero(eltype(b));
    wAi = nx+1:nx+nw
    b[1:nx] = .-T.*λb.*Pk
    # for (k,v) in enumerate(zip(well,qw))
    #     b[v[1][1]] = b[v[1][1]] + v[2]
    # end
    b[1:nx].=view(b,1:nx).-eV.*p0
    w1 = getindex.(well,1)
    w2 = getindex.(well,2)
    uff = uf[w2]
    w1f = w1[uff]
    PRDC = wct[uff].*T[w1f].*WI[uff]
    b[w1f] .= b[w1f] .- PRDC.*pw[w2[uff]]
    b[wAi[.!uf]].=qw[.!uf]
    return nothing
end

function make_fun_AG(nc,rc,dx,ds)
    r = view(rc,:,1);
    c = view(rc,:,2);
    return make_fun_AG(nc,r,c,dx,ds)
end

function make_fun_AG(nc,r,c,dx,ds)
    T = zeros(Float32, nc);
    Tr = view(T,r);
    Tc = view(T,c);
    AG = 2 .* T[r] .* T[c] ./ (T[r] .+ T[c]);
    function makeAG(kp,h)
        T.=kp.*h.*ds/dx;
        AG .= 2 .* Tr .* Tc ./ (Tr .+ Tc);
        return AG, T
    end
    return makeAG
end


function updA!(A,W1,AG,r,c,nx,nw,T,λb,w1,w2,GM, WI, wct, uft, eV=0)
    updatesp!(A,r,c,AG);
    #println("---------")
    #println(issymmetric(A))
    A2 = zeros(nx)
    W3 = zeros(nw)
    accumarray!(A2,r,AG)
    A2 .= A2 .+ T.*λb.+eV;
    WIg = WI.*view(GM,w1).*wct
    A2[w1] .= A2[w1] .+ WIg
    A2[w1[uft[w2]]] .= A2[w1[uft[w2]]] .+ WIg[uft[w2]]
    updatesp!(A,1:nx,1:nx,.-A2)
    accumarray!(W3, w2, WIg)
    updatesp!(A,nx+1:nx+nw,nx+1:nx+nw, -W3)
    updatesp!(A,w1,nx.+w2,WIg)
    updatesp!(A,nx.+w2,w1,WIg)
    updatesp!(W1,w1,w2,WIg)
    return nothing
end

function updAS!(A,AG,r,c,nx,T,λb, eV=0)
    updatesp!(A,r,c,AG);
    A2 = zeros(nx)
    accumarray!(A2,r,AG)
    A2 .= A2 .+ T.*λb.+eV;
    updatesp!(A,1:nx,1:nx,.-A2)
    return nothing
end

function make_reduce_ma3x_dims(ACL, w1, w2, nc)
    #nwc = length(w1)
    nw = length(unique(w2))
    w1f = w1
    w2f = indexin(unique(w2),w2)
    aw1 = setdiff(1:nc+nw,w1)
    aw2 = setdiff(1:nc,w1)
    Aaq = zeros(nc)
    tmp = zeros(Float32, nw+1, nc)
    #AA1 = A11 - (A22\Matrix(A12'))'*A21
    #bb1 = bb[w1] - (A22\Matrix(A12'))'*bb[aw1]

    # ia = setdiff(1:441,w1)
    # v12 = zeros(441)
    # v12[1:432] .= (T.*λbc)[ia]
    #
    # AB12 = vcat(A12, v12')'
    # AB21 = hcat(A21, v12)
    function rAdf(A, T, λbc)
        #AA = sparse(ACL)
        AA = -A;
        A11 = AA[w1,w1];
        A22 = AA[aw1,aw1];
        A12 = AA[w1f,aw1]
        A21 = AA[aw1,w1]
        Aaq[1:(nc-nw)] .= -T[aw2].*λbc[aw2]
        A22 = AA[aw1,aw1];

        A11 = vcat(hcat(A11,zeros(nw)),hcat(zeros(nw)',-sum(Aaq)))
        A12 = vcat(AA[w1f,aw1], Aaq')
        A21 = hcat(AA[aw1,w1f], Aaq)
        #A22

        #APA = vcat(hcat(A11, (T.*λbc)[w1]), hcat((T.*λbc)[w1]', sum(T.*λbc))) - (A22\Matrix(AB12))'*AB21
        tmp .= (A22\Matrix(A12'))'
        AA1 = A11 - tmp*A21
        return AA1
    end
    function rBdf(bb)
        #APA = vcat(hcat(A11, (T.*λbc)[w1]), hcat((T.*λbc)[w1]', sum(T.*λbc))) - (A22\Matrix(AB12))'*AB21
        bb1 = bb[w1f] - tmp*bb[aw1]
        return bb1
    end
    return rAdf, rBdf
end


function make_M2W(grd, well)
    #Маска для осреднения давления по площади
    nc = grd.nc
    nw = length(unique(getindex.(well,2)))

    w1 = getindex.(well,1)
    w2 = getindex.(well,2)

    XY = hcat(grd.X,grd.Y)
    kdtree = KDTree(XY')

    point = vcat(view(grd.X,w1)', view(grd.Y, w1)')
    #point = hcat(getindex.(wxy,1), getindex.(wxy,2))'
    RR = getindex.(getindex(knn(KDTree(point), point, 2),2), 1)./2;

    w1g = zeros(Int64, 0)
    w2g = zeros(Int64, 0)
    vg = zeros(Float64, 0)
    for iw = 1:nw
        RRm = minimum(RR[w2.==iw])
        idxp = zeros(Int64, 0)
        for v in eachcol(point[:, w2.==iw])
            idx = inrange(kdtree, v, RRm)
            append!(idxp, idx)
        end
        idx = unique(idxp)
        #scatter(XY[idx, 1], XY[idx, 2])
        append!(w1g, idx)
        append!(w2g, fill(iw, length(idx)))
        append!(vg, fill(inv(length(idx)), length(idx)))
    end
    M2W = sparse(w1g, w2g, vg, nc, nw)
    return M2W
end
