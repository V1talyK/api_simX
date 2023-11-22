
function make_sim(grd, gdm_prop, well, prp, nt)
    nc = grd.nc
    nw = maximum(getindex.(well,2))

    rc = grd.rc
    λbi = grd.λbi
    he = copy(prp.he)

    Paq = gdm_prop.Paq
    P0 = gdm_prop.P0
    dt = gdm_prop.dt

    λbc = zeros(nc);
    p0 = zeros(Float32, nc);
    PM = zeros(Float32, nc+nw,nt)
    pwc = zeros(nw,nt)
    pplc = zeros(nw,nt)
    qwc = zeros(nw,nt)
    pplcBt = zeros(Float32, nw)
    GM = ones(Float32, nc);
    bb = zeros(Float32,nc+nw)

    λbc[λbi] .= gdm_prop.λb;

    tM = (M2M = sparse(getindex.(well,2),getindex.(well,1),1,nw,nc),
          M2Mw = sparse(getindex.(well,1),getindex.(well,2),1,nc+nw,nw))
    actW = trues(nw,nt)
    ufl = falses(nw,nt)

    WI = 2*pi./fill(log(0.14*sqrt(2)*grd.dx/0.05),nw)
    w1 = getindex.(well,1)
    w2 = getindex.(well,2)

    A, W1 = makeA(view(rc,:,1),view(rc,:,2),nc,nw,w1,w2)
    makeAG = make_fun_AG(grd.nc,grd.rc,grd.dx,grd.ds);

    qw0 = zeros(Float32, nw, nt)
    pw0 = ones(Float32, nw, nt)
    kp0 = prp.kp
    he0 = prp.he

    function msim(; qw = qw0, pw = pw0, kp = kp0, he = he0)
        GM.=kp.*he.*10. *8.64*1e-3;
        AG, T = makeAG(kp.*10. *8.64*1e-3,he)
        updA!(A,W1,AG,view(rc,:,1),view(rc,:,2),nc,nw,T,λbc,w1,w2,GM,WI,prp.eVp)
        ACL = cholesky(-A)
        CL = make_CL_in_julia(ACL, Threads.nthreads())
        updateCL!(CL, ACL)
        p0 .= P0
        AM = transpose(convert(Array{Float32,2}, ACL\tM.M2Mw))
        rAdf, rBdf = make_reduce_ma3x_dims(ACL, w1, nc)
        AA1 = rAdf(ACL, T, λbc)

        for t=1:nt
            bb .= makeB(nc,nw,Paq,T,well,qw[:,t],λbc,p0,prp.eVp);
            #bb1 = rBdf(bb)
            #temp_ppl = -AA1\bb1;

            PM[:,t] = ACL\bb;
            PM[:,t] .= .-PM[:,t]
            p0 .= view(PM,1:nc,t)
            pwc[:,t] = view(PM,nc+1:nc+nw,t)
            pplcBt .= tM.M2M*p0;
            pplc[:,t] .= pplcBt;
            qwc[:,t] .= qw[:,t]
            #println(sum(abs,pplcBt.-temp_ppl))

        end
        rsl = (ppl = pplc, qw = qwc, pw = pwc,  PM = PM[1:nc,:], AAr = AA1)
        return rsl
    end

    return msim
end

function make_sim2f(grd, gdm_prop, well, prp, nt, satc)
    nc = grd.nc
    nw = maximum(getindex.(well,2))

    rc = grd.rc
    bnd_ind = grd.λbi
    he = copy(prp.he)

    Paq = gdm_prop.Paq
    P0 = gdm_prop.P0
    dt = gdm_prop.dt

    λbc = zeros(nc);
    p0 = zeros(Float32, nc);
    PM = zeros(Float32, nc+nw,nt)
    SW = zeros(Float32, nc,nt)
    pwc = zeros(nw,nt)
    pplc = zeros(nw,nt)
    qwc = zeros(nw,nt)
    pplcBt = zeros(Float32, nw)
    GM = ones(Float32, nc);
    Mbt = ones(Float32, nc); Mbt .= satc.fkrw.(satc.Sw0i).+satc.fkro.(1f0 .- satc.Sw0i)
    Tp = zeros(Float32, size(rc,1))
    bb = zeros(Float32,nc+nw)

    λbc[bnd_ind] .= gdm_prop.λb;

    tM = (M2M = sparse(getindex.(well,2),getindex.(well,1),1,nw,nc),
          M2Mw = sparse(getindex.(well,1),getindex.(well,2),1,nc+nw,nw))
    actW = trues(nw,nt)
    ufl = falses(nw,nt)

    WI = 2*pi./fill(log(0.14*sqrt(2)*grd.dx/0.05),nw)
    w1 = getindex.(well,1)
    w2 = getindex.(well,2)

    A, W1 = makeA(view(rc,:,1),view(rc,:,2),nc,nw,w1,w2)
    makeAG = make_fun_AG(grd.nc,grd.rc,grd.dx,grd.ds);

    qw0 = zeros(Float32, nw, nt)
    pw0 = ones(Float32, nw, nt)
    kp0 = prp.kp
    he0 = prp.he
    stream_flag = trues(length(Tp))
    Tpa = zeros(Float32, length(bnd_ind))
    Tpa1 = zeros(Float32, nc)
    bnd_stream_flag = trues(length(Tpa))

    AA1 = Vector(undef, nt)

    function msim(; qw = qw0, pw = pw0, kp = kp0, he = he0)
        GM.=kp.*he.*10. *8.64*1e-3;
        AG, T = makeAG(kp.*10. *8.64*1e-3,he)
        p0 .= P0

        updA!(A,W1,AG,view(rc,:,1),view(rc,:,2),nc,nw,T,λbc,w1,w2,GM,WI,prp.eVp)
        ACL = cholesky(-A)
        rAdf, rBdf = make_reduce_ma3x_dims(ACL, w1, nc)
        CL = make_CL_in_julia(ACL, Threads.nthreads())
        #AM = transpose(convert(Array{Float32,2}, ACL\tM.M2Mw))
        @timeit to1 "trloop" for t=1:25#nt
            stream_flag .= view(p0,view(rc,:,1)) .> view(p0,view(rc,:,2))
            Tp[stream_flag] = Mbt[view(rc,stream_flag,1)]
            stream_flag .= view(p0,view(rc,:,1)) .< view(p0,view(rc,:,2))
            Tp[stream_flag] = Mbt[view(rc,stream_flag,2)]
            stream_flag .= view(p0,view(rc,:,1)) .== view(p0,view(rc,:,2))
            Tp[stream_flag] = 2*Mbt[view(rc,stream_flag,1)].*Mbt[view(rc,stream_flag,2)]./(Mbt[view(rc,stream_flag,1)].+Mbt[view(rc,stream_flag,2)])

            WTp = ones(nw)
            WTp[qw[:,t].>0.] .= Mbt[w1[qw[:,t].>0.]]

            bnd_stream_flag .= view(p0,bnd_ind) .> Paq
            Tpa[bnd_stream_flag] = view(Mbt[bnd_ind], bnd_stream_flag)
            # println(satc.fkrp.w.(ones(length(count(.!bnd_stream_flag)))))
            # println(Tpa[.!bnd_stream_flag])
            Tpa[.!bnd_stream_flag] = satc.fkrw.(ones(Float32, count(.!bnd_stream_flag)))

            @timeit to1 "1" updA!(A,W1,AG.*Tp,view(rc,:,1),view(rc,:,2),nc,nw,T,λbc,w1,w2,GM,WI.*WTp,prp.eVp)
            @timeit to1 "2" cholesky!(ACL,-A)

            @timeit to1 "3" Tpa1[bnd_ind] .= T[bnd_ind].*Tpa;
            @timeit to1 "4" AA1[t] = rAdf(-A, Tpa1, λbc)
            @timeit to1 "6" updateCL!(CL, ACL)

            bb .= makeB(nc,nw,Paq,T,well,qw[:,t],λbc,p0,prp.eVp);
            PM[:,t] = ACL\bb;
            PM[:,t] .= .-PM[:,t]
            p0 .= view(PM,1:nc,t)
            pwc[:,t] = view(PM,nc+1:nc+nw,t)
            pplcBt .= tM.M2M*p0;
            pplc[:,t] .= pplcBt;
            qwc[:,t] .= qw[:,t]
            @timeit to1 "sat" SW[:,t], Mbt[:] = satc(p0, view(qwc,:,t), GM, AG, false)
        end
        rsl = (ppl = pplc, qw = qwc, pw = pwc,  PM = PM[1:nc,:], SW, AAr = AA1)
        return rsl
    end

    return msim
end

function make_gdm_prop()
    bet = 1e-4;
    Paq = 10;
    P0 = 10;
    dt = 30.5;
    λb = 1;
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

     rc_set = Vector{Array{Int64,1}}(undef,nc);
     rc_ind = Vector{Array{Int64,1}}(undef,nc);
     for i =1:nc
         rc_set[i] = zeros(Int64,0)
         rc_ind[i] = zeros(Int64,0)
     end
     for i in 1:length(rc[:,1])
         #if rc[i,1]>rc[i,2]
             push!(rc_set[rc[i,1]], rc[i,2])
             push!(rc_ind[rc[i,1]], i)
         #end
     end

     irc = findall(.|(rc[:,1].<=nx,rc[:,1].>nx*(ny-1),mod.(rc[:,1],nx).==1,mod.(rc[:,1],ny).==0))
     λbi = sort(unique(rc[irc,1]))

     return (nc = nc, nx = nx, ny = ny, dx=dx, ds=ds, X = X, Y=Y, rc = rc, λbi=λbi,
            rc_set = rc_set, rc_ind = rc_ind)
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
                   nt_init = 360,
                   nx_init = 21,
                   ny_init = 21,
                   Lx_init = 1000,
                   Ly_init = 1000)
    #Создаём всё что надо
    nx, ny, nt = nx_init, ny_init, nt_init;
    Lx, Ly = Lx_init, Ly_init;
    grd = make_grid(nx,ny,Lx,Ly); #кол-во ячеек x, кол-во ячеек y, размер X, размер Y
    gdm_p = make_gdm_prop()

    kp = kp_init*ones(grd.nx,grd.ny);   kp = kp[:];
    he = he_init*ones(grd.nx,grd.ny);    he = he[:];
    mp = 0.2*ones(grd.nx,grd.ny);    mp = mp[:];

    #Эффективный поровый объём ячеек (упругоёмкость)
    eVp = gdm_p.bet.*he.*grd.ds.*grd.ds/gdm_p.dt;
    #Поровый объём ячеек
    Vp = he.*grd.ds.*grd.ds.*mp

    prp = (kp = kp, he = he, eVp = eVp, Vp = Vp)

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


function makeB(nx,nw,Pk,T,well,qw,λb,p0,eV=0)
    b = zeros(nx+nw);
    b[1:nx] = .-T.*λb.*Pk
    # for (k,v) in enumerate(zip(well,qw))
    #     b[v[1][1]] = b[v[1][1]] + v[2]
    # end
    b[1:nx].=view(b,1:nx)-eV.*p0
    b[nx+1:nx+nw].=qw
    return b
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


function updA!(A,W1,AG,r,c,nx,nw,T,λb,w1,w2,GM, WI, eV=0)
    updatesp!(A,r,c,AG);
    #println("---------")
    #println(issymmetric(A))
    A2 = zeros(nx)
    accumarray!(A2,r,AG)
    A2 .= A2 .+ T.*λb.+eV;
    WIg = WI.*view(GM,w1)
    A2[w1] .= A2[w1] .+ WIg
    updatesp!(A,1:nx,1:nx,.-A2)

    updatesp!(A,nx+1:nx+nw,nx+1:nx+nw,.-WIg)
    updatesp!(A,w1,nx.+w2,WIg)
    updatesp!(A,nx.+w2,w1,WIg)
    updatesp!(W1,w1,w2,WIg)
    return nothing
end

function make_reduce_ma3x_dims(ACL, w1, nc)
    nw = length(w1)
    aw1 = setdiff(1:nc+nw,w1)
    aw2 = setdiff(1:nc,w1)
    Aaq = zeros(nc)
    #AA1 = A11 - (A22\Matrix(A12'))'*A21
    #bb1 = bb[w1] - (A22\Matrix(A12'))'*bb[aw1]

    # ia = setdiff(1:441,w1)
    # v12 = zeros(441)
    # v12[1:432] .= (T.*λbc)[ia]
    #
    # AB12 = vcat(A12, v12')'
    # AB21 = hcat(A21, v12)
    function rAdf(ACL, T, λbc)
        AA = sparse(ACL)
        A11 = AA[w1,w1];
        A22 = AA[aw1,aw1];
        A12 = AA[w1,aw1]
        A21 = AA[aw1,w1]
        Aaq[1:(nc-nw)] .= -T[aw2].*λbc[aw2]
        A22 = AA[aw1,aw1];

        A11 = vcat(hcat(A11,zeros(nw)),hcat(zeros(nw)',-sum(Aaq)))
        A12 = vcat(AA[w1,aw1], Aaq')
        A21 = hcat(AA[aw1,w1], Aaq)
        #A22

        #APA = vcat(hcat(A11, (T.*λbc)[w1]), hcat((T.*λbc)[w1]', sum(T.*λbc))) - (A22\Matrix(AB12))'*AB21
        AA1 = A11 - (A22\Matrix(A12'))'*A21
        return AA1
    end
    function rBdf(bb)
        #APA = vcat(hcat(A11, (T.*λbc)[w1]), hcat((T.*λbc)[w1]', sum(T.*λbc))) - (A22\Matrix(AB12))'*AB21
        bb1 = bb[w1] - (A22\Matrix(A12'))'*bb[aw1]
        return bb1
    end
    return rAdf, rBdf
end
