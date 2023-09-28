function x2prm!(V, x, xind;
               iFun = iFun,
               prp=prp)
    if prp in propertynames(xind)
        V.well .= view(x,getproperty(xind,prp))
    end;
    iFun.makeGW!(V.cell,V.well)
end

function calc_wtc(Sw, qw, well,fkr)
    w1 = getindex.(well,1)
    w2 = getindex.(well,2)

    wtc = fkr.w.inj(view(Sw,w1,:))./(fkr.w.inj(view(Sw,w1,:)).+fkr.o.well(view(Sw,w1,:)))
    return wtc
end

function faz(sw)
    a = 1.
    return a*sw
end

function make_calc_krw(n, alp, xo, xw, mu)
    #Генератор функции для расчёта фозовой проницаемости
    #n - степень кривизны
    xm = 1-xo-xw; #Доля подвижного объёма
    xo1 = 1-xo
    fun = function calc_krw(x, pwr = n, coef = alp)
        x = clamp(x,xw,xo1)
        @fastmath coef*clamp(((x-xw)/xm) ^pwr,0,1)/mu
    end

    dfun = function calc_dkrw(x, pwr = n, coef = alp)
        x = clamp(x,xw,xo1)
        @fastmath coef*pwr/(mu*xm)*clamp(((x-xw)/xm) ^(pwr - 1),0,1)
    end
    return fun, dfun
end

function calc_sat_step(prp, grd, gdm_prop, fkr, gdm_s, well, nt, LF;
                       step_t = Inf)
    Pa, Δt = gdm_prop.Paq, gdm_prop.dt;
    nc = grd.nc

    Sw = zeros(Float32, nc, nt);  Sw.=gdm_s.Sw0;
    Sw0 = zeros(Float32, nc);     Sw0 .= gdm_s.Sw0

    r = view(grd.rc,:,1)
    c = view(grd.rc,:,2)

    Tw = zeros(Float32, length(r))
    bw = zeros(Float32, nc)
    sbi = grd.λbi
    sbv = zeros(Float32, nc);
    sbv[sbi] .= gdm_p.λb;

    bw_aq = zeros(Float32, length(sbi))
    bw_aq .= fkr.w.aq(gdm_s.Swaq)

    makeAG = make_fun_AG(grd.nc,grd.rc,grd.dx,grd.ds);
    GG = zeros(Float64,length(r));

    A2=zeros(Float64,nc);
    AW=zeros(Float64,nc);

    w1 = getindex.(well,1)
    w2 = getindex.(well,2)

    bbw = zeros(nc)
    krw = zeros(Float32, nc)
    kro = zeros(Float32, nc)
    AWP = zeros(nc)
    dPaq = zeros(Float32, length(sbi))
    PM0 = zeros(nc); PM0 .= gdm_p.P0;

    TAGP = zeros(length(r))
    AGP = zeros(Float32, length(r))
    dP = zeros(Float32, length(r))

    bale = zeros(Float32, nc)
    Sw0i = copy(Sw0)
    qc_oil =  zeros(Float32, length(unique(w2)))

    function calc_tstep(Pt, qt, t, GM, AG, dP_du, dq_du, flags)
        #вверх по потоку ячейки
        dP .= view(Pt,r) .- view(Pt,c)
        sdp = sign.(dP)

        fp = findall(sdp.==1)
        bp = findall(sdp.==-1)

        rfp = r[fp]
        cfp = c[bp]
        #вверх по потоку ГУ
        dPaq .= Pa .- view(Pt, sbi)
        bfl = dPaq .> 0f0
        nbfl = findall(dPaq .< 0f0)

        cum_Δt = 0.;
        k=0
        flag = true

        AGP .= .*(AG,dP)
        AWP = accumarray(r, abs.(AGP))./2
        AWP[w1] *= 2
        alpt, alpi = findmax(AWP.*Δt./prp.Vp)
        alp0 = 1/alpt

        vkrwr = view(krw,rfp)
        vkrwc = view(krw,cfp)
        vkrwb = view(krw,view(sbi,nbfl))
        vbwi = view(bw,view(sbi,nbfl))

        baler = view(bale,rfp)
        balec = view(bale,cfp)
        baleb = view(bale,view(sbi,nbfl))

        injI = findall(qt.<0)
        prodI = findall(qt.>0)

        w1I =  view(w1,injI)
        w1P =  view(w1,prodI)
        Sw00 = copy(Sw0i)
        Sw00i = copy(Sw0i)

        @inbounds while flag
            k+=1
            krw .= fkr.w.cell(Sw0i)
            kro .= fkr.o.cell(Sw0i)
            bale .= 1 ./(krw.+kro)

            Tw[fp] = vkrwr.*baler
            Tw[bp] = vkrwc.*balec
            mn = bale[alpi]*krw[alpi]
            mn = clamp(mn,0.01,1)
            alp = alp0/4
            #println(mn)

            TAGP .= .*(Tw,AGP)
            accumarray!(AWP, r, TAGP)
            #A2[sbi] .= view(A2,sbi) .+ sbv.*view(GM,sbi)
            A2 .= prp.eVp.*Sw0i.*(Pt .- PM0)./alp;

            bw_aq .= sbv[sbi].*view(GM,sbi)
            bw_aq[bfl] = view(bw,bfl).*1 .*view(bale,view(sbi,bfl))
            bw_aq[nbfl] = .*(vbwi,vkrwb,baleb)
            bw[sbi] .= bw_aq.*dPaq

            bbw .= 0.0;

            bbw[w1I] .= view(qt,injI)
            bbw[w1P] .= krw[w1P].*view(qt,prodI).*bale[w1P]
            #PF.updatesp!(Aw,1:nc,1:nc,A2)
            AWP .= AWP .+A2
            #println("www=",sum(AWP))
            AWP .= AWP.+bbw
            #println(extrema(Sw0i))
            AWP .= AWP.- bw
            #println(sum(AWP))
            AWP .= AWP./prp.Vp

            cum_Δt+=alp
            # println("t=$t, k=$k")
            if alp == 0
                wer
            end
            #println("wee=",extrema(alp.*Δt .*AWP))
            Sw00i .= Sw0i
            Sw0i .= Sw0i .- alp.*Δt .*AWP
            Sw0i .= clamp.(Sw0i,0.0,1.0)
            #Sw0i .= Swi[:,i]
            flag = (cum_Δt < 1) & (k<1000)

            kro .= fkr.o.cell(Sw0i)
            qc_oil .= view(kro,w1).*qt
            qc_oil[injI] .= 0.0

        end

        Sw[:,t] .= Sw0i
        PM0 .= Pt
        return Sw0i
    end

    return calc_tstep
end

function make_gdm_prop_sat()
    bet = 1e-4;
    Swaq = 1.0;
    Sw0 = 0.0;
    return (bet = bet, Swaq = Swaq, Sw0 = Sw0)
end
