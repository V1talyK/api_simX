function x2prm!(V, x, xind;
               iFun = iFun,
               prp=prp)
    if prp in propertynames(xind)
        V.well .= view(x,getproperty(xind,prp))
    end;
    iFun.makeGW!(V.cell,V.well)
end

function calc_wtc(Sw, fkrp, well)
    w1 = getindex.(well,1)
    w2 = getindex.(well,2)

    wtc = fkrp.w.(view(Sw,w1,:))./(fkrp.w.(view(Sw,w1,:)).+fkrp.o.(1 .- view(Sw,w1,:)))
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
        @fastmath coef*clamp(((x-xw)/xm) ^pwr, 0.0, 1.0)/mu
    end

    dfun = function calc_dkrw(x, pwr = n, coef = alp)
        x = clamp(x,xw,xo1)
        @fastmath coef*pwr/(mu*xm)*clamp(((x-xw)/xm) ^(pwr - 1),0,1)
    end
    return fun, dfun
end

function calc_sat_step(prp, grd, gdm_prop, gdm_sat, well, nt;
                       step_t = Inf)
    Pa, Δt = gdm_prop.Paq, gdm_prop.dt;
    nc = grd.nc

    Sw = zeros(Float32, nc, nt);  Sw.=gdm_sat.Sw0;
    Sw0 = zeros(Float32, nc);     Sw0 .= gdm_sat.Sw0
    fkrp = gdm_sat.fkrp;

    r = view(grd.rc,:,1)
    c = view(grd.rc,:,2)

    Tw = zeros(Float32, length(r))
    bw = zeros(Float32, nc)
    sbi = grd.λbi
    sbv = zeros(Float32, nc);
    sbv[sbi] .= gdm_prop.λb;

    bw_aq = zeros(Float32, length(sbi))
    bw_aq .= fkrp.w(gdm_sat.Swaq)

    A2=zeros(Float64,nc);
    AW=zeros(Float64,nc);

    w1 = getindex.(well,1)
    w2 = getindex.(well,2)

    bbw = zeros(nc)
    krw = zeros(Float32, nc)
    kro = zeros(Float32, nc)
    AWP = zeros(nc)
    dPaq = zeros(Float32, length(sbi))
    PM0 = zeros(nc); PM0 .= gdm_prop.P0;

    TAGP = zeros(length(r))
    AGP = zeros(Float32, length(r))
    dP = zeros(Float32, length(r))

    bale = zeros(Float32, nc)
    Sw0i = copy(Sw0)

    function calc_tstep(Pt, qt, GM, AG, flags)
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
        alp0 = clamp(alp0,0,1)

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
            krw .= fkrp.w.(Sw0i)
            kro .= fkrp.o.(1.0 .- Sw0i)
            bale .= 1 #./(krw.+kro)

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

            bw[sbi] .= sbv[sbi].*view(GM,sbi)
            bw_aq[bfl] = view(bw,view(sbi,bfl)).*1 .*view(bale,view(sbi,bfl))
            bw_aq[nbfl] = .*(vbwi,vkrwb,baleb)
            bw[sbi] .= bw_aq.*dPaq

            bbw .= 0.0;

            bbw[w1I] .= view(qt,injI)
            bbw[w1P] .= krw[w1P].*view(qt,prodI).*bale[w1P]
            #PF.updatesp!(Aw,1:nc,1:nc,A2)
            #println("www=",sum(krw))
            AWP .= AWP .+A2
            #println("www=",sum(kro))
            AWP .= AWP.+bbw
            AWP .= AWP.- bw
            AWP .= AWP./prp.Vp
            #wer
            cum_Δt+=alp
            # println("t=$t, k=$k")
            if alp == 0
                wer
            end
            #println("wee=",extrema(alp.*Δt .*AWP))
            Sw00i .= Sw0i
            Sw0i .= Sw0i .- alp.*Δt .*AWP
                    Sw0i .= clamp.(Sw0i,0.0,1.0)
            #println(Sw0i)
            #wer
            flag = (cum_Δt < 1) & (k<1000)

        end
        PM0 .= Pt
        #println(cum_Δt)
        #println(sum(Sw0i.*prp.Vp) - sum(Sw00.*prp.Vp),"  ",krw[w1P].*view(qt,prodI).*bale[w1P]*30.5,
        #                                              "  ",kro[w1P].*view(qt,prodI).*bale[w1P]*30.5,
    #                                                  "  ",sum(bw))
        return Sw0i, krw .+ kro
    end

    return calc_tstep
end

function make_gdm_prop_sat(;mu_o = 1., mu_w = 1.)
    bet = 1e-4;
    Swaq = 1.0;
    Sw0 = 0.0;

    mu = (o = mu_o, w = mu_w)
    xo = 0.;
    xw = 0.;
    n_oil = 1# - степень фазовой нефть
    n_wather = 1# - степень фазовой вода
    fkrw, fdkrw = make_calc_krw(n_wather, 1, xo, xw, mu.w)
    fkro, fdkro = make_calc_krw(n_oil, 1, xw, xo, mu.o)

    return (bet = bet, Swaq = Swaq, Sw0 = Sw0, fkrp = (w = fkrw, o = fkro))
end
