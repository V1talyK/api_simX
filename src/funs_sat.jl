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

function make_calc_krw(n::Int64, alp::Float32, xo::Float32, xw::Float32, mu::Float32)
    #Генератор функции для расчёта фозовой проницаемости
    #n - степень кривизны
    xm = 1f0-xo-xw; #Доля подвижного объёма
    xo1 = 1f0-xo
    inv_mu = 1f0/mu
    inv_xm = 1f0/xm

    xw1 = xw
    n1 = n
    alp1 = alp
    fun = function calc_krw(x::Float32)
        #x = clamp(x,xw,xo1)
        #@fastmath
        # @timeit to1 "krw2_2"  krw = (x-xw)*inv_xm
        # @timeit to1 "krw2_3"  krw = krw^n
        # @timeit to1 "krw2_4"  krw = clamp(krw, 0f0, 1f0)
        # @timeit to1 "krw2_5"  krw = alp*inv_mu*krw
        #@fastmath
        krw = (x-xw1)*inv_xm
        krw = krw^n1
        krw = clamp(krw, 0f0, 1f0)
        krw::Float32 = alp1*inv_mu*krw
        #krw = 1f0
        return krw
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

    fkrp::FKR = getproperty(gdm_sat,:fkrp);
    fkrw::Function = getproperty(fkrp,:w);
    fkro::Function = getproperty(fkrp,:o);

    r = view(grd.rc,:,1)
    c = view(grd.rc,:,2)
    rcs = grd.rc_set
    rci = grd.rc_ind

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
    AWP = zeros(Float32, nc)
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
        accumarray!(AWP, r, abs.(AGP))
        AWP.=AWP./2
        AWP[w1] *= 2
        alpt, alpi = findmax(AWP.*Δt./prp.Vp)
        alp0 = 1/alpt
        #println(alpt)
        alp0 = clamp(alp0,0,1)
        alp = alp0/4

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
            AWP .= 0.0
            @inbounds for i = 1:nc
                Sw_temp = Sw0i[i]
                krw[i] = fkrw(Sw_temp)
                kro[i] = fkro(1f0 - Sw0i[i])
                bale[i] = 1f0/(krw[i]+kro[i])
                A2[i] = prp.eVp[i]*Sw0i[i]*(Pt[i] - PM0[i])#./alp;

                for (kj,j) in enumerate(rcs[i])
                    dP1 = Pt[i] - Pt[j]
                    sdp = sign(dP1)
                    if sdp==1
                        AWP[i] += dP1*krw[i]*AG[rci[i][kj]]
                    #elseif sdp==-1
                        AWP[j] -= dP1*krw[i]*AG[rci[i][kj]]
                    end

                end
            end

            #Tw[fp] = vkrwr#.*baler
            #Tw[bp] = vkrwc#.*balec

            #TAGP .= Tw.*AGP
            #accumarray!(AWP, r, TAGP)
            #A2[sbi] .= view(A2,sbi) .+ sbv.*view(GM,sbi)

            bw[sbi] .= sbv[sbi].*view(GM,sbi)
            bw_aq[bfl] = view(bw,view(sbi,bfl)).*1 .*view(bale,view(sbi,bfl))
            bw_aq[nbfl] = .*(vbwi,vkrwb)#.*(vbwi,vkrwb,baleb)
            bw[sbi] .= bw_aq.*dPaq

            bbw .= 0.0;

            bbw[w1I] .= view(qt,injI)
            bbw[w1P] .= krw[w1P].*view(qt,prodI).*bale[w1P]
            #PF.updatesp!(Aw,1:nc,1:nc,A2)
            #println("www=",sum(krw))
            # @timeit to1 "s1"  begin
            # AWP .= AWP .+A2
            # #println("www=",sum(kro))
            # AWP .= AWP.+bbw
            # AWP .= AWP.- bw
            # AWP .= AWP./prp.Vp
            # end
            #@timeit to1 "s11"
            for i1 = 1:nc
                AWP[i1] = AWP[i1] + A2[i1]
                AWP[i1] = AWP[i1] + bbw[i1]
                AWP[i1] = AWP[i1] - bw[i1]
                AWP[i1] = AWP[i1]/prp.Vp[i1]
            end
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

function make_gdm_prop_sat(;mu_o = 1f0, mu_w = 1f0)
    bet = 1e-4;
    Swaq = 1f0;
    Sw0 = 0f0;

    mu = (o = Float32(mu_o), w = Float32(mu_w))
    xo = 0f0;
    xw = 0f0;
    n_oil = 1# - степень фазовой нефть
    n_wather = 1# - степень фазовой вода
    fkrw, fdkrw = make_calc_krw(n_wather, 1f0, xo, xw, mu.w)
    fkro, fdkro = make_calc_krw(n_oil, 1f0, xw, xo, mu.o)
    fkrp = FKR(fkrw, fkro)
    return GDMS(bet, Swaq, Sw0, fkrp)
end

struct FKR
    w::Function
    o::Function
end

struct GDMS
    bet::Float32
    Swaq::Float32
    Sw0::Float32
    fkrp::FKR
end
