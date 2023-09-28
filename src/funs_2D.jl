
function make_gdm(;he_init = 1.)
    #Создаём всё что надо
    nx, ny, nt = 21, 21, 360;
    Lx, Ly = 1000, 1000;
    grd = make_grid(nx,ny,Lx,Ly); #кол-во ячеек x, кол-во ячеек y, размер X, размер Y
    gdm_p = make_gdm_prop()

    kp = 0.2*ones(grd.nx,grd.ny);   kp = kp[:];
    he = he_init*ones(grd.nx,grd.ny);    he = he[:];
    mp = 0.2*ones(grd.nx,grd.ny);    mp = mp[:];

    #Эффективный поровый объём ячеек (упругоёмкость)
    eVp = gdm_p.bet.*he.*grd.ds.*grd.ds/gdm_p.dt;
    #Поровый объём ячеек
    Vp = he.*grd.ds.*grd.ds*0.2

    prp = (kp = kp, he = he, eVp = eVp, Vp = Vp)

    x = make_well_grid(grd, 0.2, 3)
    return grd, gdm_p, prp, x, nt
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
