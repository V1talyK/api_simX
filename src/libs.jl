function empty_fun(args...)
    nothing
end

function updatesp!(sp, I, J, newK)
    #Обновляем элементы матрицы A
    @inbounds for (k,(i,j)) in enumerate(zip(I, J));
        sp[i,j] = newK[k];
    end
end

function updatesp!(sp, I, J, newK::Number)
    #Обновляем элементы матрицы A
    @inbounds for (k,(i,j)) in enumerate(zip(I, J));
        sp[i,j] = newK;
    end
end

function accumarray(subs, val, sz=(maximum(subs),))
    A = zeros(eltype(val), sz...)
    accumarray!(A, subs, val)
    return A
end
function accumarray!(A, subs, val)
    A .= eltype(val)(0)
    for i = 1:length(val)
        @inbounds A[subs[i]] += val[i]
    end
end

function gen_rand_qw(nw, nt)
    qw = zeros(nw, nt)
    for i=1:nw
        il = sort(unique(vcat(1,nt,rand(1:nt,7))))
        for k = 1:length(il)-1
            qw[i,il[k]:il[k+1]] .= rand()
        end
    end
    return qw
end
