function back_slash_slvr!(x,
                         A::NamedTuple{(:L, :U, :p, :x_temp), Tuple{SparseMatrixCSC{Float64, Int64},
                                                                SparseMatrixCSC{Float64, Int64},
                                                                Vector{Int64},
                                                                Matrix{Float64}}},
                        b)

    ldiv_cl!(x,A,b)
    return nothing
end

function back_slash_slvr!(x,A::SuiteSparse.CHOLMOD.Factor{Float64},b)
    x .= A\b;
    return nothing
end

function make_CL_in_julia(ACL, nth = 1)
    LL = sparse(ACL.L)
    UU = copy(LL')
    x_temp = zeros(LL.n, nth)
    return (L = LL, U = UU, p = ACL.p, x_temp)
end

function updateCL!(CL, ACL)
    CL.L.=sparse(ACL.L)
    CL.U.=copy(CL.L')
end

function ldiv_cl!(x,
                  CL::NamedTuple{(:L, :U, :p, :x_temp),
                      Tuple{SparseMatrixCSC{Float64, Int64},
                      SparseMatrixCSC{Float64, Int64},
                      Vector{Int64},
                      Matrix{Float64}}},
                  b)
    @inbounds bp = view(b,CL.p)
    @inbounds xp = view(x,CL.p)
    x_temp = view(CL.x_temp,:,Threads.threadid())
    forward_substit!(x_temp,CL.L, bp)
    backward_substit!(xp,CL.U, x_temp)
    #invpermute!(x,CL.p)
end

function forward_substit!(x, S, b)
    x .= 0;
    @fastmath @inbounds for col = 1:S.n
        idx = S.colptr[col]+1 : S.colptr[col + 1] - 1
        #println(idx," ",S.rowval[idx])
        x[col] = (b[col] + x[col])/S.nzval[S.colptr[col]]
        #v1 = view(S.rowval,idx)
        for v in idx
             x[S.rowval[v]] -=  S.nzval[v] * x[col]
        end
    end
end

function backward_substit!(x, UU, b)
    x .= 0;
    @fastmath @inbounds for col = UU.n:-1:1
        idx = UU.colptr[col + 1]-2 :-1 : UU.colptr[col]
        #println(idx," ",S.rowval[idx])
        x[col] = (b[col] + x[col])/UU.nzval[UU.colptr[col+1]-1]

        for v in idx
             x[UU.rowval[v]] -=  UU.nzval[v] * x[col]
        end
    end
end

# function sam_solve!(X::Dense{Tv}, Yref::Ref{Ptr{C_Dense{Tv}}}, Eref::Ref{Ptr{C_Dense{Tv}}}, F::Factor{Tv}, B::Dense{Tv}) where Tv<:VTypes
#   # Pointer to pre-allocated dense matrix
#   Xref = Ptr{C_Dense{Tv}}(pointer(X))
#   sys = CHOLMOD_A # Solve type
#   Bset = C_NULL   # Unused parameter
#   Xset = C_NULL   # Unused parameter
#
#   if size(F,1) != size(B,1)
#     throw(DimensionMismatch("LHS and RHS should have the same number of rows. " *
#     "LHS has $(size(F,1)) rows, but RHS has $(size(B,1)) rows."))
#   end
#
#   if !issuccess(F)
#     s = unsafe_load(pointer(F))
#     if s.is_ll == 1
#       throw(LinearAlgebra.PosDefException(s.minor))
#     else
#       throw(LinearAlgebra.ZeroPivotException(s.minor))
#     end
#   end
#
#   res = ccall((@cholmod_name("solve2"), :libcholmod), Cint,
#   (Cint, Ptr{C_Factor{Tv}}, Ptr{C_Dense{Tv}}, Ptr{C_Sparse{Tv}}, Ref{Ptr{C_Dense{Tv}}},  Ref{Ptr{C_Sparse{Tv}}},  Ref{Ptr{C_Dense{Tv}}},  Ref{Ptr{C_Dense{Tv}}}, Ptr{UInt8}),
#   sys,   F, B,Bset,Xref,Xset,Yref,Eref,SuiteSparse.CHOLMOD.common_struct[Threads.threadid()])
#
#   if (res != 1)
#     throw(ErrorException("CHOLMOD solve failure"))
#   end
#
#   return nothing
# end
#
# function copyto1!(dest::Dense{T},src::Vector{T}) where {T<:AbstractFloat}
#   #T = eltype(dest)
#   GC.@preserve dest begin
#     s = unsafe_load(pointer(dest));
#     pt = Ptr{T}(s.x)
#     @inbounds for (i, c) in enumerate(eachindex(src))
#         #@time x[c];
#         unsafe_store!(pt, src[c], i);
#     end
#   end
# end
