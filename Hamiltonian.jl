#This section aims to calculate the Hamiltonian and dampin matrix in Example I and II
using LinearAlgebra

#Example I
#Construct SSH model with hopping between unit cell t, N is number of unit cells.
function SSH(t::Float64, N::Int)
    H = zeros(2*N,2*N)
    for i in 1:N-1
        H[2*i-1,2*i] = 1
        H[2*i,2*i-1] = 1
        H[2*i,2*i+1] = t
        H[2*i+1,2*i] = t
    end
    H[2*N-1,2*N] = 1
    H[2*N,2*N-1] = 1
    H[1,2*N] = t
    H[2*N,1] = t
    H = H .+ μ*diagm(ones(2*N))
    return H
end

#Damping matrix.
function damp1(N::Int, μ::Float64, γ1::Float64, γ2::Float64)
    H = diagm(ones(2*N))
    Ml = zeros(2*N,2*N)
    for i in 1:N
        H[2*i-1,2*i-1] = -μ
        H[2*i,2*i] = μ
    end
    
    for i in 1:N-1
        Ml[2*i,2*i] = γ1+γ2
        Ml[2*i,2*i+1] = γ2
        Ml[2*i+1,2*i] = γ2
        Ml[2*i+1,2*i+1] = γ1+γ2
        Ml[2*i,2*i-1] = γ1
        Ml[2*i-1,2*i] = γ1
    end
    US = zeros(2*N,2*N)
    for i in 1:N
        US[2*i-1,2*i-1] = 1.0+0im
        US[2*i,2*i] = -1
    end
    Mg = US*Ml*US
    X = 1im*Matrix(transpose(H)) - Matrix(transpose(Ml))-Mg
    return X, Mg
end

#Example II
#Construct p+ip superconductor Hamiltonian. Hopping term set to be unit. Chemical potential is μ and pairing strength is Δ. System size is NxN.
function SC(N::Int, μ::Float64, Δ::Float64)
    S = N*N
    G = zeros(ComplexF64,4*S,4*S)
    H0 = zeros(S,S)
    Du = zeros(ComplexF64,S,S)
    Dd = zeros(ComplexF64,S,S)
    
    nearx = isnearx(N)
    neary = isneary(N)
    for i in 1:S
        for j in 1:S
            if haskey(nearx,(i,j))
                H0[i,j] = 1.0
                Du[i,j] = nearx[(i,j)]*Δ
                Dd[i,j] = -nearx[(i,j)]*Δ
            elseif haskey(neary,(i,j))
                H0[i,j] = 1.0
                Du[i,j] = neary[(i,j)]*1im*Δ
                Dd[i,j] = neary[(i,j)]*1im*Δ
            end
        end
        H0[i,i] = -μ
    end
    #H0 = H0.-μ*diagm(ones(S))
    
    G[1:S,1:S] = H0/2
    G[S+1:2*S,S+1:2*S] = H0/2
    G[2*S+1:3*S,2*S+1:3*S] = -H0/2
    G[3*S+1:4*S,3*S+1:4*S] = -H0/2
    
    G[1:S,2*S+1:3*S] = Du'
    G[S+1:2*S,3*S+1:4*S] = Dd'
    G[2*S+1:3*S,1:S] = Du
    G[3*S+1:4*S,S+1:2*S] = Dd
    return G
end

function isnearx(N::Int)
    nearx = Dict{Tuple, Int}()
    for i in 1:N
        k = N*i-N
        for j in 1:N
            nearx[(k+j,k+mod(j,N)+1)] = 1#direction, from left bottom to right top.
            nearx[(k+j,k+mod(j-2,N)+1)] = -1
        end
    end
    return nearx
end

function isneary(N::Int)
    neary = Dict{Tuple, Int}()
    for i in 1:N
        k = N*i-N
        k1 = mod(i,N)*N
        k2 = mod(i-2,N)*N
        for j in 1:N
            neary[(j+k,j+k1)] = 1
            neary[(j+k,j+k2)] = -1
        end
    end
    return neary
end

#Damping matrix
function damp6(N::Int, γ1::Float64, γ2::Float64, δ::Float64)
    S = N*N
    M = zeros(ComplexF64, 4*S, 4*S)
    H = zeros(ComplexF64, 4*S, 4*S)
    H0 = zeros(ComplexF64, 2*S, 2*S)
    Mp = zeros(ComplexF64, 4*S, 4*S)
    Ml = zeros(ComplexF64, 2*S, 2*S)
    Mg = zeros(ComplexF64, 2*S, 2*S)
    Mm = zeros(ComplexF64, 2*S, 2*S)
    nearx = isnearx(N)
    neary = isneary(N)
    for i in 1:S
        for j in 1:S
            if haskey(nearx,(i,j))
                H0[i,j] = δ
                H0[S+i,S+j] = -δ
                
                Mg[i,j] = γ2
                Mg[S+i,S+j] = γ2
        
                Mm[i,j] = sqrt(γ1*γ2)*(nearx[(i,j)]+1)/2
                Mm[S+i,S+j] = -sqrt(γ1*γ2)*(nearx[(i,j)]+1)/2
            elseif haskey(neary,(i,j))
                H0[i,j] = δ
                H0[S+i,S+j] = -δ
                
                Mg[i,j] = γ2*1im*neary[(i,j)]
                Mg[S+i,S+j] = -γ2*1im*neary[(i,j)]
                
                Mm[i,j] = 1im*sqrt(γ1*γ2)*(neary[(i,j)]+1)/2
                Mm[S+i,S+j] = 1im*sqrt(γ1*γ2)*(neary[(i,j)]+1)/2
            end
        end
        Mm[i,i] = 2*sqrt(γ1*γ2)
        Mm[S+i,S+i] = -2*sqrt(γ3*γ4)
    end
    Ml = Ml.+2*γ1*diagm(ones(2*S))
    Mg = Mg.+4*γ2*diagm(ones(2*S))
    
    H[1:2*S,1:2*S] = H0
    H[2*S+1:4*S,2*S+1:4*S] = -Matrix(transpose(H0))
    
    M[1:2*S,1:2*S] = Ml
    M[2*S+1:4*S,2*S+1:4*S] = Mg
    M[1:2*S,2*S+1:4*S] = Mm
    M[2*S+1:4*S,1:2*S] = Mm'
    
    Mp[1:2*S,1:2*S] = Mg
    Mp[2*S+1:4*S,2*S+1:4*S] = Ml
    Mp[1:2*S,2*S+1:4*S] = Mm'
    Mp[2*S+1:4*S,1:2*S] = Mm
    return 1im*Matrix(transpose(H))-Matrix(transpose(M))-Mp, Mp
end



