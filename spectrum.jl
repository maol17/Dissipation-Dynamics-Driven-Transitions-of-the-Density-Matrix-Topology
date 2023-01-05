#This section aims to calculate the lowest eigenvalues of density matrix
using LinearAlgebra
using Combinatorics

#Compute the largest eigenvalue of density matrix, as well as the excitations of modular Hamiltonian matrix.
function info(C::Matrix, N::Int)
    C0 = Matrix(transpose(Hermitian(C)))
    vs = eigen(C0).values
    es = sort(log.(1 ./vs .-1))
    G = 1.0
    for i in 1:N
        G = G*exp(-es[i])/(1+exp(-es[i]))/(1+exp(-es[2*N+1-i]))
    end
    return G, es
end
    
#Compute the lowest M eigenvalues of density matrix.
function spec_rho(C::Matrix, N::Int, M::Int)
    G, es = info(C, N)
    #gap = 0 - es[N]
    λs = zeros(M)
    λs[1] = G
    Es = es[N+1:M-1+N]

    for i in 2:M-1
        if sum(es[1:i])>Es[M-1]
            break
        end
        z = append!(es[N+1:M-1+N],es[N+1:M-1+N])
        x = sum.(collect(combinations(z,i)))
        Es = sort(append!(x,Es))[1:M-1]
    end
        
    for i in 2:M
        λs[i] = G*exp(-(Es[i-1]))
    end
    λs
end
