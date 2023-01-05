#This section aims to calculate the topological invariants.
using LinearAlgebra

#From modular Hamiltonian to correlation matrix.
function correlation(G::Matrix)
    G = Matrix(transpose(G))
    v, U = eigen(Hermitian(G))
    #return inv(exp(G1)+diagm(ones(l)))
    return U*diagm(1 ./(exp.(2*v).+1))*U'
end

#From correlation matrix to modular Hamiltonian.
function modular(C::Matrix)
    C = Hermitian(Matrix(transpose(C)))
    v, U = eigen(C)
    G = U*diagm(0.5*log.(1 ./v.-1))*U'
    return G
end

#Calculate the Ensemble Geometric Phase.
function T(Cor::Matrix)
    l = Int(size(Cor)[1]/2)
    dk = 2*pi/l
    x = zeros(2*l)
    for i in 1:l
        x[2*i-1] = i
        x[2*i] = i
    end
    t = diagm(exp.(1im*dk*x))
    Id = diagm(ones(2*l))
    return det(Id.-Cor.+Cor*t)
end

function EGP(Cor::Matrix)
    return imag(log(T(Cor)))
end

#Calculate the Fu-Kane invariant.
#Fourier transformation
function Gk(G::Matrix, N::Int, kx::Float64, ky::Float64)
    S = N*N
    vs = zeros(ComplexF64, 4*S, 4)
    for i in 1:N
        for j in 1:N
            p = N*(i-1)+j
            z = exp(1im*kx*j+1im*ky*i)
            vs[p,1] = z
            vs[p+S,2] = z
            vs[p+2*S,3] = z
            vs[p+3*S,4] = z
        end
    end
    gk = zeros(ComplexF64, 4, 4)
    for i in 1:4
        for j in 1:4
            gk[i,j] = vs[:,i]'*G*vs[:,j]
        end
    end
    return gk
end
     
#Construct Q matrix, given modular Hamiltonian
function Qk(G::Matrix, N::Int, kx::Float64, ky::Float64)
    #G = modular(C)
    gk = Hermitian(Gk(G,N,kx,ky))
    v, U = eigen(gk)
    Q = U*diagm(sign.(v))*U'
    return Q
end

#Block_offdiagonalize
function trans(Q::Matrix)
    #U = [1 0 1 0;0 -1 0 1;0 1 0 1;1 0 -1 0]
    U = [0 1 1 0;-1 0 0 1;1 0 0 1;0 1 -1 0]
    return U'*Q*U./2
end

#Calculate Fu-Kane invariant, given modular Hamiltonian
function fk(G::Matrix, N::Int)
    c = 1
    for i in 0:1
        for j in 0:1
            qk = trans(Qk(G,N,float(pi*i),float(pi*j)))
            c = c*sign(real(qk[1,4]))#Pfaffian
        end
    end
    return c
end