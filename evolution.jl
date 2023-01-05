#This section aims to calculate the time evolution.
using LinearAlgebra

#forward
f(u, p, t) = p[1]*u+u*(p[1])'+2*p[2]

#Runge-Kutta 4 method
function rk4(C0::Matrix, tspan::Vector, p)
    l = length(tspan)
    Cs = Array{Matrix}(undef, l)
    Cs[1] = C0
    dt = tspan[2] - tspan[1]
    for i in 2:l
        Ci = Cs[i-1]
        k1 = f(Ci, p, 0)
        k2 = f(Ci+dt/2*k1, p, 0)
        k3 = f(Ci+dt/2*k2, p, 0)
        k4 = f(Ci+dt*k3, p, 0)
        Cs[i] = Ci + dt/6*(k1+2*k2+2*k3+k4)
    end
    Cs
end

#Example I
function evolution1(N::Int, t::Float64, tspan::Array, γs...)
    X, Mg = damp1(N, γs...)
    C0 = ComplexF64.(correlation(SSH(t,N)))
    Cs = rk4(C0, tspan, [X, Mg])
    
    ps = zeros(length(tspan))
    for i in 1:length(tspan)
        ps[i] = EGP(Cs[i])
    end
    return Cs, fs
end

#Example II
function evolution2(N::Int, μ::Float64, Δ::Float64, tspan::Array, γs...)
    X, Mp = damp2(N, γs...)
    C0 = ComplexF64.(correlation(SC(N,μ,Δ)))
    Cs = rk4(C0, tspan, [X, Mp])
    
    fs = zeros(length(tspan))
    for i in 1:length(tspan_s)
        G = modular(Cs[T*(i-1)+1])
        fs[i] = fk(G,N)
    end
    return Cs, fs
end