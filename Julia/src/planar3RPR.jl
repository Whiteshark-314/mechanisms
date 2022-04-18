using LinearAlgebra
using PolynomialRoots
using Statistics

include("utils.jl")

struct Planar3RPR
    ρₘᵢₙ::Vector{Real}
    e::Vector{Real}
    f::Vector{Real}
    T::Matrix{Real}
    Θₑ::Vector{Real}
    Θ::Vector{Real}
    Fx
    Fy
    function Planar3RPR(ρₘᵢₙ, e, f, rot=0, P=[0, 0])
        @assert length(ρₘᵢₙ)==3 "Provide length of 3 Prismatic Links"
        @assert length(e)==3 "Provide length of 3 sides of End-effector triangle"
        @assert length(f)==3 "Provide length of 3 sides of Fixed/Base triangle"
        @assert ρₘᵢₙ.+√3e/2 > √3f/2 "ρₘᵢₙ + √3 E/2 should be > √3 F/2"
        T = [cosd(rot) -sind(rot) P[1];
             sind(rot)  cosd(rot) P[2];
                     0          0    1]
        Ae = acosd((e[1]^2 + e[3]^2 - e[2]^2) / (2 * e[1] * e[3]))
        Be = acosd((e[1]^2 + e[2]^2 - e[3]^2) / (2 * e[1] * e[2]))
        Ce = acosd((e[2]^2 + e[3]^2 - e[1]^2) / (2 * e[2] * e[3]))
        Θₑ = [Ae, Be, Ce]
        A = acosd((f[1]^2 + f[3]^2 - f[2]^2) / (2 * f[1] * f[3]))
        B = acosd((f[1]^2 + f[2]^2 - f[3]^2) / (2 * f[1] * f[2]))
        C = acosd((f[2]^2 + f[3]^2 - f[1]^2) / (2 * f[2] * f[3]))
        Θ = [A, B, C]
        F = T * [0.0 f[1] f[3] * cosd(Ae);
                 0.0  0.0 f[3] * sind(Ae);
                 1.0  1.0             1.0]
        Fx = F[1, :]
        Fy = F[2, :]
        return new(a, p, e, f, T, Θₑ, Θ, Fx, Fy)
    end
end

function Planar3RPR(ρₘᵢₙ::N, e::N, f::N, rot=0, P=[0, 0]) where N<:Real
    return Planar3RPR([ρₘᵢₙ, ρₘᵢₙ, ρₘᵢₙ], [e, e, e], [f, f, f], rot, P)
end