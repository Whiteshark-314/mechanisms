using LinearAlgebra
using PolynomialRoots
using Statistics

include("utils.jl")

struct Planar3RRR
    a::Vector{Real}
    p::Vector{Real}
    e::Vector{Real}
    f::Vector{Real}
    T::Matrix{Real}
    Θₑ::Vector{Real}
    Θ::Vector{Real}
    Fx
    Fy
    function Planar3RRR(a, p, e, f, rot=0, P=[0, 0])
        @assert length(a)==3 "Provide length of 3 Active Links"
        @assert length(p)==3 "Provide length of 3 Passive Links"
        @assert length(e)==3 "Provide length of 3 sides of End-effector triangle"
        @assert length(f)==3 "Provide length of 3 sides of Fixed/Base triangle"
        @assert a.+p.+√3e/2 > √3f/2 "A + P + √3 E/2 should be > √3 F/2"
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

function Planar3RRR(a::N, p::N, e::N, f::N, rot=0, P=[0, 0]) where N<:Real
    return Planar3RRR([a, a, a], [p, p, p], [e, e, e], [f, f, f], rot, P)
end

function fk(m::Planar3RRR, θ1, θ2, θ3,
    δ::Real=0, λ::NTuple{9, Real}=ntuple(i -> 0, 9))
    ζ = π - m.Θₑ[2]

    ΔC₁ = δ * (cosd(λ[2]) + cosd(λ[5]) + cosd(λ[8]) -
    cosd(λ[1]) - cosd(λ[4]) - cosd(λ[7]))
    ΔS₁ = δ * (sind(λ[2]) + sind(λ[5]) + sind(λ[8]) -
    sind(λ[1]) - sind(λ[4]) - sind(λ[7]))
    X₁ = (m.Fx[2] - m.Fx[1] + m.a[2] * cosd(θ2) - m.a[1] * cosd(θ1) + ΔC₁) ./ m.p[1]
    Y₁ = (m.Fy[2] - m.Fy[1] + m.a[2] * sind(θ2) - m.a[1] * sind(θ1) + ΔS₁) ./ m.p[1]
    pp₁ = m.p[2] ./ m.p[1]
    ep₁ = -m.e[1] ./ m.p[1]

    ΔC₃ = δ * (cosd(λ[2]) + cosd(λ[5]) + cosd(λ[8]) -
    cosd(λ[3]) - cosd(λ[6]) - cosd(λ[9]))
    ΔS₃ = δ * (sind(λ[2]) + sind(λ[5]) + sind(λ[8]) -
    sind(λ[3]) - sind(λ[6]) - sind(λ[9]))
    X₃ = (m.Fx[2] - m.Fx[3] + m.a[2] * cosd(θ2) - m.a[3] * cosd(θ3) + ΔC₃) ./ m.p[3]
    Y₃ = (m.Fy[2] - m.Fy[1] + m.a[2] * sind(θ2) - m.a[3] * sind(θ3) + ΔS₃) ./ m.p[3]
    pp₃ = m.p[2] ./ m.p[3]
    ep₃ = m.e[2] ./ m.p[3]

    a11 = 2 * X₁ * pp₁ + 2 * pp₁ * ep₁
    a12 = 0
    a13 = 2 * X₁ * pp₁ - 2 * pp₁ * ep₁
    b11 = 2 * Y₁ * pp₁
    b12 = 4 * pp₁ * ep₁
    b13 = 2 * Y₁ * pp₁
    c11 = X₁^2 + Y₁^2 + pp₁^2 + ep₁^2 + 2 * X₁ * ep₁ - 1
    c12 = 4 * Y₁ * ep₁
    c13 = X₁^2 + Y₁^2 + pp₁^2 + ep₁^2 - 2 * X₁ * ep₁ - 1

    a31 = 2 * X₃ * pp₃ + 2 * pp₃ * ep₃ * cosd(ζ)
    a32 = -4 * pp₃ * ep₃ * sind(ζ)
    a33 = 2 * X₃ * pp₃ - 2 * pp₃ * ep₃ * cosd(ζ)
    b31 = 2 * Y₃ *pp₃ + 2 * pp₃ * ep₃ * sind(ζ)
    b32 = 4 * pp₃ * ep₃ * cosd(ζ)
    b33 = 2 * Y₃ * pp₃ - 2 * pp₃ * ep₃ * sind(ζ)
    c31 = X₃^2 + Y₃^2 + pp₃^2 + ep₃^2 + 2 * X₃ * ep₃ * cosd(ζ) + 2 * Y₃ * ep₃ * sind(ζ) - 1
    c32 = -4 * X₃ * ep₃ * sind(ζ) + 4 * Y₃ * ep₃ * cosd(ζ)
    c33 = X₃^2 + Y₃^2 + pp₃^2 + ep₃^2 - 2 * X₃ * ep₃ * cosd(ζ) - 2 * Y₃ * ep₃ * sind(ζ) - 1

    p0 = ((a13 * c33 - a33 * c13)^2 - (a13 * b33 - a33 * b13)^2 +
    (b13 * c33 - b33 * c13)^2)
    p1 = (2*(a13*c33 - a33*c13)*(a12*c33 + a13*c32 - a32*c13 - a33*c12)
    - 2*(a13*b33 - a33*b13)*(a12*b33 + a13*b32 - a32*b13 - a33*b12)
    + 2*(b13*c33 - b33*c13)*(b12*c33 + b13*c32 - b32*c13 - b33*c12))
    p2 = (2*(a13*c33 - a33*c13)*(a11*c33 + a12*c32 + 
    a13*c31 - a31*c13 - a32*c12 - a33*c11)
    - 2*(a13*b33 - a33*b13)*(a11*b33 + a12*b32 + 
    a13*b31 - a31*b13 - a32*b12 - a33*b11)
    + 2*(b13*c33 - b33*c13)*(b11*c33 + b12*c32 + 
    b13*c31 - b31*c13 - b32*c12 - b33*c11)
    - (a12 * b33 + a13 * b32 - a32 * b13 - a33 * b12)^2
    + (a12 * c33 + a13 * c32 - a32 * c13 - a33 * c12)^2
    + (b12 * c33 + b13 * c32 - b32 * c13 - b33 * c12)^2)
    p3 = (2*(a12 * c33 + a13 * c32 - a32 * c13 - a33 * c12) * (a11 * c33 + a12 * c32+ a13 * c31 - a31 * c13 - a32 * c12 - a33 * c11)
    - 2*(a12 * b33 + a13 * b32 - a32 * b13 - a33 * b12) * (a11 * b33 + a12 * b32 + a13 * b31 - a31 * b13 - a32 * b12 - a33 * b11)
    + 2*(b12 * c33 + b13 * c32 - b32 * c13 - b33 * c12) * (b11 * c33 + b12 * c32 + b13 * c31 - b31 * c13 - b32 * c12 - b33 * c11)
    - 2*(a13 * b33 - a33 * b13) * (a11 * b32 + a12 * b31 - a31 * b12 - a32 * b11)
    + 2*(a13 * c33 - a33 * c13) * (a11 * c32 + a12 * c31 - a31 * c12 - a32 * c11)
    + 2*(b13 * c33 - b33 * c13) * (b11 * c32 + b12 * c31 - b31 * c12 - b32 * c11))
    p4 = (2*(a11 * c32 + a12 * c31 - a31 * c12 - a32 * c11) * (a12 * c33 + a13 * c32 - a32 * c13 - a33 * c12)
    - 2*(a11 * b32 + a12 * b31 - a31 * b12 - a32 * b11) * (a12 * b33 + a13 * b32 - a32 * b13 - a33 * b12)
    + 2*(b11 * c32 + b12 * c31 - b31 * c12 - b32 * c11) * (b12 * c33 + b13 * c32 - b32 * c13 - b33 * c12)
    - (a11 * b33 + a12 * b32 + a13 * b31 - a31 * b13 - a32 * b12 - a33 * b11)^2
    + (a11 * c33 + a12 * c32 + a13 * c31 - a31 * c13 - a32 * c12 - a33 * c11)^2
    + (b11 * c33 + b12 * c32 + b13 * c31 - b31 * c13 - b32 * c12 - b33 * c11)^2
    - 2*(a11 * b31 - a31 * b11) * (a13 * b33 - a33 * b13) + 2*(a11 * c31 - a31 * c11) * (a13 * c33 - a33 * c13)
    + 2*(b11 * c31 - b31 * c11) * (b13 * c33 - b33 * c13))
    p5 = (2*(a11 * c32 + a12 * c31 - a31 * c12 - a32 * c11) * (a11 * c33 + a12 * c32 + a13 * c31 - a31 * c13 - a32 * c12 - a33 * c11)
    - 2*(a11 * b32 + a12 * b31 - a31 * b12 - a32 * b11) * (a11 * b33 + a12 * b32 + a13 * b31 - a31 * b13 - a32 * b12 - a33 * b11)
    + 2*(b11 * c32 + b12 * c31 - b31 * c12 - b32 * c11) * (b11 * c33 + b12 * c32 + b13 * c31 - b31 * c13 - b32 * c12 - b33 * c11)
    - 2*(a11 * b31 - a31 * b11) * (a12 * b33 + a13 * b32 - a32 * b13 - a33 * b12)
    + 2*(a11 * c31 - a31 * c11) * (a12 * c33 + a13 * c32 - a32 * c13 - a33 * c12)
    + 2*(b11 * c31 - b31 * c11) * (b12 * c33 + b13 * c32 - b32 * c13 - b33 * c12))
    p6 = (2*(a11 * c31 - a31 * c11) * (a11 * c33 + a12 * c32 + a13 * c31 - a31 * c13 - a32 * c12 - a33 * c11)
    - 2*(a11 * b31 - a31 * b11) * (a11 * b33 + a12 * b32 + a13 * b31 - a31 * b13 - a32 * b12 - a33 * b11)
    + 2*(b11 * c31 - b31 * c11) * (b11 * c33 + b12 * c32 + b13 * c31 - b31 * c13 - b32 * c12 - b33 * c11)
    - (a11 * b32 + a12 * b31 - a31 * b12 - a32 * b11)^2 + (a11 * c32 + a12 * c31 - a31 * c12 - a32 * c11)^2
    + (b11 * c32 + b12 * c31 - b31 * c12 - b32 * c11)^2)
    p7 = (2*(a11 * c31 - a31 * c11) * (a11 * c32 + a12 * c31 - a31 * c12 - a32 * c11)
    - 2*(a11 * b31 - a31 * b11) * (a11 * b32 + a12 * b31 - a31 * b12 - a32 * b11)
    + 2*(b11 * c31 - b31 * c11) * (b11 * c32 + b12 * c31 - b31 * c12 - b32 * c11))
    p8 = (a11 * c31 - a31 * c11)^2 - (a11 * b31 - a31 * b11)^2 + (b11 * c31 - b31 * c11)^2

    T = roots([p0, p1, p2, p3, p4, p5, p6, p7, p8])
    T = real(T[imag(T) .== 0.0])

    α = sort(2*atand.(T))
    ϕ = repeat(zero(α), 1, 3)

    A1 =@. 2 * X₁ * pp₁ + 2 * pp₁ * ep₁ * cosd(α);
    B1 =@. 2 * Y₁ * pp₁ + 2 * pp₁ * ep₁ * sind(α);
    C1 =@. X₁^2 + Y₁^2 + pp₁^2 + ep₁^2 + 2 * ep₁ * (X₁ * cosd(α) + Y₁ * sind(α)) - 1;
    A3 =@. 2 * X₃ * pp₃ + 2 * pp₃ * ep₃ * cosd(α + ζ);
    B3 =@. 2 * Y₃ * pp₃ + 2 * pp₃ * ep₃ * sind(α + ζ);
    C3 =@. X₃^2 + Y₃^2 + pp₃^2 + ep₃^2 + 2 * ep₃ * (X₃ * cosd(α + ζ) + Y₃ * sind(α + ζ)) - 1;

    CP_2 =@.  (B3 * C1 - B1 * C3) / (A3 * B1 - A1 * B3);
    SP_2 =@. -(A3 * C1 - A1 * C3) / (A3 * B1 - A1 * B3);
    ϕ[:, 2] = atand.(SP_2, CP_2);

    CP_1 =@. X₁ + pp₁ * cosd(ϕ[:, 2]) + ep₁ * cosd(α);
    SP_1 =@. Y₁ + pp₁ * sind(ϕ[:, 2]) + ep₁ * sind(α);
    ϕ[:, 1] = atand.(SP_1, CP_1);

    CP_3 =@. X₃ + pp₃ * cosd(ϕ[:, 2]) + ep₃ * cosd(α + ζ);
    SP_3 =@. Y₃ + pp₃ * sind(ϕ[:, 2]) + ep₃ * sind(α + ζ);
    ϕ[:, 3] = atand.(SP_3, CP_3);

    E = repeat(zero(α), 1, 3, 2)

    E[:, 1, 1] =@. m.Fx[1] + m.a[1] * cosd(θ1) + m.p[1] * cosd(ϕ[:, 1]) 
    + δ * (cosd(λ[1]) + cosd(λ[4]) + cosd(λ[7]));
    E[:, 2, 1] =@. m.Fx[2] + m.a[2] * cosd(θ2) + m.p[2] * cosd(ϕ[:, 2]) 
    + δ * (cosd(λ[2]) + cosd(λ[5]) + cosd(λ[8]));
    E[:, 3, 1] =@. m.Fx[3] + m.a[3] * cosd(θ3) + m.p[3] * cosd(ϕ[:, 3]) 
    + δ * (cosd(λ[3]) + cosd(λ[6]) + cosd(λ[9]));
    E[:, 1, 2] =@. m.Fx[1] + m.a[1] * sind(θ1) + m.p[1] * sind(ϕ[:, 1]) 
    + δ * (sind(λ[1]) + sind(λ[4]) + sind(λ[7]));
    E[:, 2, 2] =@. m.Fx[2] + m.a[2] * sind(θ2) + m.p[2] * sind(ϕ[:, 2]) 
    + δ * (sind(λ[2]) + sind(λ[5]) + sind(λ[8]));
    E[:, 3, 2] =@. m.Fx[3] + m.a[3] * sind(θ3) + m.p[3] * sind(ϕ[:, 3]) 
    + δ * (sind(λ[3]) + sind(λ[6]) + sind(λ[9]));

    O = mean(E, dims=2);

    return (O, α, E)
end

function ik(m::Planar3RRR, O, α, 
    δ::Real=0, λ::NTuple{9, Real}=ntuple(i -> 0, 9))
    ι = m.Θₑ[1]
    A₁ =[1 1 1; -1 1 0; -1 0 1];
    A₂=[1 1 1; -1 1 0; -1 0 1];
    b₁=[3 * O[2]; m.e[1] * cosd(α); m.e[3] * cosd(ι + α)];
    b₂=[3 * O[2]; m.e[1] * sind(α); m.e[3] * sind(ι + α)];
    ex = A₁ \ b₁;
    ey = A₂ \ b₂;
    L=reshape(λ, 3, 3);
    C=cos(L);
    S=sin(L);
    eFx=transpose(ex .- m.Fx) + C[:, 2]-C[:, 1];
    eFy=transpose(ey .- m.Fy) + S[:, 2]-S[:, 1];

    A=eFx.^2+eFy.^2+m.a.^2+2*eFx.*m.a-m.p.^2;
    B=-4*eFy.*m.a;
    C=eFx.^2+eFy.^2+m.a.^2-2*eFx.*m.a-m.p.^2;

    T1=roots([A[1] B[1] C[1]]);
    T2=roots([A[2] B[2] C[2]]);
    T3=roots([A[3] B[3] C[3]]);
    T1=T1(imag(T1)==0);
    T2=T2(imag(T2)==0);
    T3=T3(imag(T3)==0);
    if isempty(T1)||isempty(T2)||isempty(T3)
        θ=[];
    else
        θ=[2*atand(T1)';2*atand(T2)';2*atand(T3)'];
        ϕ[1, :]=atand((m.Ex[1]-m.a[1]*cos(θ[1, :]))./(m.Ey[1]-m.a[1]*cos(θ[1, :])));
        ϕ[2, :]=atand((m.Ex[2]-m.a[2]*cos(θ[2, :]))./(m.Ey[2]-m.a[2]*cos(θ[2, :])));
        ϕ[3, :]=atand((m.Ex[3]-m.a[3]*cos(θ[3, :]))./(m.Ey[3]-m.a[3]*cos(θ[3, :])));
        θ[θ<0] = 2*pi + θ[θ<0];
        sort!(θ);
        unique_list=unique(nchoosek([1,2,1,2,1,2],3));
        θ_combi=[θ(1,unique_list(:,1))', θ(2,unique_list(:,2))', θ(3,unique_list(:,3))'];
    end
    return θ_combi
end