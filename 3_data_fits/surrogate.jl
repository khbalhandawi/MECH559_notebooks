module surrogate_models
export Surrogate, train, basis, predict, Guassian

using LinearAlgebra
import LinearAlgebra.norm

abstract type AbstractSurrogate{T} end

mutable struct LS{T} <: AbstractSurrogate{T}
    X::Matrix{T}
    Y::Matrix{T}
    r::T
    d::Int
    n_centers::Int
    dim_i::Int
    dim_o::Int
    isempty::Bool
    B::Matrix{T}
    W::Matrix{T}
    J::Matrix{T}

    function LS(::Type{T}=Float64; X::Matrix{T}, Y::Matrix{T}, r::T, d::Int) where T
        
        ls = new{T}()
        ls.X = X
        ls.Y = Y
        ls.r = r
        ls.d = d
        
        # Dimensionality
        ls.n_centers = size(X,1)
        ls.dim_i = size(X,2)
        ls.dim_o = size(Y,2)

        # Weights
        ls.isempty = true
        ls.B = zeros(T,(size(X,1),size(X,2)*d+1))
        ls.W = zeros(T,(size(X,2)*d+1,size(Y,2)))

        # reguralization
        ls.J = I(size(X,2)*d+1) .* r
        ls.J[1,1] = 0
        return ls
    end
end

function Polynomial(d::Int, ζ::Vector{T}) where T
    n = length(ζ)
    b = zeros(T, d*n+1)
    b[1] = 1
    for i = 1:d
        b[2+(i-1)*n:1+i*n] = ζ.^i
    end
    return b
end

function basis(d::Int, Z::Matrix{T}, f::Function) where T

    B = zeros(T, size(Z, 1), size(Z,2)*d+1)

    for k = 1:size(Z,1)
        B[k,:] = f(d,Z[k,:])
    end
    return B
end

function train(m::LS)
    B = basis(m.d,m.X,Polynomial)
    m.B = B
    m.W = pinv(B' * B + m.J) * B' * m.Y
    m.isempty = false
end

function predict(m::LS{T}, Z::Matrix{T}) where T

    if size(Z,2) != size(m.X,2)
        throw(ArgumentError("dimenions of prediction and training sites are not equal"))
    end

    if m.isempty
        throw(DomainError(m,"model is not trained!"))
    end
    basis(m.d,Z,Polynomial)*m.W
end

mutable struct RBF{T} <: AbstractSurrogate{T}
    X::Matrix{T}
    Y::Matrix{T}
    r::T
    λ::T
    kernel::Function
    n_centers::Int
    dim_i::Int
    dim_o::Int
    isempty::Bool
    B::Matrix{T}
    W::Matrix{T}
    J::Matrix{T}

    function RBF(::Type{T}=Float64; X::Matrix{T}, Y::Matrix{T}, r::T, λ::T, kernel::Function) where T
        
        rbf = new{T}()
        rbf.X = X
        rbf.Y = Y
        rbf.r = r
        rbf.λ = λ
        rbf.kernel = kernel
        
        # Dimensionality
        rbf.n_centers = size(X,1)
        rbf.dim_i = size(X,2)
        rbf.dim_o = size(Y,2)

        # Weights
        rbf.isempty = true
        rbf.B = zeros(T,(size(X,1),size(X,1)))
        rbf.W = zeros(T,(size(X,1),size(Y,2)))

        # reguralization
        rbf.J = I(size(X,1)) .* r
        rbf.J[1,1] = 0
        return rbf
    end
end

# basis functions
function norm(X::Matrix{T}, dims::Integer) where T
    sqrt.(sum(X.^2, dims=dims))
end

function Guassian(λ::T, ζ::Vector{T}, X::Matrix{T}) where T
    exp.(-λ*norm(ζ' .- X, 2).^2)
end

function basis(λ::T, Z::Matrix{T}, X::Matrix{T}, f::Function) where T

    B = zeros(T, size(Z, 1), size(X,1))

    for k = 1:size(Z,1)
        B[k,:] = f(λ,Z[k,:],X)
    end
    return B
end

function train(m::RBF)
    B = basis(m.λ,m.X,m.X,m.kernel)
    m.B = B
    m.W = pinv(B' * B + m.J) * B' * m.Y
    m.isempty = false
end

function predict(m::RBF{T}, Z::Matrix{T}) where T

    if size(Z,2) != size(m.X,2)
        throw(ArgumentError("dimenions of prediction and training sites are not equal"))
    end

    if m.isempty
        throw(DomainError(m,"model is not trained!"))
    end
    basis(m.λ,Z,m.X,m.kernel)*m.W
end
    
mutable struct Surrogate{T}
    X::Matrix{T}
    Y::Matrix{T}
    name::String
    type::String
    model::AbstractSurrogate{T}

    function Surrogate(::Type{T}=Float64; X::Matrix{T}, Y::Matrix{T}, type::String, r::T=1e-6, λ::T=15.0, d::Int=1, kernel::Function=Guassian, name::String = "") where T

        if size(X,1) != size(Y,1)
            throw(ArgumentError("number of training inputs does not match number of training outputs"))
        end

        sur = new{T}()
        sur.type = type
        sur.name = name

        if type == "RBF"
            sur.model = RBF(T, X=X , Y=Y, r=r, λ=λ, kernel=kernel)
        end

        if type == "LS"
            sur.model = LS(T, X=X , Y=Y, r=r, d=d)
        end

        return sur
    end
end

end # module