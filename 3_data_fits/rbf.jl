module surrogate_models
export RBF, train, kernel, predict, Guassian

using LinearAlgebra
import LinearAlgebra.norm

mutable struct RBF{T}
    X::Matrix{T}
    Y::Matrix{T}
    r::T
    λ::T
    name::String
    n_centers::Int
    dim_i::Int
    dim_o::Int
    isempty::Bool
    B::Matrix{T}
    W::Matrix{T}
    J::Matrix{T}

    function RBF(::Type{T} = Float64; X::Matrix{T}, Y::Matrix{T}, r::T = 1e-6, λ::T = 0.1, name::String = "", ) where T
        sur = new{T}()
        sur.X = X
        sur.Y = Y 
        sur.r = r
        sur.λ = λ
        sur.name = name
    
        if size(X,1) != size(Y,1)
            throw(ArgumentError("number of training inputs does not match number of training outputs"))
        end
        
        # Dimensionality
        sur.n_centers = size(X,1)
        sur.dim_i = size(X,2)
        sur.dim_o = size(Y,2)

        # Weights
        sur.isempty = true
        sur.B = zeros(T,(size(X,1),size(X,1)))
        sur.W = zeros(T,(size(X,1),size(Y,2)))

        # reguralization
        sur.J = I(size(X,1)) .* r
        sur.J[1,1] = 0
        return sur
    end
end

function norm(X::Matrix{T}, dims::Integer) where T
    sqrt.(sum(X.^2, dims=dims))
end

function Guassian(λ::T, ζ::Vector{T}, X::Matrix{T})where T
    exp.(-λ*norm(ζ' .- X, 2).^2)
end

function kernel(λ::T, Z::Matrix{T}, X::Matrix{T}, f::Function) where T

    if size(Z,2) != size(X,2)
        throw(ArgumentError("dimenions of prediction and training sites are not equal"))
    end

    B = zeros(T, size(Z, 1), size(X,1))

    for k = 1:size(Z,1)
        B[k,:] = f(λ,Z[k,:],X)
    end
    return B
end

function train(m::RBF, kernel_f::Function)
    B = kernel(m.λ,m.X,m.X,kernel_f)
    m.B = B
    m.W = pinv(B' * B + m.J) * B' * m.Y
    m.isempty = false
end

function predict(m::RBF{T}, Z::Matrix{T}, kernel_f::Function) where T
    if m.isempty
        throw(DomainError(m,"model is not trained!"))
    end
    kernel(m.λ,Z,m.X,kernel_f)*m.W
end
    
end