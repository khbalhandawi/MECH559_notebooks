function consecutive(f, A::Matrix{Float64})::Matrix{Float64}
    dims = size(A)
    df = Matrix{Float64}(undef,dims[1]-1,dims[2])
    for j = 1:dims[2]
        df[:,j] = [ f(A[i+1,j], A[i,j]) for i = 1:dims[1]-1 ]
    end
    return df
end

function gridsamp(bounds::Matrix{Float64}, q::Vector{Int64})
    mr, n = size(bounds)
    dr = consecutive(-,bounds)
    
    if mr != 2 || any(x->x<0, dr)
        throw(DomainError(bounds,"bounds must be an array with two rows and bounds(1,:) <= bounds(2,:)"))
    end
    
    if ndims(q) > 1 || any(x->x<=0, q)
        throw(DomainError(q,"q must be a vector with non-negative elements"))
    end
    
    p = length(q)
    if p == 1
        q = repeat(q, outer = (1, n))
    elseif p != n
        throw(DomainError(q,"length of q must be either 1 or $n"))
    end
    
    # Check for degenerate intervals
    i = findall(x->x==0, dr)
    if length(i) > 0
        for index in i
            q[index[2]] = 0 * q[index[2]]
        end
    end

    # Recursive computation
    if n > 1
        a = gridsamp(bounds[:, 2:end], q[2:end])  # Recursive call
        m,_ = size(a)
        q = q[1]

        s = hcat(zeros(m * q, 1), repeat(a, outer = (q, 1)))
        y = LinRange(bounds[1, 1], bounds[2, 1], q)

        k = 1:m
        for i = 1:q
            aug = fill(y[i], (m, 1))
            aug = reshape(aug, size(s[k, 1]))

            s[k, 1] = aug
            k = k .+ m
        end
    else
        s = LinRange(bounds[1, 1], bounds[2, 1], q[end])
        s = reshape(s, (:,1))
    end

    return s
end