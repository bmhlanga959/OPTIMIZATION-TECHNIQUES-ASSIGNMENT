# simplex.jl -- two-phase tableau simplex for  min c'x  s.t.  A x = b, x >= 0
# Pivot rules: :dantzig (largest coefficient) and :bland (smallest index).
# NOTE: written to mirror the verified Python reference implementation; run it yourself and
# compare with the numbers reported in the paper.
using LinearAlgebra

function simplex_tableau(A0::AbstractMatrix, b0::AbstractVector, c::AbstractVector;
                         rule::Symbol = :dantzig, max_iter::Int = 100_000, eps::Float64 = 1e-9)
    A = Matrix{Float64}(A0); b = Vector{Float64}(b0)
    m, n = size(A)
    for i in 1:m                                   # make b >= 0
        if b[i] < 0
            A[i, :] .*= -1; b[i] *= -1
        end
    end
    # ---- Phase I: artificials, minimise their sum
    T = zeros(m + 1, n + m + 1)
    T[1:m, 1:n] = A
    T[1:m, n+1:n+m] = Matrix(1.0I, m, m)
    T[1:m, end] = b
    basis = collect(n+1:n+m)
    T[m+1, n+1:n+m] .= 1.0
    for i in 1:m
        T[m+1, :] .-= T[i, :]
    end
    iters = Ref(0)

    function run!(T, basis, ncols)
        visited = Set{Vector{Int}}()
        while true
            key = sort(basis)
            key in visited && return :cycling
            push!(visited, key)
            row = T[end, 1:ncols]
            j = 0
            if rule == :bland
                idx = findfirst(r -> r < -eps, row)
                idx === nothing && return :optimal
                j = idx
            else
                j = argmin(row)
                row[j] >= -eps && return :optimal
            end
            col = T[1:end-1, j]
            pos = col .> eps
            any(pos) || return :unbounded
            ratios = fill(Inf, length(col))
            ratios[pos] = T[1:end-1, end][pos] ./ col[pos]
            rmin = minimum(ratios)
            ties = findall(r -> r <= rmin + eps, ratios)
            i = rule == :bland ? ties[argmin(basis[ties])] : ties[1]
            T[i, :] ./= T[i, j]
            for r in 1:size(T, 1)
                if r != i && T[r, j] != 0
                    T[r, :] .-= T[r, j] .* T[i, :]
                end
            end
            basis[i] = j
            iters[] += 1
            iters[] > max_iter && return :iterlimit
        end
    end

    st = run!(T, basis, n + m)
    st in (:cycling, :iterlimit) && return (status = st, iterations = iters[])
    T[m+1, end] < -1e-7 && return (status = :infeasible, iterations = iters[])
    for i in 1:m                                   # drive artificials out of the basis
        if basis[i] > n
            nz = findall(v -> abs(v) > eps, T[i, 1:n])
            if !isempty(nz)
                j = nz[1]
                T[i, :] ./= T[i, j]
                for r in 1:m+1
                    r != i && (T[r, :] .-= T[r, j] .* T[i, :])
                end
                basis[i] = j
            end
        end
    end
    # ---- Phase II
    T2 = zeros(m + 1, n + 1)
    T2[1:m, 1:n] = T[1:m, 1:n]; T2[1:m, end] = T[1:m, end]
    T2[m+1, 1:n] = c
    for i in 1:m
        basis[i] <= n && (T2[m+1, :] .-= c[basis[i]] .* T2[i, :])
    end
    st = run!(T2, basis, n)
    st in (:cycling, :iterlimit) && return (status = st, iterations = iters[])
    st == :unbounded && return (status = :unbounded, iterations = iters[])
    x = zeros(n)
    for i in 1:m
        basis[i] <= n && (x[basis[i]] = T2[i, end])
    end
    duals = all(basis .<= n) ? (A[:, basis]' \ c[basis]) : nothing
    return (status = :optimal, x = x, obj = dot(c, x), basis = copy(basis),
            reduced_costs = T2[m+1, 1:n], duals = duals, iterations = iters[])
end
