# compare_highs.jl -- JuMP + HiGHS side of the comparison.
# Packages:  import Pkg; Pkg.add(["JuMP", "HiGHS"])
using JuMP, HiGHS, Random, LinearAlgebra
include("simplex.jl")

# ---------------------------------------------------------------- 1. crop-mix LP (worked example)
c = [320.0, 410.0, 280.0]                      # profit per ha: maize, cotton, soya (USD)
A = [1.0 1.0 1.0; 3.0 5.0 2.0; 60.0 110.0 45.0]
b = [100.0, 400.0, 9000.0]                     # land (ha), labour (days), fertiliser (USD)

model = Model(HiGHS.Optimizer); set_silent(model)
@variable(model, x[1:3] >= 0)
@constraint(model, land,   sum(A[1, j] * x[j] for j in 1:3) <= b[1])
@constraint(model, labour, sum(A[2, j] * x[j] for j in 1:3) <= b[2])
@constraint(model, fert,   sum(A[3, j] * x[j] for j in 1:3) <= b[3])
@objective(model, Max, sum(c[j] * x[j] for j in 1:3))
optimize!(model)
println("HiGHS objective = ", objective_value(model), "   x = ", value.(x))
println("shadow prices   = ", shadow_price.([land, labour, fert]))
println("reduced costs   = ", reduced_cost.(x))
report = lp_sensitivity_report(model)          # ranges are changes (delta_low, delta_high)
for (j, name) in enumerate(["maize", "cotton", "soya"])
    println("cost range  $name: ", report[x[j]])
end
for (con, name) in [(land, "land"), (labour, "labour"), (fert, "fertiliser")]
    println("RHS  range  $name: ", report[con])
end

# own simplex on the same model (min form, slack columns appended)
As = hcat(A, Matrix(1.0I, 3, 3)); cs = vcat(-c, zeros(3))
r = simplex_tableau(As, b, cs; rule = :dantzig)
println("own simplex: ", r.status, "  objective = ", -r.obj, "  x = ", r.x[1:3], "  pivots = ", r.iterations)

# ---------------------------------------------------------------- 2. Klee-Minty cube (Dantzig rule)
function klee_minty(n)
    cc = [2.0^(n - j) for j in 1:n]
    AA = zeros(n, n); bb = zeros(n)
    for i in 1:n
        for j in 1:i-1
            AA[i, j] = 2.0^(i - j + 1)
        end
        AA[i, i] = 1.0
        bb[i] = 5.0^i
    end
    return cc, AA, bb
end
println("\nKlee-Minty pivots (own Dantzig / own Bland / expected 2^n-1)")
for n in 2:10
    cc, AA, bb = klee_minty(n)
    Ast = hcat(AA, Matrix(1.0I, n, n)); cst = vcat(-cc, zeros(n))
    rd = simplex_tableau(Ast, bb, cst; rule = :dantzig)
    rb = simplex_tableau(Ast, bb, cst; rule = :bland)
    println(n, "  ", rd.iterations, "  ", rb.iterations, "  ", 2^n - 1)
end

# ---------------------------------------------------------------- 3. scaling on random sparse LPs
println("\nScaling: m, n, own pivots, own s, HiGHS simplex s, HiGHS IPM s, relative gap")
Random.seed!(2026)
function solve_highs(cc, AA, bb, solver)
    md = Model(HiGHS.Optimizer); set_silent(md)
    set_attribute(md, "solver", solver)                # "simplex" or "ipm"
    @variable(md, y[1:length(cc)] >= 0)
    @constraint(md, AA * y .<= bb)
    @objective(md, Min, dot(cc, y))
    t = @elapsed optimize!(md)
    return objective_value(md), t
end
for m in (50, 100, 200, 400)
    n = 2m
    AA = rand(m, n) .* (rand(m, n) .< 0.1) .+ 1e-3 .* rand(m, n)
    bb = rand(m) .* 10 .+ 1
    cc = -(rand(n) .* 5 .+ 1)
    Ast = hcat(AA, Matrix(1.0I, m, m)); cst = vcat(cc, zeros(m))
    t_own = @elapsed ro = simplex_tableau(Ast, bb, cst; rule = :dantzig)
    o_ds, t_ds = solve_highs(cc, AA, bb, "simplex")
    o_ip, t_ip = solve_highs(cc, AA, bb, "ipm")
    println(m, ", ", n, ", ", ro.iterations, ", ", round(t_own; digits = 3), ", ",
            round(t_ds; digits = 3), ", ", round(t_ip; digits = 3), ", ",
            abs(ro.obj - o_ds) / max(1, abs(o_ds)))
end
