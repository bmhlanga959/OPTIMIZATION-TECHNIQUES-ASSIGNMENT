# special_cases.jl -- the four "awkward case" tests (Table 7 of the paper), in Julia.
# This is a direct translation of python/simplex_ref2.py's Beale/alt-optima/unbounded/infeasible
# section. It has NOT been executed for the paper -- run it yourself and compare against Table 7.
using LinearAlgebra
include("simplex.jl")

println("--- Beale cycling example ---")
Ab = [1.0 0 0 .25 -8 -1 9; 0 1.0 0 .5 -12 -.5 3; 0 0 1.0 0 0 1 0]
bb = [0.0, 0.0, 1.0]
cb = [0.0, 0, 0, -.75, 20, -.5, 6.0]
rd = simplex_tableau(Ab, bb, cb; rule = :dantzig, max_iter = 200)
rl = simplex_tableau(Ab, bb, cb; rule = :bland, max_iter = 200)
println("Dantzig rule: ", rd.status, "  iterations = ", rd.iterations)
println("Bland rule:   ", rl.status, "  objective = ", rl.obj, "  iterations = ", rl.iterations)

println("\n--- Alternative optima: max x1 + x2 s.t. x1+x2<=4, x1<=3 ---")
c2 = [1.0, 1.0]
A2 = [1.0 1.0; 1.0 0.0]
b2 = [4.0, 3.0]
As = hcat(A2, Matrix(1.0I, 2, 2))
cs = vcat(-c2, [0.0, 0.0])
r2 = simplex_tableau(As, b2, cs; rule = :dantzig)
println("objective = ", -r2.obj, "  x = ", r2.x[1:2], "  reduced costs = ", r2.reduced_costs)
println("(a zero reduced cost on a nonbasic column indicates a second optimal vertex)")

println("\n--- Unbounded: max x1 + x2 s.t. x1 - x2 <= 1 ---")
Au = reshape([1.0, -1.0], 1, 2)
bu = [1.0]
cu = [1.0, 1.0]
ru = simplex_tableau(hcat(Au, Matrix(1.0I, 1, 1)), bu, vcat(-cu, [0.0]); rule = :dantzig)
println("status = ", ru.status)

println("\n--- Infeasible: x1+x2<=1, x1+x2>=3 ---")
Ai = [1.0 1.0 1.0 0.0; 1.0 1.0 0.0 -1.0]
bi = [1.0, 3.0]
ci = [1.0, 1.0, 0.0, 0.0]
ri = simplex_tableau(Ai, bi, ci; rule = :dantzig)
println("status = ", ri.status)
