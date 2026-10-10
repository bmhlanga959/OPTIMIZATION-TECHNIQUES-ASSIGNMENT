using LinearAlgebra, Plots
println("===== PART A: GRAPHICAL/LINEAR METHOD =====\n")

constraints = [
    (1,   1,   40,    "Land"),
    (120, 40,  3600,  "Labour"),
    (1200,300, 30000, "Capital"),
    (1,   0,   0,     "x1-axis"),
    (0,   1,   0,     "x2-axis")
]

c1, c2 = 1800, 700
corners = Tuple{Float64,Float64}[]

n = length(constraints)
for i in 1:n, j in (i+1):n
    a1,b1,_,_ = constraints[i]
    a2,b2,_,_ = constraints[j]
    A = [a1 b1; a2 b2]
    if abs(det(A)) > 1e-9
        rhs = [constraints[i][3], constraints[j][3]]
        pt = A \ rhs
        push!(corners, (pt[1], pt[2]))
    end
end

function feasible(x1, x2)
    x1 >= -1e-6 && x2 >= -1e-6 &&
    x1 + x2 <= 40 + 1e-6 &&
    120x1 + 40x2 <= 3600 + 1e-6 &&
    1200x1 + 300x2 <= 30000 + 1e-6
end

println("Feasible corners and their profit:")
best_z = -Inf
best_pt = (0.0, 0.0)
seen = Set{Tuple{Float64,Float64}}()
feasible_pts = Tuple{Float64,Float64}[]   # collect these for the plot

for (x1, x2) in corners
    x1r, x2r = round(x1, digits=4), round(x2, digits=4)
    if feasible(x1r, x2r) && !((x1r,x2r) in seen)
        push!(seen, (x1r,x2r))
        push!(feasible_pts, (x1r, x2r))
        z = c1*x1r + c2*x2r
        println("  (x1=$x1r, x2=$x2r)  ->  Z = \$$(round(z,digits=2))")
        if z > best_z
            best_z = z
            best_pt = (x1r, x2r)
        end
    end
end

println("\nOptimal corner: x1 = $(best_pt[1]) ha, x2 = $(best_pt[2]) ha")
println("Maximum profit: \$$(round(best_z, digits=2))")

# ---------------- PLOTTING SECTION ----------------

# x-axis range for drawing lines
xs = 0:0.5:40

# Rearrange each constraint into "x2 = ..." form so we can plot it as
# a line: a1*x1 + b1*x2 = c  ->  x2 = (c - a1*x1) / b1
land_line    = [(40 - 1*x)      for x in xs]
labour_line  = [(3600 - 120x)/40   for x in xs]
capital_line = [(30000 - 1200x)/300 for x in xs]

plt = plot(xs, land_line,    label="Land: x1+x2=40",         linewidth=2)
plot!(xs, labour_line,  label="Labour: 120x1+40x2=3600",  linewidth=2)
plot!(xs, capital_line, label="Capital: 1200x1+300x2=30000", linewidth=2)

# Shade the feasible region by sorting corner points around their
# centre (so the polygon draws correctly instead of a crossed shape)
cx = sum(p[1] for p in feasible_pts) / length(feasible_pts)
cy = sum(p[2] for p in feasible_pts) / length(feasible_pts)
sorted_pts = sort(feasible_pts, by = p -> atan(p[2]-cy, p[1]-cx))
region_x = [p[1] for p in sorted_pts]
region_y = [p[2] for p in sorted_pts]
plot!(Shape(region_x, region_y), fillalpha=0.25, label="Feasible region", color=:green)

# Mark every feasible corner with a dot
scatter!([p[1] for p in feasible_pts], [p[2] for p in feasible_pts],
          label="Corners", color=:black, markersize=5)

# Highlight the optimal point in red
scatter!([best_pt[1]], [best_pt[2]], label="Optimal (Z=\$$(round(best_z,digits=0)))",
          color=:red, markersize=8, markershape=:star5)

xlabel!("x1 (Tobacco, ha)")
ylabel!("x2 (Maize, ha)")
title!("Part A: Feasible Region and Optimal Solution")
xlims!(0, 45)
ylims!(0, 45)

display(plt)
savefig(plt, "part_a_graph.png")   # also saves a copy you can insert into your report

