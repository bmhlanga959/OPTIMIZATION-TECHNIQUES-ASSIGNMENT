################################# PART B, SIMPLEX METHOD ############################
using JuMP, HiGHS

crops  = ["Flue-cured", "Burley", "Maize", "Soyabean", "Groundnuts"]
margin = [1800, 1500, 700, 950, 1100]   # profit $/ha
labr   = [ 120,  100,  40,  45,   70]   # person-days/ha
cap    = [1200, 1000, 300, 400,  500]   # capital $/ha
n = length(crops)

model = Model(HiGHS.Optimizer)
@variable(model, y[1:n] >= 0)

@objective(model, Max, sum(margin[i]*y[i] for i in 1:n))

@constraint(model, land,    sum(y[i] for i in 1:n) <= 40)
@constraint(model, labour,  sum(labr[i]*y[i] for i in 1:n) <= 3600)
@constraint(model, capital, sum(cap[i]*y[i] for i in 1:n) <= 30000)
@constraint(model, barn,    y[1] + y[2] <= 18)   # only tobacco needs barns
@constraint(model, food,    y[3] >= 6)           # minimum maize required

optimize!(model)

for i in 1:n
    println(rpad(crops[i], 12), " => ", round(value(y[i]), digits=2), " ha")
end
println("Max profit Z: \$", round(objective_value(model), digits=2))

println("Shadow prices:")
println("  Land:    ", round(shadow_price(land), digits=2))
println("  Labour:  ", round(shadow_price(labour), digits=2))
println("  Capital: ", round(shadow_price(capital), digits=2))
println("  Barn:    ", round(shadow_price(barn), digits=2))
println("  Food:    ", round(shadow_price(food), digits=2))