include("gridsamp.jl")
include("surrogate.jl")
using Random, DelimitedFiles, LaTeXStrings, Plots, LinearAlgebra
using .sampling
using .surrogate_models

data = readdlm("./MECH559_notebooks/3_data_fits/sample_data_1.csv", ',', Float64, skipstart=1)
X_data = data[:,3:4]
X_data[:,2] = X_data[:,2]
Y_data = data[:,2:2]

# plot bounds
lb = [1.0 20.0]
ub = [12.0 30.0]

bounds_reg = Matrix{Float64}(undef,2,2)
bounds_reg[1,:] = lb
bounds_reg[2,:] = ub

n_grids = 200
n_reg = [n_grids,n_grids]
X_reg = sampling.gridsamp(bounds_reg,n_reg)
x_reg = LinRange(bounds_reg[1,1],bounds_reg[2,1],n_grids)
y_reg = LinRange(bounds_reg[1,2],bounds_reg[2,2],n_grids)

# train models
ridge = 0.0
degree = 1
ls = surrogate_models.Surrogate(X=X_data, Y=Y_data, type="LS", lb=vec(lb), ub=vec(ub), r=ridge, d=degree, scale=false)
surrogate_models.train(ls.model);

# Predictions
y_hat = surrogate_models.predict(ls.model,X_reg)

println("unscaled: ", cond(ls.model.B'*ls.model.B + ls.model.J))
println("unscaled design matrix: ")
println(ls.model.B'*ls.model.B + ls.model.J)

# Surface plot
p1 = plot(y_reg, x_reg, vec(y_hat),st=:surface,camera=(30,40),zlabel=L"\hat{y}(\mathbf{x})",legend=:none,size=(400, 370))
scatter3d!(X_data[:,2], X_data[:,1], vec(Y_data), label = "samples", markersize=2);
xlabel!(L"x_2", xguidefontsize=12)
ylabel!(L"x_1", yguidefontsize=12)
savefig(p1,"./MECH559/images/4_cont/ls_no_scaling.pdf")

# train models
ridge = 0.0
degree = 1
ls = surrogate_models.Surrogate(X=X_data, Y=Y_data, type="LS", lb=vec(lb), ub=vec(ub), r=ridge, d=degree, scale=true)
surrogate_models.train(ls.model);

# Predictions
y_hat = surrogate_models.predict(ls.model,X_reg)

println("scaled: ", cond(ls.model.B'*ls.model.B + ls.model.J))
println("scaled design matrix: ")
println(ls.model.B'*ls.model.B + ls.model.J)

# Surface plot
y_plot = vec(surrogate_models.scaling(X_reg[:,2:2],vec([lb[2]]),vec([ub[2]]),1))
x_plot = vec(surrogate_models.scaling(X_reg[:,1:1],vec([lb[1]]),vec([ub[1]]),1))
y_data = vec(surrogate_models.scaling(X_data[:,2:2],vec([lb[2]]),vec([ub[2]]),1))
x_data = vec(surrogate_models.scaling(X_data[:,1:1],vec([lb[1]]),vec([ub[1]]),1))

p1 = plot(x_plot, y_plot, vec(y_hat),st=:surface,camera=(30,40),zlabel=L"\hat{y}(\mathbf{x})",legend=:none,size=(400, 370))
scatter3d!(y_data, x_data, vec(Y_data), label = "samples", markersize=2);
xlabel!(L"x_2^\mathrm{scaled}", xguidefontsize=12)
ylabel!(L"x_1^\mathrm{scaled}", yguidefontsize=12)
savefig(p1,"./MECH559/images/4_cont/ls_scaling.pdf")