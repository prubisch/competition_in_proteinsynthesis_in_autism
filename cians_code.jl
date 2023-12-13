using Distributions, Random, DifferentialEquations,Plots

# Set parameters
N = 2 # number of genes
l_mean = 1
l_sd = 1.
l_alpha = (l_mean/l_sd)^2 # parameters for gamma distribution
l_theta = l_sd^2/l_mean # parameters for gamma distribution
α_mean = 0.05
# α_sd = 0.2
β_mean = 5
β_sd = β_mean
β_alpha = (β_mean/β_sd)^2 # parameters for gamma distribution
β_theta = β_sd^2/β_mean # parameters for gamma distribution
γ_mean = .3
γ_sd = 0.0001
x0_mean = 1
x0_sd = 0.001
Y_max = 1*N
R = 1
# β_mean = (1/l_mean)*L_max*x0_mean*R

# calculate parameters
l = sort(rand(Gamma(l_alpha,l_theta),N))
α = α_mean ./ l .+ 0.001*randn(N)
# β = rand(Gamma(β_alpha,β_theta),N)
β = β_mean ./ l .+ 0.001*randn(N)
γ = γ_mean .+ γ_sd.*randn(N)
x0 = x0_mean .+ x0_sd.*randn(N)


# α[1] = 2*α[1]
paramvec = [α;β;γ;x0;l;Y_max;R]

# derivative function
function ddt(z,params,t)

    #unpack params
    params[1:N] = α
    params[N+1:2*N] = β
    params[2*N+1:3*N] = γ
    params[3*N+1:4*N] = x0
    params[4*N+1:5*N] = l
    params[end-1] = Y_max
    params[end] = R

    x_vec = z[1:N] # free mRNA concentrations
    xr_vec = z[N+1:2*N] # bound mRNA-ribosome concentrations
    y_vec = z[2*N+1:3*N] # protein concentrations

    # R = y_vec[1] # total concentration of ribosomes
    r = R - sum(xr_vec) # concentration of free ribosomes

    # free protein resource
    y_free = (Y_max-sum(y_vec))

    dxdt = (-α .* y_free .* r .* x_vec)  .+ (β .* y_free .* xr_vec )
    dxrdt = - dxdt
    dydt = (β .* y_free .* xr_vec ) .- (γ .* y_vec)

    return [dxdt; dxrdt; dydt]

end

# simulate
z_init = [x0; zeros(N); (Y_max/N)*ones(N)]
tspan = (0.0,100.0)
prob = ODEProblem(ddt,z_init,tspan,paramvec)

sol = solve(prob,KenCarp4())


# extract different elements
x_sol = zeros(length(sol),N)
xr_sol = zeros(length(sol),N)
y_sol = zeros(length(sol),N)
for i = 1:length(sol)
    x_sol[i,:] = sol.u[i][1:N]
    xr_sol[i,:] = sol.u[i][N+1:2*N]
    y_sol[i,:] = sol.u[i][2*N+1:3*N]
end
# r_sol = y_sol[:,1] .- sum(xr_sol,dims=2)[:]
r_sol = R .- sum(xr_sol,dims=2)[:]

######
# PLOTS
######
plot(sol.t,x_sol,ylim=(0,Inf),xlabel="Time",ylabel="Free mRNA")
plot(sol.t,xr_sol,ylim=(0,Inf),xlabel="Time",ylabel="Bound mRNA-ribosomes")
plot(sol.t,y_sol,ylim=(0,Inf),xlabel="Time",ylabel="Protein")
plot(sol.t,r_sol,ylim=(0,Inf),xlabel="Time",ylabel="Free ribosomes")

scatter(l,y_sol[end,:],ylim=(0,Inf),xlabel="mRNA length",ylabel="Protein concentration (a.u.)",yaxis=:linear)
scatter(α,y_sol[end,:],ylim=(0,Inf),xlabel="Initiation rate",ylabel="Protein concentration (a.u.)",yaxis=:linear)
scatter(β,y_sol[end,:],ylim=(0,Inf),xlabel="Elongation rate",ylabel="Protein concentration (a.u.)")

sum((y_sol[end,:]))/Y_max

