# Import necessary packages

```julia
using YardSale
using Plots, Measures, LaTeXStrings, ColorSchemes
```

# Description of the notebook
In this notebook we'll explore the behavior of the Extended Yard-Sale Model (EYSM), simulating the (stochastic) evolution of an economical systems of $N$ agents. These agents are characterized by their wealth $w_i, \;i=1,...,N$. The rules of the model have been described in the `ode_systems.ipynb` notebook. 
Here, we add the stochasticity to the model, in the framework of Monte Carlo simulations. With this addition, we expect to see not only the phase transition of the model, but also the fluctuations around the critical point. 

Remember that the exchange rule of the EYSM is given by:
```math
\Delta w_{ij} = \chi \Delta t \bigg(\frac{W}{N} - w_i\bigg) + \eta_{ij} \sqrt{\gamma \Delta t} \min(w_i,w_j),
```
where ``\chi`` is the redistribution rate, ``\sqrt{\gamma \Delta t}`` represents the fraction ``f`` of the wealth that is exchanged between agents ``i`` and ``j`` at each time step, and ``\eta_{ij}`` is a random number that can be ``\pm 1``. This stochastic variable has the properties

```math
\langle \eta_{ij} \rangle = \zeta N \sqrt{\frac{\Delta t}{\gamma}} \frac{(w_i - w_j)}{W}, \qquad \langle \eta_{ij}^2 \rangle = 1.
```
The parameter ``\zeta`` represents the Wealth-Attained Advantage (WAA) parameter, which is a bias that favors the richer agent in the exchange.
Finally, the wealth of the agents is updated as
```math
w_i \to w_i + \Delta w_{ij}, \qquad w_j \to w_j - \Delta w_{ij}.
```

In this notebook, we''ll explore the behavior of the EYSM model using Monte Carlo simulations. We''ll study the dynamics and the phase transition of the model measuring the order parameter and its fluctuations at the steady state.

# 1. Time evolution of the wealth distribution
For a small system, we can visualize the evolution of the wealth of each agent over time.


```julia
# Parameters of the system
N = 16
W_N = 1.0

# Temperature of the system: χ/ζ
chi = 1.0
zeta = 1.0

# The fraction of exchange determines in some way the intensity of the noise,
# Larger values of f will speed up the system, but it will also make the noise
# more intense, an the dynamics will be less smooth. 
# I found that f < 0.1 is a rule of thumb.
f = 0.01

# Steps to run the simulation
steps = 10000 # MC steps
# Seed
seed = 42 # The answer to the ultimate question of life, the universe and everything

# Optional
save_every = 1 
;
```


```julia
# Run 
@time w_t = EYSM_base_full(N, W_N, chi, zeta, f, steps, seed; save_every=save_every)

# Sort the columns according to the final state
w_t = w_t[:,sortperm(w_t[end,:])]
```
      0.026488 seconds (350.02 k allocations: 14.344 MiB)
    10001×16 Matrix{Float64}:
     1.0       1.0       1.0       1.0       …  1.0       1.0       1.0
     1.00988   1.00011   0.990101  1.01         0.990199  1.0       1.01979
     1.00988   0.999997  0.970359  0.989991     0.990195  0.999996  1.02948
     1.02985   0.980294  0.980077  0.999894     0.970505  0.989931  1.04963
     1.01964   0.970732  0.979959  0.999896     0.970441  0.989934  1.04964
     1.00953   1.00018   0.979958  1.02938   …  0.960755  0.999838  1.02978
     0.999528  0.981201  0.989694  1.02911      0.97039   1.01993   1.02015
     0.999133  0.971524  1.00831   1.02907      0.960703  1.0204    1.02003
     0.989324  0.99047   0.998984  1.03846      0.97033   1.03038   1.03996
     0.979914  0.990466  0.998689  1.02873      0.950995  1.04025   1.03979
     0.969973  0.990473  0.979881  1.01851   …  0.950961  1.0313    1.03885
     0.970002  0.98056   0.960972  1.00877      0.970246  1.04145   1.02892
     0.969923  0.970407  0.98012   1.00876      0.970233  1.07141   1.01859
     ⋮                                       ⋱                      ⋮
     0.440844  0.552208  0.643806  0.75115      1.56723   1.76156   1.9168
     0.445601  0.547458  0.646276  0.743793  …  1.56723   1.73571   1.89882
     0.450387  0.559005  0.652935  0.751363     1.54405   1.74003   1.89271
     0.451083  0.582863  0.659337  0.73784      1.56042   1.74721   1.8839
     0.448141  0.582775  0.652858  0.730516     1.59655   1.74695   1.8786
     0.45294   0.571665  0.653567  0.736467     1.57927   1.76105   1.8678
     0.452921  0.555469  0.647222  0.757503  …  1.57921   1.75976   1.85909
     0.448712  0.555443  0.648287  0.757467     1.5646    1.78162   1.87482
     0.440498  0.561309  0.642063  0.757512     1.57233   1.78772   1.86648
     0.436412  0.555931  0.642407  0.757888     1.59429   1.76861   1.84066
     0.428482  0.567727  0.642484  0.76592      1.58645   1.77934   1.84088
     0.42847   0.579662  0.612112  0.765899  …  1.58641   1.77612   1.83163




```julia
p1 = plot(
    ylabel = L"x_i",
    xlabel = L"t"*" [MC steps]",
    fontfamily = "Computer Modern",
    legend_title = "MC - Wealth distribution",
    xguidefontsize = 17,
    yguidefontsize = 17,
    legendfontsize = 15,
    legendtitlefontsize = 17,
    xtickfontsize = 17,
    ytickfontsize = 17,
    palette = :tol_nightfall,
    size=(640,480),
    fmt=:png,
    ylims = (0, 0.23),
    xlims = (0,:auto),
    minorticks = true,
    right_margin = 8mm,
    legend = :topright,
)    

for i in 1:N
    plot!(p1,w_t[:,i]/N, label = "", lw = 2)
end

display(p1)
```
![MC time evolution](output_8_0.png)


# 2. Phase transition properties

Now, we would like to study the phase transition of the model. For that, we need to make some statistics over a certain number of realizations of the model. Then we can calculate the order parameter and its fluctuations.
In Julia, we can use the `Threads.@threads` macro to parallelize the independent simulations. This will speed up the calculations. Then the order parameter and its fluctuations can be calculated with the `Statistics.jl` package.


```julia
using ProgressMeter
```
```julia
# Number of runs
N = 64
W_N = 1.0f0
zeta = 1.0f0
chi = 0.0f0:0.04f0:2.0f0
steps = 30000 * N # MC steps: 30000 * 128 = 3.84e6
seed = 42
# Run
sims = zeros(typeof(W_N), length(chi), (steps ÷ N) + 1, N)

# Safe lock 
lk = ReentrantLock()
# ProgressMeter
p = Progress(length(chi), showspeed=true)
Threads.@threads for i in 1:length(chi)
    w_t = EYSM_base_full(N, W_N, chi[i], zeta, f, steps, seed)
    lock(lk) do 
        sims[i,:,:] = w_t
        next!(p)
    end
end

# Save with JLD2
#using JLD2
#save("sims.jld2", "sims", sims)
# Load with JLD2
#sims = load("sims.jld2")["sims"]
``` 
Now at each temperature, we can calculate the mean value of the order parameter and its fluctuations using simple statistics.


```julia
# Each row is x1(T,t) for a given T=χ/ζ
x1 = dropdims(maximum(sims,dims=3),dims=3)
# Calculate mean and variance for all temperatures, in the steady state
using Statistics
start = 20000
mean_x1 = mean(x1[:,start:end]/N,dims=2)
var_x1 = var(x1[:,start:end]/N,dims=2)
;
```

```julia
p1 = plot(
    ylabel = L"x_1",
    xlabel = L"T",
    fontfamily = "Computer Modern",
    xguidefontsize = 17,
    yguidefontsize = 17,
    legendfontsize = 15,
    legendtitlefontsize = 17,
    xtickfontsize = 17,
    ytickfontsize = 17,
    minorticks = true,
    palette = :julia,
    size=(640,480),
    y_foreground_color_text = ColorSchemes.julia[1],
    y_foreground_color_axis = ColorSchemes.julia[1]
)

scatter!(p1, chi /zeta, mean_x1, yerr=sqrt.(var_x1),label="MC, " *L"N="*"$N",ms=5)
p2 = twinx(p1)
scatter!(
    p2,
    chi /zeta, 
    1e4*var_x1, 
    label="",
    ms=5,
    palette=:julia,
    color=2,
    minorticks=true,
    ylabel=L"\sigma^2(x_1)\,[\,\times 10^{-4}]",
    yguidefontsize = 17,
    ytickfontsize = 17,
    y_foreground_color_text = ColorSchemes.julia[2],
    y_foreground_color_axis = ColorSchemes.julia[2]
    )
```
![MC phase transition](output_14_0.png)

As it can be seen, the two curves are consistent, and the order parameter shows a clear phase transition. The fluctuations of the order parameter are also shown, and they are maximum around the critical point.