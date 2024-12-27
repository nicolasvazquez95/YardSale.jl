# Monte Carlo Simulations
```@autodocs
Modules = [YardSale]
Pages = ["mc/mc_simulator.jl"]
Order   = [:function, :type]
```

# Callbacks

In the context of Monte Carlo simulations, a callback is a function that is called at the end of each iteration of the simulation. This can be used to monitor the progress of the simulation with any desired metrics. The callbacks functions take the state of the simulation as an argument, i.e. `f(w)` where `w` is the wealth distribution at the current iteration. It is easy to define custom callbacks by defining a function that takes the wealth distribution as an argument and returns something (a number, a vector, etc.). The function can then be passed to the `mc_simulator` function as a callback.

Some predefined callbacks are available in the `YardSale` module, which are commonly used in the study of wealth distributions. These are:

```@autodocs
Modules = [YardSale]
Pages = ["mc/callbacks.jl"]
Order   = [:function, :type]
```