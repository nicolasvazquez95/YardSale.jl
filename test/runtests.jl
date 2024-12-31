using YardSale
using Test, SafeTestsets

@time begin
      # Utils
   @safetestset "utils/kappa_beta" begin
        include("unit/utils/kappa_beta.jl")
      end

   @safetestset "utils/mc_auxiliary_functions" begin
        include("unit/utils/mc_auxiliary_functions.jl")
      end
      # ODE
    @safetestset "ode/ode_solvers" begin
          include("unit/ode/ode_solvers.jl")
    end
      # MC
      ## Callbacks
    @safetestset "mc/callbacks" begin
        include("unit/mc/callbacks.jl")
      end
end
