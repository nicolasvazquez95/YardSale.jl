using YardSale
using Test, SafeTestsets

@time begin
   @safetestset "utils/kappa_beta" begin
        include("unit/utils/kappa_beta.jl")
   end

    @safetestset "ode/ode_solvers" begin
          include("unit/ode/ode_solvers.jl")
    end
end
