using YardSale
using Test, SafeTestsets

@time begin
   @safetestset "utils/kappa_beta" begin
        include("unit/utils/kappa_beta.jl")
   end
end
