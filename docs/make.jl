using Documenter
using YardSale

makedocs(
    sitename = "YardSale",
    format = Documenter.HTML(),
    modules = [YardSale]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
