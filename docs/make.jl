# Add the package source directory to the LOAD_PATH
#push!(LOAD_PATH,"../src/")
# Load the packages
using Documenter
using YardSale
# Documenter can automatically generate documentation for a package.
makedocs(
    sitename = "YardSale",
    format = Documenter.HTML(),
    modules = [YardSale],
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Tutorials" => [
            "ODE systems" => "tutorials/ode_systems/ode_systems.md",
        ]
    ],
    checkdocs = :none,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(;
    repo = "github.com/nicolasvazquez95/YardSale.jl.git",
    branch = "gh-pages",
    target = "build"
)
