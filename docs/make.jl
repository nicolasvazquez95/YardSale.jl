# Add the package source directory to the LOAD_PATH
#push!(LOAD_PATH,"../src/")
# Load the packages
using Documenter
using YardSale
# Documenter can automatically generate documentation for a package.
makedocs(
    sitename = "YardSale.jl",
    format = Documenter.HTML(; assets = ["assets/custom.css"],
    prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [YardSale],
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "ODE systems" => "tutorials/ode_systems/ode_systems.md",
            "Monte Carlo" => [
                "EYSM Base" => "tutorials/mc_eysm_base/mc_eysm_base.md",
            ]
        ],
        "References" => [
            "ODE systems" => "references/ode.md",
            "Monte Carlo" => "references/mc.md",
            "Utils" => "references/utils.md",
        ]
    ],
    checkdocs = :none,

)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/nicolasvazquez95/YardSale.jl.git",
    branch = "gh-pages",
    target = "build",
)
