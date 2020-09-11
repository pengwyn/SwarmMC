using Documenter
using SwarmMC

makedocs(
    sitename = "SwarmMC",
    format = Documenter.HTML(),
    modules = [SwarmMC]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/pengwyn/SwarmMC.git"
)
