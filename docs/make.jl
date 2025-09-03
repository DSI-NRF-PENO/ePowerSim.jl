using ePowerSim
using Documenter

DocMeta.setdocmeta!(ePowerSim, :DocTestSetup, :(using ePowerSim); recursive=true)

makedocs(;
    modules=[ePowerSim],
    authors="Adedayo Ademola Yusuff dsi.nrf.peno@gmail.com",
    sitename="ePowerSim.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Adedayo Ademola Yusuff/ePowerSim.jl",
    devbranch="main",
)
