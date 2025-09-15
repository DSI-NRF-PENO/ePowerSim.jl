using Documenter
using ePowerSim

push!(LOAD_PATH,"../src/")

DocMeta.setdocmeta!(ePowerSim, :DocTestSetup, :(using ePowerSim); recursive=true)

makedocs(;
    modules=[ePowerSim],
    authors="dsi.nrf.peno@gmail.com",
    sitename="ePowerSim.jl Documentation",
    format=Documenter.HTML(;
        edit_link="main",
        size_threshold_warn = 800 * 1024, # 819 KiB
        size_threshold = 1000 * 1024,    # 1 MiB
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "Components" => "components.md",
            "Modelling Concepts" => "modelling-concepts.md",
            "Static Modelling and Simulation" => "static-modelling-and-simulation.md",
            "Dynamic Modelling and Simulation" => "dynamic-modelling-and-simulation.md",
        ],
        "Api" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/DSI-NRF-PENO/ePowerSim.jl",
    devbranch="main",
)
