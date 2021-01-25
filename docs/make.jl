# Workaround for GR warnings
ENV["GKSwstype"] = "100"

using Documenter, Meshes
using DocumenterTools: Themes

istravis = "TRAVIS" ∈ keys(ENV)

Themes.compile(joinpath(@__DIR__,"src/assets/meshes-light.scss"), joinpath(@__DIR__,"src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__,"src/assets/meshes-dark.scss"), joinpath(@__DIR__,"src/assets/themes/documenter-dark.css"))

makedocs(
  format = Documenter.HTML(
    assets = ["assets/favicon.ico", asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css)],
    prettyurls = istravis,
    mathengine = KaTeX(Dict(
      :macros => Dict(
        "\\x" => "\\boldsymbol{x}",
        "\\z" => "\\boldsymbol{z}",
        "\\l" => "\\boldsymbol{\\lambda}",
        "\\c" => "\\boldsymbol{c}",
        "\\C" => "\\boldsymbol{C}",
        "\\g" => "\\boldsymbol{g}",
        "\\G" => "\\boldsymbol{G}",
        "\\f" => "\\boldsymbol{f}",
        "\\F" => "\\boldsymbol{F}",
        "\\R" => "\\mathbb{R}",
        "\\1" => "\\mathbb{1}"
      )
    ))
  ),
  sitename = "Meshes.jl",
  authors = "Júlio Hoffimann and contributors",
  pages = [
    "Home" => "index.md",
    "Reference guide" => [
      "Points" => "points.md",
      "Vectors" => "vectors.md",
      "Angles" => "angles.md",
      "Geometries" => [
        "geometries/primitives.md",
        "geometries/polytopes.md"
      ],
      "Meshes" => "meshes.md",
      "Algorithms" => [
        "algorithms/sampling.md",
        "algorithms/discretization.md",
        "algorithms/boundbox.md"
      ]
    ]
  ]
)

deploydocs(repo="github.com/JuliaGeometry/Meshes.jl.git")
