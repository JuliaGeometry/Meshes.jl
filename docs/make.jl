using Documenter, Meshes
using DocumenterTools: Themes

istravis = "TRAVIS" ∈ keys(ENV)

Themes.compile(
  joinpath(@__DIR__, "src/assets/light.scss"),
  joinpath(@__DIR__, "src/assets/themes/documenter-light.css")
)
Themes.compile(joinpath(@__DIR__, "src/assets/dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

makedocs(
  format=Documenter.HTML(
    assets=[
      "assets/favicon.ico",
      asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css)
    ],
    prettyurls=istravis,
    mathengine=KaTeX(
      Dict(
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
      )
    )
  ),
  sitename="Meshes.jl",
  authors="Júlio Hoffimann and contributors",
  pages=[
    "Home" => "index.md",
    "Reference guide" => [
      "Vectors" => "vectors.md",
      "Geometries" => ["geometries/primitives.md", "geometries/polytopes.md"],
      "Domains" => ["domains/sets.md", "domains/meshes.md", "domains/trajectories.md"],
      "Predicates" => "predicates.md",
      "Algorithms" => [
        "algorithms/sampling.md",
        "algorithms/partitioning.md",
        "algorithms/discretization.md",
        "algorithms/refinement.md",
        "algorithms/simplification.md",
        "algorithms/intersection.md",
        "algorithms/clipping.md",
        "algorithms/merging.md",
        "algorithms/winding.md",
        "algorithms/sideof.md",
        "algorithms/orientation.md",
        "algorithms/neighborsearch.md",
        "algorithms/boundingbox.md",
        "algorithms/hulls.md"
      ],
      "Transforms" => "transforms.md",
      "Visualization" => "visualization.md",
      "Input/Output" => "io.md"
    ],
    "Contributing" => ["contributing/guidelines.md"],
    "About" => ["License" => "about/license.md"],
    "Index" => "links.md"
  ]
)

repo = "github.com/JuliaGeometry/MeshesDocs.git"

withenv("GITHUB_REPOSITORY" => repo) do
  deploydocs(; repo, versions=["stable" => "v^", "dev" => "dev"])
end
