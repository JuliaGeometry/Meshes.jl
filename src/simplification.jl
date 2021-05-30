# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimplificationMethod

A method for simplifying geometries.
"""
abstract type SimplificationMethod end

"""
    simplify(geometry, method)

Simplify `geometry` with given `method`.
"""
function simplify end

function simplify(multi::Multi, method::SimplificationMethod)
  Multi([simplify(geom, method) for geom in multi])
end

function simplify(polygon::Polygon, method::SimplificationMethod)
  c = [simplify(chain, method) for chain in chains(polygon)]
  PolyArea(c[1], c[2:end])
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("simplification/douglaspeucker.jl")
