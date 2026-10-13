# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

getdifference(p₁::Union{Polygon,MultiPolygon}, p₂::Union{Polygon,MultiPolygon}) = _ghclip(p₁, p₂, :-)
