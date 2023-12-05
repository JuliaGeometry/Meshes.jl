# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

intersection(f, chain₁::Chain, chain₂::Chain) =
  intersection(f, GeometrySet(segments(chain₁)), GeometrySet(segments(chain₂)))
