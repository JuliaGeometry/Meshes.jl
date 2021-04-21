# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(mesh::Mesh)
  seriescolor --> :black
  fillcolor --> :gray90
  fill --> true

  for element in elements(mesh)
    @series begin
      primary --> false
      element
    end
  end
end
