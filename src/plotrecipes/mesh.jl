# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(mesh::Mesh)
  seriescolor --> :black
  fill --> true
  fillcolor --> :gray90
  primary --> false

  for elm in elements(mesh)
    @series begin
      elm
    end
  end
end
