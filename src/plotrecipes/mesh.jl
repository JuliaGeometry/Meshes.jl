# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(mesh::Mesh)
  seriescolor --> :black
  legend --> false

  for elm in elements(mesh)
    @series begin
      primary --> false
      elm
    end
  end
end
