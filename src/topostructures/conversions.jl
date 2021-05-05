# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Base.convert(::Type{HalfEdgeStructure}, s::TopologicalStructure)
  HalfEdgeStructure(collect(elements(s)))
end

function Base.convert(::Type{FullStructure}, s::TopologicalStructure)
  # TODO: add all faces, not just the elements
  FullStructure(collect(elements(s)))
end
