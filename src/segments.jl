# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    segments(chain)

Segments connecting vertices of the `chain`.
"""
function segments(c::Chain)
  v = vertices(c)
  n = length(v) - !isclosed(c)
  @inbounds (Segment(v[i], v[i + 1]) for i in 1:n)
end

"""
    segments(mesh)

Segments connecting vertices of the `mesh`.
"""
segments(m::Mesh) = _segments(topology(m), m)

# segments are just faces of rank 1
_segments(::Topology, m) = faces(m, 1)

# try harder in case of simple topology
# because it doesn't implement faces
function _segments(::SimpleTopology, m)
  try
    # attempt to convert to a topology that
    # provides an efficient implementation
    m′ = topoconvert(HalfEdgeTopology, m)
    _segments(topology(m′), m′)
  catch
    error("not implemented")
  end
end
