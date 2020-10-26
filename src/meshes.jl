# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    mesh(geometry::Meshable{N,T}; facetype=TriangleFace)

Creates a mesh from `geometry`.
"""
function mesh(geometry::Geometry{Dim,T}; facetype=TriangleFace) where {Dim,T}
    points = decompose(Point{Dim,T}, geometry)
    faces  = decompose(facetype, geometry)
    if faces === nothing
        # try to triangulate
        faces = decompose(facetype, points)
    end
    return Mesh(points, faces)
end
