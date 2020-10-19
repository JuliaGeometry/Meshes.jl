# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const FaceMesh{Dim,T,Element} = Mesh{Dim,T,Element,<:FaceView{Element}}
coordinates(mesh::FaceMesh) = coordinates(getfield(mesh, :simplices))
faces(mesh::FaceMesh) = faces(getfield(mesh, :simplices))

"""
    mesh(geometry::Meshable{N,T}; facetype=TriangleFace)

Creates a mesh from `geometry`.
"""
function mesh(geometry::Meshable{N,T}; facetype=TriangleFace) where {N,T}
    positions = decompose(Point{N,T}, geometry)
    faces = decompose(facetype, geometry)
    if faces === nothing
        # try to triangulate
        faces = decompose(facetype, positions)
    end
    return Mesh(positions, faces)
end

function mesh(polygon::AbstractVector{P}; facetype=TriangleFace) where {P<:Point}
    return mesh(Polygon(polygon); facetype=facetype)
end

function mesh(polytope::Polytope{Dim,T}; facetype=TriangleFace) where {Dim,T}
    faces = decompose(facetype, polytope)
    positions = decompose(Point{Dim,T}, polytope)
    return Mesh(positions, faces)
end

function triangle_mesh(geometry::Meshable)
    return mesh(geometry; facetype=TriangleFace)
end

function triangle_mesh(points::AbstractVector{P}) where {P<:Point}
    return triangle_mesh(Polygon(points))
end

"""
    volume(mesh)

Calculate the signed volume of all tetrahedra. Be sure the orientation of your
surface is right.
"""
function volume(mesh::Mesh) where {VT,FT}
    return sum(volume, mesh)
end

function Base.merge(meshes::AbstractVector{<:Mesh})
    return if isempty(meshes)
        error("No meshes to merge")
    elseif length(meshes) == 1
        return meshes[1]
    else
        m1 = meshes[1]
        ps = copy(coordinates(m1))
        fs = copy(faces(m1))
        for mesh in Iterators.drop(meshes, 1)
            append!(fs, map(f -> f .+ length(ps), faces(mesh)))
            append!(ps, coordinates(mesh))
        end
        return Mesh(ps, fs)
    end
end
