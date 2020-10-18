const FaceMesh{Dim,T,Element} = Mesh{Dim,T,Element,<:FaceView{Element}}

coordinates(mesh::FaceMesh) = coordinates(getfield(mesh, :simplices))
faces(mesh::FaceMesh) = faces(getfield(mesh, :simplices))

"""
    TriangleMesh{Dim, T, PointType}

Abstract Mesh with triangle elements of eltype `T`.
"""
const TriangleMesh{Dim,T,PointType} = AbstractMesh{TriangleP{Dim,T,PointType}}

"""
    PlainMesh{Dim, T}

Triangle mesh (points + triangle faces)
"""
const PlainMesh{Dim,T} = TriangleMesh{Dim,T,Point{Dim,T}}

"""
    mesh(primitive::Meshable{N,T}; pointtype=Point{N,T}, facetype=TriangleFace)

Creates a mesh from `primitive`.
Uses the element types from the keyword arguments to create the attributes.
The attributes that have their type set to nothing are not added to the mesh.
Note, that this can be an `Int` or `Tuple{Int, Int}``, when the primitive is grid based.
It also only losely correlates to the number of vertices, depending on the algorithm used.
#TODO: find a better number here!
"""
function mesh(primitive::Meshable{N,T}; pointtype=Point{N,T}, facetype=TriangleFace) where {N,T}

    positions = decompose(pointtype, primitive)
    faces = decompose(facetype, primitive)
    if faces === nothing
        # try to triangulate
        faces = decompose(facetype, positions)
    end
    return Mesh(positions, faces)
end

function mesh(polygon::AbstractVector{P}; pointtype=P, facetype=TriangleFace) where {P<:AbstractPoint}

    return mesh(Polygon(polygon); pointtype=pointtype, facetype=facetype)
end

function mesh(polygon::Polytope{Dim,T}; pointtype=Point{Dim,T}, facetype=TriangleFace) where {Dim,T}
    faces = decompose(facetype, polygon)
    positions = decompose(pointtype, polygon)
    return Mesh(positions, faces)
end

function triangle_mesh(primitive::Meshable{N,T}; nvertices=nothing) where {N,T}
    if nvertices !== nothing
        @warn("nvertices argument deprecated. Wrap primitive in `Tesselation(primitive, nvertices)`")
        primitive = Tesselation(primitive, nvertices)
    end
    return mesh(primitive; pointtype=Point{N,T}, facetype=TriangleFace)
end

function triangle_mesh(points::AbstractVector{P}; nvertices=nothing) where {P<:AbstractPoint}
    triangle_mesh(Polygon(points), nvertices=nvertices)
end

"""
    volume(triangle)

Calculate the signed volume of one tetrahedron. Be sure the orientation of your
surface is right.
"""
function volume(triangle::Triangle) where {VT,FT}
    v1, v2, v3 = coordinates.(triangle)
    sig = sign(orthogonal_vector(v1, v2, v3) ⋅ v1)
    return sig * abs(v1 ⋅ (v2 × v3)) / 6
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
