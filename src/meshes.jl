# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function decompose(::Type{T}, primitive) where {T}
    return collect_with_eltype(T, primitive)
end

function decompose(::Type{P}, primitive) where {P<:Point}
    return convert.(P, coordinates(primitive))
end

function decompose(::Type{F}, primitive) where {F<:AbstractFace}
    f = faces(primitive)
    f === nothing && return nothing
    return collect_with_eltype(F, f)
end

function collect_with_eltype(::Type{T}, iter) where {T}
    result = T[]
    for element in iter
        for telement in convert_simplex(T, element)
            push!(result, telement)
        end
    end
    return result
end

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
