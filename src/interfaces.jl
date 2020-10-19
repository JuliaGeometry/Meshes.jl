# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# TODO: review these
"""
    coordinates(geometry)
Returns the vertices of a geometry.
"""
function coordinates(points::AbstractVector{<:Point})
    return points
end

"""
    faces(geometry)
Returns the face connections of a geometry. Is allowed to return lazy iterators!
Use `decompose(ConcreteFaceType, geometry)` to get `Vector{ConcreteFaceType}` with
`ConcreteFaceType` to be something like `TriangleFace`.
"""
function faces(f::AbstractVector{<:AbstractFace})
    return f
end

function faces(primitive, nvertices=nothing; kw...)
    # doesn't have any specific algorithm to generate faces
    # so will try to triangulate the coordinates!
    return nothing
end

"""
    Tesselation(primitive, nvertices)
For abstract geometries, when we generate
a mesh from them, we need to decide how fine grained we want to mesh them.
To transport this information to the various decompose methods, you can wrap it
in the Tesselation object e.g. like this:

```julia
sphere = Sphere(Point3f(0,0,0), 1.0f0)
m1 = mesh(sphere) # uses a default value for tesselation
m2 = mesh(Tesselation(sphere, 64)) # uses 64 for tesselation
length(coordinates(m1)) != length(coordinates(m2))
```
For grid based tesselation, you can also use a tuple.
"""
struct Tesselation{Dim,T,Primitive,NGrid}
    primitive::Primitive
    nvertices::NTuple{NGrid,Int}
end

function Tesselation(primitive::Primitive{Dim,T},
                     nvertices::NTuple{N,<:Integer}) where {Dim,T,N}
    return Tesselation{Dim,T,typeof(primitive),N}(primitive, Int.(nvertices))
end

Tesselation(primitive, nvertices::Integer) = Tesselation(primitive, (nvertices,))

# This is a bit lazy, I guess we should just refactor these methods
# to directly work on Tesselation - but this way it's backward compatible and less
# refactor work :D
nvertices(tesselation::Tesselation) = tesselation.nvertices
nvertices(tesselation::Tesselation{T,N,P,1}) where {T,N,P} = tesselation.nvertices[1]

function coordinates(tesselation::Tesselation)
    return coordinates(tesselation.primitive, nvertices(tesselation))
end

faces(tesselation::Tesselation) = faces(tesselation.primitive, nvertices(tesselation))

# Types that can be converted to a mesh via the functions below
const Meshable{Dim,T} = Union{Mesh{Dim,T},
                              Tesselation{Dim,T},
                              Polytope{Dim,T},
                              Primitive{Dim,T}}

"""
    decompose(T, meshable)

Decompose a `meshable` object (e.g. Polygon) into elements of type `T`.
"""
function decompose(::Type{T}, primitive) where {T}
    return collect_with_eltype(T, primitive)
end

# Specializations

function decompose(::Type{P}, primitive) where {P<:Point}
    convert.(P, coordinates(primitive))
end

function decompose(::Type{F}, primitive) where {F<:AbstractFace}
    f = faces(primitive)
    f === nothing && return nothing
    return collect_with_eltype(F, f)
end

# TODO: review these
function decompose(::Type{Point}, primitive::LineString{Dim,T}) where {Dim,T}
    return collect_with_eltype(Point{Dim,T}, coordinates(primitive))
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
