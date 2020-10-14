# TODO: review these
"""
    coordinates(geometry)
Returns the vertices of a geometry.
"""
function coordinates(points::AbstractVector{<:AbstractPoint})
    return points
end

"""
    faces(geometry)
Returns the face connections of a geometry. Is allowed to return lazy iterators!
Use `decompose(ConcreteFaceType, geometry)` to get `Vector{ConcreteFaceType}` with
`ConcreteFaceType` to be something like `TriangleFace{Int}`.
"""
function faces(f::AbstractVector{<:AbstractFace})
    return f
end

function normals(primitive, nvertices=nothing; kw...)
    # doesn't have any specific algorithm to generate normals
    # so will be generated from faces + positions
    # which we indicate by returning nothing!
    # Overload normals(primitive::YourPrimitive), to calcalute the normals
    # differently
    return nothing
end

function faces(primitive, nvertices=nothing; kw...)
    # doesn't have any specific algorithm to generate faces
    # so will try to triangulate the coordinates!
    return nothing
end

texturecoordinates(primitive, nvertices=nothing) = nothing

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
For grid based tesselation, you can also use a tuple:
```julia
rect = Rect2D(0, 0, 1, 1)
Tesselation(rect, (5, 5))
"""
struct Tesselation{Dim,T,Primitive,NGrid}
    primitive::Primitive
    nvertices::NTuple{NGrid,Int}
end

function Tesselation(primitive::GeometryPrimitive{Dim,T},
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
normals(tesselation::Tesselation) = normals(tesselation.primitive, nvertices(tesselation))
function texturecoordinates(tesselation::Tesselation)
    return texturecoordinates(tesselation.primitive, nvertices(tesselation))
end

# Types that can be converted to a mesh via the functions below
const Meshable{Dim,T} = Union{Mesh{Dim,T},
                              Tesselation{Dim,T},
                              AbstractPolygon{Dim,T},
                              GeometryPrimitive{Dim,T}}

"""
    decompose(T, meshable)

Decompose a `meshable` object (e.g. Polygon) into elements of type `T`

## Example

```julia
decompose(Point3, Rect3D())
```
"""
function decompose(::Type{T}, primitive) where {T}
    return collect_with_eltype(T, primitive)
end

# Specializations

function decompose(::Type{P}, primitive) where {P<:AbstractPoint}
    convert.(P, metafree(coordinates(primitive)))
end

function decompose(::Type{F}, primitive) where {F<:AbstractFace}
    f = faces(primitive)
    f === nothing && return nothing
    return collect_with_eltype(F, f)
end

# TODO: review these
function decompose(::Type{Point}, primitive::LineString{Dim,T}) where {Dim,T}
    return collect_with_eltype(Point{Dim,T}, metafree(coordinates(primitive)))
end

# TODO: review these
struct UV{T} end
UV(::Type{T}) where {T} = UV{T}()
UV() = UV(Vec2f)

struct UVW{T} end
UVW(::Type{T}) where {T} = UVW{T}()
UVW() = UVW(Vec3f)

struct Normal{T} end
Normal(::Type{T}) where {T} = Normal{T}()
Normal() = Normal(Vec3f)

decompose_uv(primitive) = decompose(UV(), primitive)
decompose_uvw(primitive) = decompose(UVW(), primitive)
decompose_normals(primitive) = decompose(Normal(), primitive)

function decompose(NT::Normal{T}, primitive) where {T}
    n = normals(primitive)
    if n === nothing
        return collect_with_eltype(T, normals(coordinates(primitive), faces(primitive)))
    end
    return collect_with_eltype(T, n)
end

function decompose(UVT::Union{UV{T},UVW{T}}, primitive::Meshable{Dim,V}) where {Dim,T,V}
    # This is the fallback for texture coordinates if a primitive doesn't overload them
    # We just take the positions and normalize them
    uv = texturecoordinates(primitive)
    if uv === nothing
        # If the primitive doesn't even have coordinates, we're out of options and return
        # nothing, indicating that texturecoordinates aren't implemented
        positions = decompose(Point{Dim,V}, primitive)
        positions === nothing && return nothing
        # Let this overlord do the work
        return decompose(UVT, positions)
    end
    return collect_with_eltype(T, uv)
end

function decompose(UVT::Union{UV{T},UVW{T}}, positions::AbstractVector{<:AbstractPoint}) where {T}
    N = length(T)
    bb = boundingbox(positions)
    return map(positions) do p
        return (coordinates(p) - minimum(bb)) ./ widths(bb)
    end
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
