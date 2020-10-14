const FaceMesh{Dim,T,Element} = Mesh{Dim,T,Element,<:FaceView{Element}}

coordinates(mesh::FaceMesh) = coordinates(getfield(mesh, :simplices))
faces(mesh::FaceMesh) = faces(getfield(mesh, :simplices))

function texturecoordinates(mesh::AbstractMesh)
    hasproperty(mesh, :uv) && return mesh.uv
    hasproperty(mesh, :uvw) && return mesh.uvw
    return nothing
end

function normals(mesh::AbstractMesh)
    hasproperty(mesh, :normals) && return mesh.normals
    return nothing
end

const GLTriangleElement = Triangle{3,Float32}
const GLTriangleFace = TriangleFace{GLIndex}
const PointWithUV{Dim,T} = PointMeta{Dim,T,Point{Dim,T},(:uv,),Tuple{Vec{2,T}}}
const PointWithNormal{Dim,T} = PointMeta{Dim,T,Point{Dim,T},(:normals,),Tuple{Vec{3,T}}}
const PointWithUVNormal{Dim,T} = PointMeta{Dim,T,Point{Dim,T},(:normals, :uv),
                                           Tuple{Vec{3,T},Vec{2,T}}}
const PointWithUVWNormal{Dim,T} = PointMeta{Dim,T,Point{Dim,T},(:normals, :uvw),
                                            Tuple{Vec{3,T},Vec{3,T}}}

"""
    TriangleMesh{Dim, T, PointType}

Abstract Mesh with triangle elements of eltype `T`.
"""
const TriangleMesh{Dim,T,PointType} = AbstractMesh{TriangleP{Dim,T,PointType}}

"""
    PlainMesh{Dim, T}

Triangle mesh with no meta information (just points + triangle faces)
"""
const PlainMesh{Dim,T} = TriangleMesh{Dim,T,Point{Dim,T}}
const GLPlainMesh{Dim} = PlainMesh{Dim,Float32}
const GLPlainMesh2D = GLPlainMesh{2}
const GLPlainMesh3D = GLPlainMesh{3}

"""
    UVMesh{Dim, T}

PlainMesh with texture coordinates meta at each point.
`uvmesh.uv isa AbstractVector{Vec2f}`
"""
const UVMesh{Dim,T} = TriangleMesh{Dim,T,PointWithUV{Dim,T}}
const GLUVMesh{Dim} = UVMesh{Dim,Float32}
const GLUVMesh2D = UVMesh{2}
const GLUVMesh3D = UVMesh{3}

"""
    NormalMesh{Dim, T}

PlainMesh with normals meta at each point.
`normalmesh.normals isa AbstractVector{Vec3f}`
"""
const NormalMesh{Dim,T} = TriangleMesh{Dim,T,PointWithNormal{Dim,T}}
const GLNormalMesh{Dim} = NormalMesh{Dim,Float32}
const GLNormalMesh2D = GLNormalMesh{2}
const GLNormalMesh3D = GLNormalMesh{3}

"""
    NormalUVMesh{Dim, T}

PlainMesh with normals and uv meta at each point.
`normalmesh.normals isa AbstractVector{Vec3f}`
`normalmesh.uv isa AbstractVector{Vec2f}`
"""
const NormalUVMesh{Dim,T} = TriangleMesh{Dim,T,PointWithUVNormal{Dim,T}}
const GLNormalUVMesh{Dim} = NormalUVMesh{Dim,Float32}
const GLNormalUVMesh2D = GLNormalUVMesh{2}
const GLNormalUVMesh3D = GLNormalUVMesh{3}

"""
    NormalUVWMesh{Dim, T}

PlainMesh with normals and uvw (texture coordinates in 3D) meta at each point.
`normalmesh.normals isa AbstractVector{Vec3f}`
`normalmesh.uvw isa AbstractVector{Vec3f}`
"""
const NormalUVWMesh{Dim,T} = TriangleMesh{Dim,T,PointWithUVWNormal{Dim,T}}
const GLNormalUVWMesh{Dim} = NormalUVWMesh{Dim,Float32}
const GLNormalUVWMesh2D = GLNormalUVWMesh{2}
const GLNormalUVWMesh3D = GLNormalUVWMesh{3}

"""
    mesh(primitive::Meshable{N,T};
         pointtype=Point{N,T}, facetype=GLTriangle,
         uvtype=nothing, normaltype=nothing)

Creates a mesh from `primitive`.
Uses the element types from the keyword arguments to create the attributes.
The attributes that have their type set to nothing are not added to the mesh.
Note, that this can be an `Int` or `Tuple{Int, Int}``, when the primitive is grid based.
It also only losely correlates to the number of vertices, depending on the algorithm used.
#TODO: find a better number here!
"""
function mesh(primitive::Meshable{N,T}; pointtype=Point{N,T}, facetype=GLTriangleFace,
              uv=nothing, normaltype=nothing) where {N,T}

    positions = decompose(pointtype, primitive)
    faces = decompose(facetype, primitive)
    # If faces returns nothing for primitive, we try to triangulate!
    if faces === nothing
        # triangulation.jl
        faces = decompose(facetype, positions)
    end

    # We want to preserve any existing attributes!
    attrs = attributes(primitive)
    # Make sure this doesn't contain position, we'll add position explicitely via meta!
    delete!(attrs, :position)

    if uv !== nothing
        # this may overwrite an existing :uv, but will only create a copy
        # if it has a different eltype, otherwise it should replace it
        # with exactly the same instance - which is what we want here
        attrs[:uv] = decompose(UV(uv), primitive)
    end

    if normaltype !== nothing
        primitive_normals = normals(primitive)
        if primitive_normals !== nothing
            attrs[:normals] = convert.(normaltype, primitive_normals)
        else
            # Normals not implemented for primitive, so we calculate them!
            n = normals(positions, faces; normaltype=normaltype)
            if n !== nothing # ok jeez, this is a 2d mesh which cant have normals
                attrs[:normals] = n
            end
        end
    end
    return Mesh(meta(positions; attrs...), faces)
end

"""
    mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace,
         normaltype=nothing)
Polygon triangluation!
"""
function mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace,
              normaltype=nothing) where {P<:AbstractPoint{2}}

    return mesh(Polygon(polygon); pointtype=pointtype, facetype=facetype,
                normaltype=normaltype)
end

function mesh(polygon::AbstractPolygon{Dim,T}; pointtype=Point{Dim,T},
              facetype=GLTriangleFace, normaltype=nothing) where {Dim,T}

    faces = decompose(facetype, polygon)
    positions = decompose(pointtype, polygon)

    if normaltype !== nothing
        n = normals(positions, faces; normaltype=normaltype)
        positions = meta(positions; normals=n)
    end
    return Mesh(positions, faces)
end

function triangle_mesh(primitive::Meshable{N,T}; nvertices=nothing) where {N,T}
    if nvertices !== nothing
        @warn("nvertices argument deprecated. Wrap primitive in `Tesselation(primitive, nvertices)`")
        primitive = Tesselation(primitive, nvertices)
    end
    return mesh(primitive; pointtype=Point{N,T}, facetype=GLTriangleFace)
end

function triangle_mesh(points::AbstractVector{P}; nvertices=nothing) where {P<:AbstractPoint}
    triangle_mesh(Polygon(points), nvertices=nvertices)
end

function uv_mesh(primitive::Meshable{N,T}) where {N,T}
    mesh(primitive; pointtype=Point{N,T}, uv=Vec{2,T}, facetype=GLTriangleFace)
end

function uv_normal_mesh(primitive::Meshable{N,T}) where {N,T}
    mesh(primitive; pointtype=Point{N,T}, uv=Vec{2,T}, normaltype=Vec{3,T}, facetype=GLTriangleFace)
end

function normal_mesh(points::AbstractVector{<:AbstractPoint},
                     faces::AbstractVector{<:AbstractFace})
    _points = decompose(Point3f, points)
    _faces = decompose(GLTriangleFace, faces)
    return Mesh(meta(_points; normals=normals(_points, _faces)), _faces)
end

function normal_mesh(primitive::Meshable{N}; nvertices=nothing) where {N}
    if nvertices !== nothing
        @warn("nvertices argument deprecated. Wrap primitive in `Tesselation(primitive, nvertices)`")
        primitive = Tesselation(primitive, nvertices)
    end
    return mesh(primitive; pointtype=Point{N,Float32}, normaltype=Vec3f, facetype=GLTriangleFace)
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

"""
    pointmeta(mesh::Mesh; meta_data...)

Attaches metadata to the coordinates of a mesh
"""
function pointmeta(mesh::Mesh; meta_data...)
    points = coordinates(mesh)
    attr = attributes(points)
    delete!(attr, :position) # position == metafree(points)
    # delete overlapping attributes so we can replace with `meta_data`
    foreach(k -> delete!(attr, k), keys(meta_data))
    return Mesh(meta(metafree(points); attr..., meta_data...), faces(mesh))
end

function pointmeta(mesh::Mesh, uv::UV)
    return pointmeta(mesh; uv=decompose(uv, mesh))
end

function pointmeta(mesh::Mesh, normal::Normal)
    return pointmeta(mesh; normals=decompose(normal, mesh))
end

"""
    pop_pointmeta(mesh::Mesh, property::Symbol)
Remove `property` from point metadata.
Returns the new mesh, and the property!
"""
function pop_pointmeta(mesh::Mesh, property::Symbol)
    points = coordinates(mesh)
    attr = attributes(points)
    delete!(attr, :position) # position == metafree(points)
    # delete overlapping attributes so we can replace with `meta_data`
    m = pop!(attr, property)
    return Mesh(meta(metafree(points); attr...), faces(mesh)), m
end

"""
    facemeta(mesh::Mesh; meta_data...)

Attaches metadata to the faces of a mesh
"""
function facemeta(mesh::Mesh; meta_data...)
    return Mesh(coordinates(mesh), meta(faces(mesh); meta_data...))
end

function attributes(hasmeta::Mesh)
    return Dict{Symbol,Any}((name => getproperty(hasmeta, name)
                             for name in propertynames(hasmeta)))
end
