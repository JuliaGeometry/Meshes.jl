# Helper functions that works around the fact, that there is no generic
# Table interface for this functionality. Once this is in e.g. Tables.jl,
# it should be removed from Meshes

"""
    attributes(hasmeta)
Returns all attributes of meta as a Dict{Symbol, Any}.
Needs to be overloaded, and returns empty dict for non overloaded types!
Gets overloaded by default for all Meta types.
"""
function attributes(hasmeta)
    return Dict{Symbol,Any}()
end

function attributes(hasmeta::StructArray)
    return Dict{Symbol,Any}((name => getproperty(hasmeta, name)
                             for name in propertynames(hasmeta)))
end

"""
    getcolumns(t, colnames::Symbol...)

Gets a column from any Array like (Table/AbstractArray).
For AbstractVectors, a column will be the field names of the element type.
"""
function getcolumns(tablelike, colnames::Symbol...)
    return getproperty.((tablelike,), colnames)
end

getcolumn(t, colname::Symbol) = getcolumns(t, colname)[1]

"""
    MetaType(::Type{T})

Returns the Meta Type corresponding to `T`
E.g:
```julia
MetaType(Point) == PointMeta
```
"""
MetaType(::Type{T}) where {T} = error("No Meta Type for $T")

"""
    MetaFree(::Type{T})

Returns the original type containing no metadata for `T`
E.g:
```julia
MetaFree(PointMeta) == Point
```
"""
MetaFree(::Type{T}) where {T} = error("No meta free Type for $T")

"""
    meta(x::MetaObject)

Returns the metadata of `x`
"""
meta(x::T) where {T} = error("$T has no meta!")

metafree(x::T) where {T} = x

macro meta_type(name, mainfield, supertype, params...)
    MetaName = Symbol("$(name)Meta")
    field = QuoteNode(mainfield)
    NoParams = Symbol("$(MetaName)NoParams")

    params_sym = map(params) do param
        param isa Symbol && return param
        param isa Expr && param.head == :(<:) && return param.args[1]
        return error("Unsupported type parameter: $(param)")
    end

    expr = quote
        struct $MetaName{$(params...),Typ<:$supertype{$(params_sym...)},Names,Types} <:
               $supertype{$(params_sym...)}
            main::Typ
            meta::NamedTuple{Names,Types}
        end

        const $NoParams{Typ,Names,Types} = $MetaName{$(params_sym...),Typ,Names,
                                                     Types} where {$(params_sym...)}

        function Base.getproperty(x::$MetaName{$(params_sym...),Typ,Names,Types},
                                  field::Symbol) where {$(params...),Typ,Names,Types}
            field === $field && return getfield(x, :main)
            field === :main && return getfield(x, :main)
            Base.sym_in(field, Names) && return getfield(getfield(x, :meta), field)
            return error("Field $field not part of Element")
        end

        function Meshes.MetaType(XX::Type{<:$supertype{$(params_sym...)} where {$(params...)}})
            return $MetaName
        end

        function Meshes.MetaType(ST::Type{<:$supertype{$(params_sym...)}},
                                         ::Type{NamedTuple{Names,Types}}) where {$(params...),
                                                                                 Names,
                                                                                 Types}
            return $MetaName{$(params_sym...),ST,Names,Types}
        end

        Meshes.MetaFree(::Type{<:$MetaName{Typ}}) where {Typ} = Typ
        Meshes.MetaFree(::Type{<:$MetaName}) = $name
        Meshes.metafree(x::$MetaName) = getfield(x, :main)
        Meshes.metafree(x::AbstractVector{<:$MetaName}) = getproperty(x, $field)
        Meshes.meta(x::$MetaName) = getfield(x, :meta)
        Meshes.meta(x::AbstractVector{<:$MetaName}) = getproperty(x, :meta)

        function Meshes.meta(main::$supertype{$(params_sym...)};
                                     meta...) where {$(params...)}
            isempty(meta) && return elements # no meta to add!
            return $MetaName(main; meta...)
        end

        function Meshes.meta(elements::AbstractVector{XX};
                                     meta...) where {XX<:$supertype{$(params_sym...)}} where {$(params...)}
            isempty(meta) && return elements # no meta to add!
            n = length(elements)
            for (k, v) in meta
                if v isa AbstractVector
                    mn = length(v)
                    mn != n && error("Metadata array needs to have same length as data.
                    Found $(n) data items, and $mn metadata items")
                else
                    error("Metadata needs to be an array with the same length as data items. Found: $(typeof(v))")
                end
            end
            nt = values(meta)
            # get the first element to get the per element named tuple type
            ElementNT = typeof(map(first, nt))

            return StructArray{MetaType(XX, ElementNT)}(($(mainfield)=elements, nt...))
        end

        function Meshes.attributes(hasmeta::$MetaName)
            return Dict{Symbol,Any}((name => getproperty(hasmeta, name)
                                     for name in propertynames(hasmeta)))
        end

        function (MT::Type{<:$MetaName})(args...; meta...)
            nt = values(meta)
            obj = MetaFree(MT)(args...)
            return MT(obj, nt)
        end

        function (MT::Type{<:$MetaName})(main::$(name); meta...)
            nt = values(meta)
            return MT(main, nt)
        end

        function Base.propertynames(::$MetaName{$(params_sym...),Typ,Names,Types}) where {$(params...),
                                                                                          Typ,
                                                                                          Names,
                                                                                          Types}
            return ($field, Names...)
        end

        function StructArrays.staticschema(::Type{$MetaName{$(params_sym...),Typ,Names,
                                                            Types}}) where {$(params...),
                                                                            Typ,Names,Types}
            return NamedTuple{($field, Names...),Base.tuple_type_cons(Typ, Types)}
        end

        function StructArrays.createinstance(::Type{$MetaName{$(params_sym...),Typ,Names,
                                                              Types}}, metafree,
                                             args...) where {$(params...),Typ,Names,Types}
            return $MetaName(metafree, NamedTuple{Names,Types}(args))
        end
    end
    return esc(expr)
end

@meta_type(Point, position, AbstractPoint, Dim, T)
Base.getindex(x::PointMeta, idx::Int) = getindex(metafree(x), idx)

@meta_type(NgonFace, ngon, AbstractNgonFace, N, T)
Base.getindex(x::NgonFaceMeta, idx::Int) = getindex(metafree(x), idx)

@meta_type(SimplexFace, simplex, AbstractSimplexFace, N, T)
Base.getindex(x::SimplexFaceMeta, idx::Int) = getindex(metafree(x), idx)

@meta_type(Polygon, polygon, AbstractPolygon, N, T)

@meta_type(LineString, lines, AbstractVector, P <: Line)
Base.getindex(x::LineStringMeta, idx::Int) = getindex(metafree(x), idx)
Base.size(x::LineStringMeta) = size(metafree(x))

@meta_type(MultiPoint, points, AbstractVector, P <: AbstractPoint)
Base.getindex(x::MultiPointMeta, idx::Int) = getindex(metafree(x), idx)
Base.size(x::MultiPointMeta) = size(metafree(x))

@meta_type(MultiLineString, linestrings, AbstractVector, P <: LineString)
Base.getindex(x::MultiLineStringMeta, idx::Int) = getindex(metafree(x), idx)
Base.size(x::MultiLineStringMeta) = size(metafree(x))

@meta_type(MultiPolygon, polygons, AbstractVector, P <: Polygon)
Base.getindex(x::MultiPolygonMeta, idx::Int) = getindex(metafree(x), idx)
Base.size(x::MultiPolygonMeta) = size(metafree(x))

@meta_type(Mesh, mesh, AbstractMesh, Element <: Polytope)
Base.getindex(x::MeshMeta, idx::Int) = getindex(metafree(x), idx)
Base.size(x::MeshMeta) = size(metafree(x))

"""

    MetaT(geometry, meta::NamedTuple)
    MetaT(geometry; meta...)
    
Returns a `MetaT` that holds a geometry and its metadata

`MetaT` acts the same as `Meta` method.
The difference lies in the fact that it is designed to handle
geometries and metadata of different/heterogeneous types.

eg: While a Point MetaGeometry is a `PointMeta`, the MetaT representation is `MetaT{Point}`
The downside being it's not subtyped to `AbstractPoint` like a `PointMeta` is.

Example:
```julia
julia> MetaT(Point(1, 2), city = "Mumbai")
MetaT{Point{2,Int64},(:city,),Tuple{String}}([1, 2], (city = "Mumbai",))
```
"""
struct MetaT{T,Names,Types}
    main::T
    meta::NamedTuple{Names,Types}
end

MetaT(x; kwargs...) = MetaT(x, values(kwargs))

"""

    metafree(x::MetaT)
    metafree(x::Array{MetaT})

Free the MetaT from metadata
i.e. returns the geometry/array of geometries
"""
function metafree(x::MetaT)
    return getfield(x, :main)
end
metafree(x::AbstractVector{<:MetaT}) = map(metafree, x)

"""

    meta(x::MetaT)
    meta(x::Array{MetaT})
    
Returns the metadata of a `MetaT`
"""
function meta(x::MetaT)
    return getfield(x, :meta)
end
meta(x::AbstractVector{<:MetaT}) = map(meta, x)

# helper methods
function Base.getproperty(x::MetaT, field::Symbol)
    return if field == :main
        metafree(x)
    elseif field == :meta
        meta(x)
    else
        getproperty(meta(x), field)
    end
end

Base.propertynames(x::MetaT) = (:main, propertynames(meta(x))...)
getnamestypes(::Type{MetaT{T,Names,Types}}) where {T,Names,Types} = (T, Names, Types)

# explicitly give the "schema" of the object to StructArrays
function StructArrays.staticschema(::Type{F}) where {F<:MetaT}
    T, names, types = getnamestypes(F)
    return NamedTuple{(:main, names...),Base.tuple_type_cons(T, types)}
end

# generate an instance of MetaT type
function StructArrays.createinstance(::Type{F}, x, args...) where {F<:MetaT}
    T, names, types = getnamestypes(F)
    return MetaT(x, NamedTuple{names,types}(args))
end

"""
Puts an iterable of MetaT's into a StructArray 
"""
function meta_table(iter)
    cols = Tables.columntable(iter)
    return meta_table(first(cols), Base.tail(cols))
end

function meta_table(main, meta::NamedTuple{names,types}) where {names,types}
    F = MetaT{eltype(main),names,StructArrays.eltypes(types)}
    return StructArray{F}(; main=main, meta...)
end

Base.getindex(x::MetaT, idx::Int) = getindex(metafree(x), idx)
Base.size(x::MetaT) = size(metafree(x))
