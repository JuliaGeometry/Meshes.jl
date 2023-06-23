# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Multi(geoms)

A collection of geometries `geoms` seen as a single [`Geometry`](@ref).

In geographic information systems (GIS) it is common to represent
multiple polygons as a single entity (e.g. country with islands).
"""
struct Multi{Dim,T,G<:Geometry{Dim,T}} <: Geometry{Dim,T}
  geoms::Vector{G}
end

# constructor with iterator of geometries
Multi(geoms) = Multi(collect(geoms))

paramdim(m::Multi) = maximum(paramdim, m.geoms)

vertex(m::Multi, ind) = vertices(m)[ind]

vertices(m::Multi) = [vertex for geom in m.geoms for vertex in vertices(geom)]

nvertices(m::Multi) = sum(nvertices, m.geoms)

Base.unique(m::Multi) = unique!(deepcopy(m))

function Base.unique!(m::Multi)
  foreach(unique!, m.geoms)
  m
end

function centroid(m::Multi)
  cs = coordinates.(centroid.(m.geoms))
  Point(sum(cs) / length(cs))
end

measure(m::Multi) = sum(measure, m.geoms)

Base.length(m::Multi{Dim,T,<:Chain}) where {Dim,T} = measure(m)
area(m::Multi{Dim,T,<:Polygon}) where {Dim,T} = measure(m)
volume(m::Multi{Dim,T,<:Polyhedron}) where {Dim,T} = measure(m)

function boundary(m::Multi)
  bounds = [boundary(geom) for geom in m.geoms]
  valid = filter(!isnothing, bounds)
  isempty(valid) ? nothing : reduce(merge, valid)
end

rings(m::Multi{Dim,T,<:Polygon}) where {Dim,T} = [ring for poly in m.geoms for ring in rings(poly)]

Base.collect(m::Multi) = m.geoms

Base.in(point::Point, m::Multi) = any(geom -> point ∈ geom, m.geoms)

==(m₁::Multi, m₂::Multi) =
  length(m₁.geoms) == length(m₂.geoms) && all(g -> g[1] == g[2], zip(m₁.geoms, m₂.geoms))

function Base.show(io::IO, m::Multi{Dim,T}) where {Dim,T}
  name = prettyname(eltype(m.geoms))
  print(io, "Multi$name{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", m::Multi)
  println(io, m)
  print(io, io_lines(m, "  "))
end
