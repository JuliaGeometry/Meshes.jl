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

paramdim(multi::Multi) = maximum(paramdim, multi.geoms)

vertex(multi::Multi, ind) = vertices(multi)[ind]

vertices(multi::Multi) = [vertex for geom in multi.geoms for vertex in vertices(geom)]

nvertices(multi::Multi) = sum(nvertices, multi.geoms)

function centroid(multi::Multi)
  cs = coordinates.(centroid.(multi.geoms))
  Point(sum(cs) / length(cs))
end

measure(multi::Multi) = sum(measure, multi.geoms)

Base.length(multi::Multi{Dim,T,<:Chain}) where {Dim,T} = measure(multi)
area(multi::Multi{Dim,T,<:Polygon}) where {Dim,T} = measure(multi)
volume(multi::Multi{Dim,T,<:Polyhedron}) where {Dim,T} = measure(multi)

function boundary(multi::Multi)
  bounds = [boundary(geom) for geom in multi.geoms]
  valid = filter(!isnothing, bounds)
  isempty(valid) ? nothing : reduce(merge, valid)
end

rings(multi::Multi{Dim,T,<:Polygon}) where {Dim,T} = [ring for poly in multi.geoms for ring in rings(poly)]

Base.collect(multi::Multi) = multi.geoms

Base.in(point::Point, multi::Multi) = any(geom -> point ∈ geom, multi.geoms)

==(multi₁::Multi, multi₂::Multi) =
  length(multi₁.geoms) == length(multi₂.geoms) && all(g -> g[1] == g[2], zip(multi₁.geoms, multi₂.geoms))

function Base.show(io::IO, multi::Multi{Dim,T}) where {Dim,T}
  name = prettyname(eltype(multi.geoms))
  print(io, "Multi$name{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", multi::Multi)
  println(io, multi)
  print(io, io_lines(multi, "  "))
end
