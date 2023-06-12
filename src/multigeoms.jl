# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Multi(items)

A collection of points or geometries seen as a single item.

In geographic information systems (GIS) it is common to represent
multiple polygons as a single item. In this case the polygons refer
to the same object in the real world (e.g. country with islands).
"""
struct Multi{Dim,T,I<:PointOrGeometry{Dim,T}} <: Geometry{Dim,T}
  items::Vector{I}
end

# constructor with iterator of items
Multi(items) = Multi(collect(items))

paramdim(multi::Multi) = maximum(paramdim, multi.items)

vertex(multi::Multi, ind) = vertices(multi)[ind]

vertices(multi::Multi) = [vertex for geom in multi.items for vertex in vertices(geom)]

nvertices(multi::Multi) = sum(nvertices, multi.items)

function centroid(multi::Multi)
  cs = coordinates.(centroid.(multi.items))
  Point(sum(cs) / length(cs))
end

measure(multi::Multi) = sum(measure, multi.items)

Base.length(multi::Multi{Dim,T,<:Chain}) where {Dim,T} = measure(multi)
area(multi::Multi{Dim,T,<:Polygon}) where {Dim,T} = measure(multi)
volume(multi::Multi{Dim,T,<:Polyhedron}) where {Dim,T} = measure(multi)

function boundary(multi::Multi)
  bounds = [boundary(geom) for geom in multi.items]
  valid = filter(!isnothing, bounds)
  isempty(valid) ? nothing : reduce(merge, valid)
end

rings(multi::Multi{Dim,T,<:Polygon}) where {Dim,T} = [ring for poly in multi.items for ring in rings(poly)]

Base.collect(multi::Multi) = multi.items

Base.in(point::Point, multi::Multi) = any(geom -> point ∈ geom, multi.items)

==(multi₁::Multi, multi₂::Multi) =
  length(multi₁.items) == length(multi₂.items) && all(g -> g[1] == g[2], zip(multi₁.items, multi₂.items))

function Base.show(io::IO, multi::Multi{Dim,T}) where {Dim,T}
  name = prettyname(eltype(multi.items))
  print(io, "Multi$name{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", multi::Multi)
  println(io, multi)
  print(io, io_lines(multi, "  "))
end
