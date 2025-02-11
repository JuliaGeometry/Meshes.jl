# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Domain{M,CRS}

A domain is an indexable collection of geometries (e.g. mesh)
in a given manifold `M` with point coordinates specified in a
coordinate reference system `CRS`.
"""
abstract type Domain{M<:Manifold,C<:CRS} end

"""
    element(domain, ind)

Return the `ind`-th element in the `domain`.
"""
element(domain::Domain, ind::Int)

"""
    nelements(domain)

Return the number of elements in the `domain`.
"""
function nelements end

# ----------
# FALLBACKS
# ----------

==(d1::Domain, d2::Domain) = nelements(d1) == nelements(d2) && all(d1[i] == d2[i] for i in 1:nelements(d1))

Base.isapprox(d1::Domain, d2::Domain) = nelements(d1) == nelements(d2) && all(d1[i] ≈ d2[i] for i in 1:nelements(d1))

Base.getindex(d::Domain, ind::Int) = element(d, ind)

Base.getindex(d::Domain, inds::AbstractVector) = [element(d, ind) for ind in inds]

Base.firstindex(d::Domain) = 1

Base.lastindex(d::Domain) = nelements(d)

Base.length(d::Domain) = nelements(d)

Base.iterate(d::Domain, state=1) = state > nelements(d) ? nothing : (d[state], state + 1)

Base.eltype(d::Domain) = eltype([d[i] for i in 1:nelements(d)])

Base.keys(d::Domain) = 1:nelements(d)

Base.parent(d::Domain) = d

Base.parentindices(d::Domain) = 1:nelements(d)

Base.vcat(d1::Domain, d2::Domain) = GeometrySet(vcat(collect(d1), collect(d2)))

Base.vcat(ds::Domain...) = reduce(vcat, ds)

"""
    embeddim(domain)

Return the number of dimensions of the space where the `domain` is embedded.
"""
embeddim(::Type{<:Domain{M,CRS}}) where {M,CRS} = CoordRefSystems.ndims(CRS)
embeddim(d::Domain) = embeddim(typeof(d))

"""
    paramdim(domain)

Return the number of parametric dimensions of the `domain` as the number of
parametric dimensions of its elements.
"""
paramdim(d::Domain) = paramdim(first(d))

"""
    crs(domain)

Return the coordinate reference system (CRS) of the `domain`.
"""
crs(::Type{<:Domain{M,CRS}}) where {M,CRS} = CRS
crs(d::Domain) = crs(typeof(d))

"""
    manifold(domain)

Return the manifold where the `domain` is defined.
"""
manifold(::Type{<:Domain{M,CRS}}) where {M,CRS} = M
manifold(d::Domain) = manifold(typeof(d))

"""
    lentype(domain)

Return the length type of the `domain`.
"""
lentype(::Type{<:Domain{M,CRS}}) where {M,CRS} = lentype(CRS)
lentype(d::Domain) = lentype(typeof(d))

"""
    centroid(domain, ind)

Return the centroid of the `ind`-th element in the `domain`.
"""
centroid(d::Domain, ind::Int) = centroid(d[ind])

"""
    centroid(domain)

Return the centroid of the `domain`, i.e. the centroid of all
its element's centroids.
"""
function centroid(d::Domain)
  vector(i) = to(centroid(d, i))
  volume(i) = measure(element(d, i))
  n = nelements(d)
  x = vector.(1:n)
  w = volume.(1:n)
  all(iszero, w) && (w = ones(eltype(w), n))
  withcrs(d, sum(w .* x) / sum(w))
end

"""
    extrema(domain)

Return the top left and bottom right corners of the
bounding box of the `domain`.
"""
Base.extrema(d::Domain) = extrema(boundingbox(d))

"""
    topology(domain)

Return the topological structure of the `domain`.
"""
topology(d::Domain) = d.topology

# -----------
# IO METHODS
# -----------

function Base.summary(io::IO, d::Domain)
  nelm = nelements(d)
  name = prettyname(d)
  print(io, "$nelm $name")
end

Base.show(io::IO, d::Domain) = summary(io, d)

function Base.show(io::IO, ::MIME"text/plain", d::Domain)
  summary(io, d)
  println(io)
  printelms(io, d)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("domains/sets.jl")
include("domains/meshes.jl")
include("domains/trajecs.jl")
include("domains/subdomains.jl")

# ------------
# CONVERSIONS
# ------------

Base.convert(::Type{GeometrySet}, d::Domain) = GeometrySet(collect(d))

Base.convert(::Type{SimpleMesh}, m::Mesh) = SimpleMesh(vertices(m), topology(m))

Base.convert(::Type{StructuredGrid}, g::Grid) = StructuredGrid{datum(crs(g))}(XYZ(g))

Base.convert(::Type{RectilinearGrid}, g::CartesianGrid) = RectilinearGrid{datum(crs(g))}(xyz(g))
