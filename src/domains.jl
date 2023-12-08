# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Domain

A domain is an indexable collection of geometries (e.g. mesh).
"""
abstract type Domain{Dim,T} end

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

Base.getindex(domain::Domain, ind::Int) = element(domain, ind)

Base.getindex(domain::Domain, inds::AbstractVector) = [element(domain, ind) for ind in inds]

Base.firstindex(domain::Domain) = 1

Base.lastindex(domain::Domain) = nelements(domain)

Base.length(domain::Domain) = nelements(domain)

Base.iterate(domain::Domain, state=1) = state > nelements(domain) ? nothing : (domain[state], state + 1)

Base.eltype(domain::Domain) = eltype([domain[i] for i in 1:nelements(domain)])

Base.keys(domain::Domain) = 1:nelements(domain)

Base.parent(domain::Domain) = domain

Base.parentindices(domain::Domain) = 1:nelements(domain)

Base.vcat(d1::Domain, d2::Domain) = GeometrySet(vcat(collect(d1), collect(d2)))

Base.vcat(domains::Domain...) = reduce(vcat, domains)

"""
    embeddim(domain)

Return the number of dimensions of the space where the `domain` is embedded.
"""
embeddim(::Type{<:Domain{Dim,T}}) where {Dim,T} = Dim
embeddim(domain::Domain) = embeddim(typeof(domain))

"""
    paramdim(domain)

Return the number of parametric dimensions of the `domain` as the number of
parametric dimensions of its elements.
"""
paramdim(domain::Domain) = paramdim(first(domain))

"""
    coordtype(domain)

Return the machine type of each coordinate used to describe the `domain`.
"""
coordtype(::Type{<:Domain{Dim,T}}) where {Dim,T} = T
coordtype(domain::Domain) = coordtype(typeof(domain))

"""
    centroid(domain, ind)

Return the centroid of the `ind`-th element in the `domain`.
"""
centroid(domain::Domain, ind::Int) = centroid(domain[ind])

"""
    centroid(domain)

Return the centroid of the `domain`, i.e. the centroid of all
its element's centroids.
"""
function centroid(domain::Domain{Dim,T}) where {Dim,T}
  coords(i) = coordinates(centroid(domain, i))
  volume(i) = measure(element(domain, i))
  n = nelements(domain)
  x = coords.(1:n)
  w = volume.(1:n)
  all(iszero, w) && (w = ones(T, n))
  Point(sum(w .* x) / sum(w))
end

"""
    extrema(domain)

Return the top left and bottom right corners of the
bounding box of the `domain`.
"""
Base.extrema(domain::Domain) = extrema(boundingbox(domain))

"""
    topology(domain)

Return the topological structure of the `domain`.
"""
topology(domain::Domain) = domain.topology

# -----------
# IO METHODS
# -----------

function Base.summary(io::IO, domain::Domain{Dim,T}) where {Dim,T}
  nelm = nelements(domain)
  name = prettyname(domain)
  print(io, "$nelm $name{$Dim,$T}")
end

Base.show(io::IO, domain::Domain) = summary(io, domain)

function Base.show(io::IO, ::MIME"text/plain", domain::Domain)
  summary(io, domain)
  println(io)
  printelms(io, domain)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("sets.jl")
include("mesh.jl")
include("trajectories.jl")
