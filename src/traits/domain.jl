# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Domain{Dim,T}

A domain is an indexable collection of geometries (i.e. mesh) or points
where each element can be mapped to a feature vector. For example, a
collection of polygonal areas representing the states of a country can
be seen as a domain. In this case, the features can be the name of the
state, the total population of the state, or any other quantity stored
in each polygon.
"""
abstract type Domain{Dim,T} end

"""
    embeddim(domain)

Return the number of dimensions of the space where the `domain` is embedded.
"""
embeddim(::Type{<:Domain{Dim,T}}) where {Dim,T} = Dim
embeddim(domain::Domain) = embeddim(typeof(domain))

"""
    coordtype(domain)

Return the machine type of each coordinate used to describe the `domain`.
"""
coordtype(::Type{<:Domain{Dim,T}}) where {Dim,T} = T
coordtype(domain::Domain) = coordtype(typeof(domain))

"""
    getindex(domain, ind)

Return the `ind`-th element in the `domain`.
"""
Base.getindex(domain::Domain, ind::Int)

"""
    nelements(domain)

Return the number of elements in the `domain`.
"""
nelements(domain::Domain)

"""
    coordinates!(buff, domain, ind)

Compute the coordinates `buff` of the `ind`-th element in the `domain` in place.
"""
coordinates!(buff, domain::Domain, ind::Int)

# ----------
# FALLBACKS
# ----------

function coordinates!(buff, domain::Domain, inds::AbstractVector{Int})
  for j in eachindex(inds)
    coordinates!(view(buff,:,j), domain, inds[j])
  end
  buff
end

function coordinates(domain::Domain, ind::Int)
  buff = MVector{embeddim(domain),coordtype(domain)}(undef)
  coordinates!(buff, domain, ind)
end

function coordinates(domain::Domain, inds::AbstractVector{Int})
  buff = Matrix{coordtype(domain)}(undef, embeddim(domain), length(inds))
  coordinates!(buff, domain, inds)
end

Base.eltype(domain::Domain) = eltype([domain[i] for i in 1:nelements(domain)])

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, domain::Domain{Dim,T}) where {Dim,T}
  nelm = nelements(domain)
  name = nameof(typeof(domain))
  print(io, "$nelm $name{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", domain::Domain{Dim,T}) where {Dim,T}
  println(io, domain)
  N = nelements(domain)
  I, J = N > 10 ? (5, N-4) : (N, N+1)
  lines = [["  └─$(domain[i])" for i in 1:I]
           (N > 10 ? ["  ⋮"] : [])
           ["  └─$(domain[i])" for i in J:N]]
  print(io, join(lines, "\n"))
end
