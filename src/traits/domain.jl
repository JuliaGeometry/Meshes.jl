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
    getindex(domain, ind)

Return the `ind`-th element in the `domain`.
"""
Base.getindex(domain::Domain, ind::Int)

"""
    nelements(domain)

Return the number of elements in the `domain`.
"""
function nelements end

# ----------
# FALLBACKS
# ----------

==(d1::Domain, d2::Domain) =
  nelements(d1) == nelements(d2) &&
  all(d1[i] == d2[i] for i in 1:nelements(d1))

Base.firstindex(domain::Domain) = 1

Base.lastindex(domain::Domain) = nelements(domain)

Base.iterate(domain::Domain, state=1) =
  state > nelements(domain) ? nothing : (domain[state], state+1)

Base.eltype(domain::Domain) =
  eltype([domain[i] for i in 1:nelements(domain)])

Base.length(domain::Domain) = nelements(domain)

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
    centroid(domain, ind)

Compute the the centroid of the `ind`-th element in the `domain`.
"""
centroid(domain::Domain, ind::Int) = centroid(domain[ind])

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, domain::Domain{Dim,T}) where {Dim,T}
  nelm = nelements(domain)
  name = nameof(typeof(domain))
  print(io, "$nelm $name{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", domain::Domain)
  println(io, domain)
  N = nelements(domain)
  I, J = N > 10 ? (5, N-4) : (N, N+1)
  lines = [["  └─$(domain[i])" for i in 1:I]
           (N > 10 ? ["  ⋮"] : [])
           ["  └─$(domain[i])" for i in J:N]]
  print(io, join(lines, "\n"))
end
