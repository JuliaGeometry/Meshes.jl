# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Discretization{Dim,T}

A discretization is an indexable collection of geometries (i.e. mesh) or points.
"""
abstract type Discretization{Dim,T} end

"""
    embeddim(discretization)

Return the number of dimensions of the space where the `discretization` is embedded.
"""
embeddim(::Type{<:Discretization{Dim,T}}) where {Dim,T} = Dim
embeddim(d::Discretization) = embeddim(typeof(d))

"""
    coordtype(discretization)

Return the machine type of each coordinate used to describe the `discretization`.
"""
coordtype(::Type{<:Discretization{Dim,T}}) where {Dim,T} = T
coordtype(d::Discretization) = coordtype(typeof(d))

"""
    getindex(discretization, ind)

Return the `ind`-th element in the `discretization`.
"""
Base.getindex(d::Discretization, ind::Int)

"""
    nelements(discretization)

Return the number of elements in the `discretization`.
"""
nelements(d::Discretization)

"""
    coordinates!(buff, discretization, ind)

Compute the coordinates `buff` of the `ind`-th element in the `discretization` in place.
"""
coordinates!(buff, d::Discretization, ind::Int)

# ----------
# FALLBACKS
# ----------

function coordinates!(buff, d::Discretization, inds::AbstractVector{Int})
  for j in eachindex(inds)
    coordinates!(view(buff,:,j), d, inds[j])
  end
  buff
end

function coordinates(d::Discretization, ind::Int)
  buff = MVector{embeddim(d),coordtype(d)}(undef)
  coordinates!(buff, d, ind)
end

function coordinates(d::Discretization, inds::AbstractVector{Int})
  buff = Matrix{coordtype(d)}(undef, embeddim(d), length(inds))
  coordinates!(buff, d, inds)
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, d::Discretization{Dim,T}) where {Dim,T}
  nelm = nelements(d)
  name = nameof(typeof(d))
  print(io, "$nelm $name{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", d::Discretization{Dim,T}) where {Dim,T}
  println(io, d)
  N = nelements(d)
  I, J = N > 10 ? (5, N-4) : (N, N+1)
  lines = [["  └─$(d[i])" for i in 1:I]
           (N > 10 ? ["  ⋮"] : [])
           ["  └─$(d[i])" for i in J:N]]
  print(io, join(lines, "\n"))
end
