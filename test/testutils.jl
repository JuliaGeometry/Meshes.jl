# -------------
# HELPER TYPES
# -------------

# meter type
ℳ = Meshes.Met{T}

# dummy type implementing the Domain trait
struct DummyDomain{M<:Meshes.Manifold,C<:CRS} <: Domain{M,C}
  origin::Point{M,C}
end

function Meshes.element(domain::DummyDomain, ind::Int)
  ℒ = Meshes.lentype(domain)
  T = Unitful.numtype(ℒ)
  c = domain.origin + Vec(ntuple(i -> T(ind) * unit(ℒ), embeddim(domain)))
  r = oneunit(ℒ)
  Ball(c, r)
end

Meshes.nelements(d::DummyDomain) = 3

# -------------
# IO FUNCTIONS
# -------------

# helper function to read *.line files containing polygons
# generated with RPG (https://github.com/cgalab/genpoly-rpg)
function readpoly(T, fname)
  open(fname, "r") do f
    # read outer chain
    n = parse(Int, readline(f))
    outer = map(1:n) do _
      coords = readline(f)
      x, y = parse.(T, split(coords))
      Point(x, y)
    end

    # read inner chains
    inners = []
    while !eof(f)
      n = parse(Int, readline(f))
      inner = map(1:n) do _
        coords = readline(f)
        x, y = parse.(T, split(coords))
        Point(x, y)
      end
      push!(inners, inner)
    end

    # return polygonal area
    @assert first(outer) == last(outer)
    @assert all(first(i) == last(i) for i in inners)
    rings = [outer, inners...]
    PolyArea([r[begin:(end - 1)] for r in rings])
  end
end

# helper function to read *.ply files containing meshes
function readply(T, fname)
  ply = load_ply(fname)
  x = T.(ply["vertex"]["x"])
  y = T.(ply["vertex"]["y"])
  z = T.(ply["vertex"]["z"])
  points = Point.(x, y, z)
  connec = [connect(Tuple(c .+ 1)) for c in ply["face"]["vertex_indices"]]
  SimpleMesh(points, connec)
end

# --------------
# CRS FUNCTIONS
# --------------

cart(T::Type, coords...) = cart(T, coords)
cart(T::Type, coords::Tuple) = Point(T.(coords))

merc(T::Type, coords...) = merc(T, coords)
merc(T::Type, coords::Tuple) = Point(Mercator(T.(coords)...))

latlon(T::Type, coords...) = latlon(T, coords)
latlon(T::Type, coords::Tuple) = Point(LatLon(T.(coords)...))

vector(T::Type, coords...) = vector(T, coords)
vector(T::Type, coords::Tuple) = Vec(T.(coords))

cartgrid(args...) = cartgrid(T, args...)
cartgrid(T::Type, dims...) = cartgrid(T, dims)
function cartgrid(T::Type, dims::Dims{Dim}) where {Dim}
  origin = ntuple(i -> T(0.0), Dim)
  spacing = ntuple(i -> T(1.0), Dim)
  CartesianGrid(origin, spacing, GridTopology(dims))
end

randcart(T, Dim, n) = [Point(ntuple(i -> rand(T), Dim)) for _ in 1:n]

# methods with fixed T
cart(xs...) = cart(T, xs...)
merc(xs...) = merc(T, xs...)
latlon(xs...) = latlon(T, xs...)
vector(xs...) = vector(T, xs...)
randpoint1(n) = randcart(T, 1, n)
randpoint2(n) = randcart(T, 2, n)
randpoint3(n) = randcart(T, 3, n)

# ----------------
# OTHER FUNCTIONS
# ----------------

numconvert(T, x::Quantity{S,D,U}) where {S,D,U} = convert(Quantity{T,D,U}, x)

withprecision(_, x) = x
withprecision(T, v::Vec) = numconvert.(T, v)
withprecision(T, p::Point) = Meshes.withcrs(p, withprecision(T, to(p)))
withprecision(T, len::Meshes.Len) = numconvert(T, len)
withprecision(T, lens::NTuple{Dim,Meshes.Len}) where {Dim} = numconvert.(T, lens)
withprecision(T, geoms::StaticVector{Dim,<:Geometry}) where {Dim} = withprecision.(T, geoms)
withprecision(T, geoms::AbstractVector{<:Geometry}) = [withprecision(T, g) for g in geoms]
withprecision(T, geoms::CircularVector{<:Geometry}) = CircularVector([withprecision(T, g) for g in geoms])
@generated function withprecision(T, g::G) where {G<:Meshes.GeometryOrDomain}
  ctor = Meshes.constructor(G)
  names = fieldnames(G)
  exprs = (:(withprecision(T, g.$name)) for name in names)
  :($ctor($(exprs...)))
end

# helper function for type stability tests
function someornone(g1, g2)
  intersection(g1, g2) do I
    if type(I) == NotIntersecting
      "None"
    else
      "Some"
    end
  end
end

setify(lists) = Set(Set.(lists))

function equaltest(g)
  @test g == withprecision(Float64, g)
  @test g == withprecision(Float32, g)
end

isapproxtest(g::Geometry) = _isapproxtest(g, Val(embeddim(g)))

function _isapproxtest(g::Geometry, ::Val{1})
  τ64 = Meshes.atol(Float64) * u"m"
  τ32 = Meshes.atol(Float32) * u"m"
  g64 = withprecision(Float64, g)
  g32 = withprecision(Float32, g)
  @test isapprox(g, Translate(τ64)(g64), atol=1.1τ64)
  @test isapprox(g, Translate(τ32)(g32), atol=1.1τ32)
end

function _isapproxtest(g::Geometry, ::Val{2})
  τ64 = Meshes.atol(Float64) * u"m"
  τ32 = Meshes.atol(Float32) * u"m"
  g64 = withprecision(Float64, g)
  g32 = withprecision(Float32, g)
  @test isapprox(g, Translate(τ64, 0u"m")(g64), atol=1.1τ64)
  @test isapprox(g, Translate(0u"m", τ64)(g64), atol=1.1τ64)
  @test isapprox(g, Translate(τ32, 0u"m")(g32), atol=1.1τ32)
  @test isapprox(g, Translate(0u"m", τ32)(g32), atol=1.1τ32)
end

function _isapproxtest(g::Geometry, ::Val{3})
  τ64 = Meshes.atol(Float64) * u"m"
  τ32 = Meshes.atol(Float32) * u"m"
  g64 = withprecision(Float64, g)
  g32 = withprecision(Float32, g)
  @test isapprox(g, Translate(τ64, 0u"m", 0u"m")(g64), atol=1.1τ64)
  @test isapprox(g, Translate(0u"m", τ64, 0u"m")(g64), atol=1.1τ64)
  @test isapprox(g, Translate(0u"m", 0u"m", τ64)(g64), atol=1.1τ64)
  @test isapprox(g, Translate(τ32, 0u"m", 0u"m")(g32), atol=1.1τ32)
  @test isapprox(g, Translate(0u"m", τ32, 0u"m")(g32), atol=1.1τ32)
  @test isapprox(g, Translate(0u"m", 0u"m", τ32)(g32), atol=1.1τ32)
end

function eachvertexalloc(g)
  iterate(eachvertex(g)) # precompile run
  @allocated for _ in eachvertex(g)
  end
end

function vertextest(g)
  @test collect(eachvertex(g)) == vertices(g)
  @test eachvertexalloc(g) == 0
  # type stability
  @test isconcretetype(eltype(vertices(g)))
  @inferred vertices(g)
end
