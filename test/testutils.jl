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

point(T::Type, coords...) = point(T, coords)
point(T::Type, coords::Tuple) = Point(T.(coords))

vector(T::Type, coords...) = vector(T, coords)
vector(T::Type, coords::Tuple) = Vec(T.(coords))

cartgrid(T::Type, dims...) = cartgrid(T, dims)
function cartgrid(T::Type, dims::Dims{Dim}) where {Dim}
  origin = ntuple(i -> T(0.0), Dim)
  spacing = ntuple(i -> T(1.0), Dim)
  offset = ntuple(i -> 1, Dim)
  CartesianGrid(dims, origin, spacing, offset)
end

randpoint(T, Dim, n) = [Point(ntuple(i -> rand(T), Dim)) for _ in 1:n]

numconvert(T, x::Quantity{S,D,U}) where {S,D,U} = convert(Quantity{T,D,U}, x)

withprecision(_, x) = x
withprecision(T, v::Vec) = numconvert.(T, v)
withprecision(T, p::Point) = Meshes.withdatum(p, withprecision(T, to(p)))
withprecision(T, len::Meshes.Len) = numconvert(T, len)
withprecision(T, lens::NTuple{Dim,Meshes.Len}) where {Dim} = numconvert.(T, lens)
withprecision(T, geoms::NTuple{Dim,<:Geometry}) where {Dim} = withprecision.(T, geoms)
withprecision(T, geoms::AbstractVector{<:Geometry}) = [withprecision(T, g) for g in geoms]
withprecision(T, geoms::CircularVector{<:Geometry}) = CircularVector([withprecision(T, g) for g in geoms])
@generated function withprecision(T, g::G) where {G<:Meshes.GeometryOrDomain}
  ctor = Meshes.constructor(G)
  names = fieldnames(G)
  exprs = (:(withprecision(T, g.$name)) for name in names)
  :($ctor($(exprs...)))
end

function equaltest(g)
  @test g == withprecision(Float64, g)
  @test g == withprecision(Float32, g)
end

function isapproxtest(g::Geometry{1})
  τ64 = Meshes.atol(Float64) * u"m"
  τ32 = Meshes.atol(Float32) * u"m"
  g64 = withprecision(Float64, g)
  g32 = withprecision(Float32, g)
  @test isapprox(g, Translate(τ64)(g64), atol=1.1τ64)
  @test isapprox(g, Translate(τ32)(g32), atol=1.1τ32)
end

function isapproxtest(g::Geometry{2})
  τ64 = Meshes.atol(Float64) * u"m"
  τ32 = Meshes.atol(Float32) * u"m"
  g64 = withprecision(Float64, g)
  g32 = withprecision(Float32, g)
  @test isapprox(g, Translate(τ64, 0u"m")(g64), atol=1.1τ64)
  @test isapprox(g, Translate(0u"m", τ64)(g64), atol=1.1τ64)
  @test isapprox(g, Translate(τ32, 0u"m")(g32), atol=1.1τ32)
  @test isapprox(g, Translate(0u"m", τ32)(g32), atol=1.1τ32)
end

function isapproxtest(g::Geometry{3})
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
