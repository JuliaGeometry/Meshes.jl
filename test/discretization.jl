@testitem "FanTriangulation" begin
  pts = cart.([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.75, 1.5), (0.25, 1.5), (0.0, 1.0)])
  tris = [Triangle(pts[1], pts[i], pts[i + 1]) for i in 2:(length(pts) - 1)]
  hex = Hexagon(pts...)
  mesh = discretize(hex, FanTriangulation())
  @test nvertices(mesh) == 6
  @test nelements(mesh) == 4
  @test eltype(mesh) <: Triangle
  @test vertices(mesh) == pts
  @test collect(elements(mesh)) == tris
end

@testitem "DehnTriangulation" begin
  octa = Octagon(
    cart(0.2, 0.2),
    cart(0.5, -0.5),
    cart(0.8, 0.2),
    cart(1.5, 0.5),
    cart(0.8, 0.8),
    cart(0.5, 1.5),
    cart(0.2, 0.8),
    cart(-0.5, 0.5)
  )
  mesh = discretize(octa, DehnTriangulation())
  @test nvertices(mesh) == 8
  @test nelements(mesh) == 6
  @test eltype(mesh) <: Triangle

  # type stability tests
  poly = PolyArea(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  @inferred discretize(poly, DehnTriangulation())

  octa = Octagon(
    cart(0.2, 0.2, 0.0),
    cart(0.5, -0.5, 0.0),
    cart(0.8, 0.2, 0.0),
    cart(1.5, 0.5, 0.0),
    cart(0.8, 0.8, 0.0),
    cart(0.5, 1.5, 0.0),
    cart(0.2, 0.8, 0.0),
    cart(-0.5, 0.5, 0.0)
  )
  mesh = discretize(octa, DehnTriangulation())
  @test nvertices(mesh) == 8
  @test nelements(mesh) == 6
  @test eltype(mesh) <: Triangle
end

@testitem "HeldTriangulation" begin
  𝒫 = Ring(cart.([(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2)]))
  @test Meshes.earsccw(𝒫) == [2, 4, 5]

  𝒫 = Ring(cart.([(0, 0), (1, 0), (1, 1), (2, 1), (1, 2)]))
  @test Meshes.earsccw(𝒫) == [2, 4]

  𝒫 = Ring(cart.([(0, 0), (1, 0), (1, 1), (1, 2)]))
  @test Meshes.earsccw(𝒫) == [2, 4]

  𝒫 = Ring(cart.([(0, 0), (1, 1), (1, 2)]))
  @test Meshes.earsccw(𝒫) == []

  𝒫 = Ring(
    cart.([
      (0.443339268495331, 0.283757618605357),
      (0.497822414616971, 0.398142813114205),
      (0.770343126156527, 0.201815462842808),
      (0.761236456732531, 0.330085709922366),
      (0.985658085510286, 0.221530395507904),
      (0.877899962498139, 0.325516131702896),
      (0.561404274882782, 0.540334008885703),
      (0.949459768187313, 0.396227653478068),
      (0.594962560615951, 0.584927547374551),
      (0.324208409133154, 0.607290684450708),
      (0.424085089823892, 0.493532112641353),
      (0.209843417261654, 0.590030658255966),
      (0.27993878548962, 0.525162463476181),
      (0.385557753911967, 0.322338556632868)
    ])
  )
  @test Meshes.earsccw(𝒫) == [1, 3, 5, 6, 8, 10, 12, 14]

  points = cart.([(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2)])
  connec = connect.([(4, 5, 6), (3, 4, 6), (3, 6, 1), (1, 2, 3)], Triangle)
  target = SimpleMesh(points, connec)
  poly = PolyArea(points)
  mesh = discretize(poly, HeldTriangulation(shuffle=false))
  @test mesh == target
  @test Set(vertices(poly)) == Set(vertices(mesh))
  @test nelements(mesh) == length(vertices(mesh)) - 2

  # https://github.com/JuliaGeometry/Meshes.jl/issues/675
  poly = PolyArea(
    cart.([
      (1.1794224993e7, 1.7289506814e7),
      (1.1794045018e7, 1.7289446822e7),
      (1.1793985026e7, 1.7289486817e7),
      (1.1793965029e7, 1.7289586803e7),
      (1.1794105009e7, 1.7289766778e7),
      (1.1794184998e7, 1.7289866764e7),
      (1.179424499e7, 1.728996675e7),
      (1.179424499e7, 1.7290106731e7),
      (1.1794344976e7, 1.7290246711e7),
      (1.1794364973e7, 1.7290386692e7),
      (1.1794504954e7, 1.7290406689e7),
      (1.1794724923e7, 1.729018672e7),
      (1.1794624937e7, 1.7289946753e7),
      (1.1794624937e7, 1.7289806772e7),
      (1.1794564946e7, 1.7289706786e7),
      (1.1794424965e7, 1.7289626797e7)
    ])
  )
  rng = StableRNG(123)
  mesh = discretize(poly, HeldTriangulation(rng))
  @test nvertices(mesh) == 16
  @test nelements(mesh) == 14

  # https://github.com/JuliaGeometry/Meshes.jl/issues/738
  poly = PolyArea(
    cart.([
      (-0.5, 0.3296139),
      (-0.19128194, -0.5),
      (-0.37872985, 0.29592824),
      (0.21377224, -0.0076110554),
      (-0.20127837, 0.24671146)
    ])
  )
  rng = StableRNG(123)
  mesh = discretize(poly, HeldTriangulation(rng))
  @test nvertices(mesh) == 5
  @test nelements(mesh) == 3
end

@testitem "DelaunayTriangulation" begin
  rng = StableRNG(123)
  poly = Pentagon(cart(0, 0), cart(1, 0), cart(1, 1), cart(0.5, 2), cart(0, 1))
  mesh = discretize(poly, DelaunayTriangulation(rng))
  @test Set(vertices(poly)) == Set(vertices(mesh))
  @test nelements(mesh) == length(vertices(mesh)) - 2
end

@testitem "Miscellaneous triangulations" begin
  rng = StableRNG(123)
  for method in [DehnTriangulation(), HeldTriangulation(rng), DelaunayTriangulation(rng)]
    triangle = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
    mesh = discretize(triangle, method)
    @test vertices(mesh) == [cart(0, 0), cart(1, 0), cart(0, 1)]
    @test nelements(mesh) == 1
    @test mesh[1] ≗ triangle

    quadrangle = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    mesh = discretize(quadrangle, method)
    elms = collect(elements(mesh))
    @test vertices(mesh) == [cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1)]
    @test eltype(elms) <: Triangle
    @test length(elms) == 2

    q = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    t = Triangle(cart(1, 0), cart(2, 1), cart(1, 1))
    m = Multi([q, t])
    mesh = discretize(m, method)
    elms = collect(elements(mesh))
    @test vertices(mesh) == [pointify(q); pointify(t)]
    @test vertices(elms[1]) ⊆ vertices(q)
    @test vertices(elms[2]) ⊆ vertices(q)
    @test vertices(elms[3]) ⊆ vertices(t)
    @test eltype(elms) <: Triangle
    @test length(elms) == 3

    outer = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    hole1 = Ring(cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]))
    hole2 = Ring(cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]))
    poly = PolyArea([outer, reverse(hole1), reverse(hole2)])
    bpoly = poly |> Bridge(T(0.01))
    mesh = discretizewithin(boundary(bpoly), method)
    @test nvertices(mesh) == 16
    @test nelements(mesh) == 14
    @test all(t -> area(t) > zero(ℳ)^2, mesh)

    # 3D chains
    chain = Ring(cart.([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 1)]))
    mesh = discretizewithin(chain, method)
    @test vertices(mesh) == vertices(chain)
    @test eltype(mesh) <: Triangle
    @test nelements(mesh) == 2

    # latlon coordinates
    poly = PolyArea(latlon(0, 0), latlon(0, 1), latlon(1, 1), latlon(1, 0))
    mesh = discretize(poly, method)
    @test vertices(mesh) == vertices(poly)
    @test eltype(mesh) <: Triangle
    @test nelements(mesh) == 2

    # preserves order of vertices
    poly = Quadrangle(cart(0, 1, 0), cart(1, 1, 0), cart(1, 0, 0), cart(0, 0, 0))
    mesh = simplexify(poly)
    @test pointify(mesh) == pointify(poly)
  end
end

@testitem "Difficult triangulations" begin
  rng = StableRNG(123)
  for method in [DehnTriangulation(), HeldTriangulation(rng)]
    poly = readpoly(T, joinpath(datadir, "taubin.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2

    poly = readpoly(T, joinpath(datadir, "poly1.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2

    poly = readpoly(T, joinpath(datadir, "poly2.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2

    poly = readpoly(T, joinpath(datadir, "poly3.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2

    poly = readpoly(T, joinpath(datadir, "poly4.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2

    poly = readpoly(T, joinpath(datadir, "poly5.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2

    poly = readpoly(T, joinpath(datadir, "smooth1.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2

    poly = readpoly(T, joinpath(datadir, "smooth2.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2

    poly = readpoly(T, joinpath(datadir, "smooth3.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2

    poly = readpoly(T, joinpath(datadir, "smooth4.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2

    poly = readpoly(T, joinpath(datadir, "smooth5.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2

    # https://github.com/JuliaGeometry/Meshes.jl/issues/738
    poly = readpoly(T, joinpath(datadir, "hole1.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == 32

    poly = readpoly(T, joinpath(datadir, "hole2.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == 30

    poly = readpoly(T, joinpath(datadir, "hole3.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == 32

    poly = readpoly(T, joinpath(datadir, "hole4.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == 30

    poly = readpoly(T, joinpath(datadir, "hole5.line"))
    mesh = discretize(poly, method)
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == 32
  end

  if T == Float64
    poly = PolyArea(
      cart.([
        (-48.03012478813999, -18.323912004531923),
        (-48.030125176275845, -18.323904748608573),
        (-48.03017873307118, -18.323925747019675),
        (-48.03017945243984, -18.32393728592407),
        (-48.030185785831904, -18.32394021501982),
        (-48.03017951837907, -18.323938343610457),
        (-48.030124261780436, -18.32392184444903),
        (-48.0301218833633, -18.323910661117687)
      ])
    )
    mesh = discretize(poly)
    @test nvertices(mesh) == 8
    @test nelements(mesh) == 6
  end

  # degenerate triangle
  poly = PolyArea(cart.([(0, 0), (1, 1), (1, 1)]))
  mesh = discretize(poly)
  @test nvertices(mesh) == 3
  @test nelements(mesh) == 1
  @test vertices(mesh) == [cart(0, 0), cart(0, 0), cart(0, 0)]
  @test mesh[1] == Triangle(cart(0, 0), cart(0, 0), cart(0, 0))
end

@testitem "ManualDiscretization" begin
  box = Box(cart(0, 0, 0), cart(1, 1, 1))
  hexa = Hexahedron(pointify(box)...)
  bmesh = discretize(box, ManualDiscretization())
  hmesh = discretize(hexa, ManualDiscretization())
  @test bmesh == hmesh
  @test nvertices(bmesh) == 8
  @test nelements(bmesh) == 5
end

@testitem "RegularDiscretization" begin
  bezier = BezierCurve([cart(0, 0), cart(1, 0), cart(1, 1)])
  mesh = discretize(bezier, RegularDiscretization(10))
  @test nvertices(mesh) == 11
  @test nelements(mesh) == 10
  @test eltype(mesh) <: Segment
  @test nvertices.(mesh) ⊆ [2]

  box = Box(cart(0, 0), cart(2, 2))
  mesh = discretize(box, RegularDiscretization(10))
  @test mesh isa CartesianGrid
  @test nvertices(mesh) == 121
  @test nelements(mesh) == 100
  @test eltype(mesh) <: Quadrangle
  @test nvertices.(mesh) ⊆ [4]

  sphere = Sphere(cart(0, 0), T(1))
  mesh = discretize(sphere, RegularDiscretization(10))
  @test nvertices(mesh) == 10
  @test nelements(mesh) == 10
  @test eltype(mesh) <: Segment
  @test nvertices.(mesh) ⊆ [2]

  sphere = Sphere(cart(0, 0, 0), T(1))
  mesh = discretize(sphere, RegularDiscretization(10))
  @test nvertices(mesh) == 11 * 10 + 2
  @test nelements(mesh) == 10 * 10 + 2 * 10
  @test eltype(mesh) <: Ngon
  @test nvertices.(mesh) ⊆ [3, 4]

  ellips = Ellipsoid((T(3), T(2), T(1)))
  mesh = discretize(ellips, RegularDiscretization(10))
  @test nvertices(mesh) == 11 * 10 + 2
  @test nelements(mesh) == 10 * 10 + 2 * 10
  @test eltype(mesh) <: Ngon
  @test nvertices.(mesh) ⊆ [3, 4]

  ball = Ball(cart(0, 0), T(1))
  mesh = discretize(ball, RegularDiscretization(10))
  @test nvertices(mesh) == 11 * 10 + 1
  @test nelements(mesh) == 10 * 10 + 10
  @test eltype(mesh) <: Ngon
  @test nvertices.(mesh) ⊆ [3, 4]

  disk = Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(1))
  mesh = discretize(disk, RegularDiscretization(10))
  @test nvertices(mesh) == 11 * 10 + 1
  @test nelements(mesh) == 10 * 10 + 10
  @test eltype(mesh) <: Ngon
  @test nvertices.(mesh) ⊆ [3, 4]

  cyl = Cylinder(Plane(cart(0, 0, 0), vector(0, 0, 1)), Plane(cart(1, 1, 1), vector(0, 0, 1)), T(1))
  mesh = discretize(cyl, RegularDiscretization(10))
  @test nvertices(mesh) == 11 * 10 * 11 + 11
  @test nelements(mesh) == 11 * 10 * 10
  @test eltype(mesh) <: Polyhedron
  @test nvertices.(mesh) ⊆ [6, 8]

  cylsurf = CylinderSurface(Plane(cart(0, 0, 0), vector(0, 0, 1)), Plane(cart(1, 1, 1), vector(0, 0, 1)), T(1))
  mesh = discretize(cylsurf, RegularDiscretization(10))
  @test nvertices(mesh) == 10 * 11 + 2
  @test nelements(mesh) == 10 * 10 + 2 * 10
  @test eltype(mesh) <: Ngon
  @test nvertices.(mesh) ⊆ [3, 4]

  consurf = ConeSurface(Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(1)), cart(0, 0, 1))
  mesh = discretize(consurf, RegularDiscretization(10))
  @test nvertices(mesh) == 10 * 11 + 2
  @test nelements(mesh) == 10 * 10 + 2 * 10
  @test eltype(mesh) <: Ngon
  @test nvertices.(mesh) ⊆ [3, 4]

  parsurf = rand(ParaboloidSurface)
  mesh = discretize(parsurf, RegularDiscretization(10))
  @test nvertices(mesh) == 10 * (10 + 1)
  @test nelements(mesh) == 10 * 10
  @test eltype(mesh) <: Ngon
  @test nvertices.(mesh) ⊆ [3, 4]

  poly = PolyArea(cart.([(0, 0), (0, 1), (1, 2), (2, 1), (2, 0)]))
  mesh = discretize(poly, RegularDiscretization(50))
  @test mesh isa Meshes.SubGrid
  grid = parent(mesh)
  @test grid isa CartesianGrid
  @test eltype(mesh) <: Quadrangle
  @test all(intersects(poly), mesh)
end

@testitem "Discretize" begin
  ball = Ball(cart(0, 0), T(1))
  mesh = discretize(ball)
  @test !(eltype(mesh) <: Triangle)
  @test !(eltype(mesh) <: Quadrangle)
  @test nelements(mesh) == 2550

  sphere = Sphere(cart(0, 0, 0), T(1))
  mesh = discretize(sphere)
  @test !(eltype(mesh) <: Triangle)
  @test !(eltype(mesh) <: Quadrangle)
  @test nelements(mesh) == 2600

  cyl = Cylinder(T(1))
  mesh = discretize(cyl)
  @test !(eltype(mesh) <: Wedge)
  @test !(eltype(mesh) <: Hexahedron)
  @test nelements(mesh) == 300

  cylsurf = CylinderSurface(T(1))
  mesh = discretize(cylsurf)
  @test !(eltype(mesh) <: Triangle)
  @test !(eltype(mesh) <: Quadrangle)
  @test nelements(mesh) == 200

  grid = CartesianGrid(10)
  @test discretize(grid) == grid

  mesh = SimpleMesh(randpoint2(3), connect.([(1, 2, 3)]))
  @test discretize(mesh) == mesh
end

@testitem "Simplexify" begin
  # simplexify is a helper function that calls an
  # appropriate discretization method depending on
  # the geometry type that is given to it
  box = Box(cart(0), cart(1))
  msh = simplexify(box)
  @test eltype(msh) <: Segment
  @test topology(msh) == GridTopology(1)
  @test nvertices(msh) == 2
  @test nelements(msh) == 1
  @test msh[1] == Segment(cart(0), cart(1))

  seg = Segment(cart(0), cart(1))
  msh = simplexify(seg)
  @test eltype(msh) <: Segment
  @test topology(msh) == GridTopology(1)
  @test nvertices(msh) == 2
  @test nelements(msh) == 1
  @test msh[1] == Segment(cart(0), cart(1))

  chn = Rope(cart.([(0, 0), (1, 0), (1, 1)]))
  msh = simplexify(chn)
  @test eltype(msh) <: Segment
  @test nvertices(msh) == 3
  @test nelements(msh) == 2
  @test msh[1] == Segment(cart(0, 0), cart(1, 0))
  @test msh[2] == Segment(cart(1, 0), cart(1, 1))
  chn = Ring(cart.([(0, 0), (1, 0), (1, 1)]))
  msh = simplexify(chn)
  @test eltype(msh) <: Segment
  @test nvertices(msh) == 3
  @test nelements(msh) == 3
  @test msh[1] == Segment(cart(0, 0), cart(1, 0))
  @test msh[2] == Segment(cart(1, 0), cart(1, 1))
  @test msh[3] == Segment(cart(1, 1), cart(0, 0))

  sph = Sphere(cart(0, 0), T(1))
  msh = simplexify(sph)
  @test eltype(msh) <: Segment
  @test nvertices(msh) == nelements(msh)

  bez = BezierCurve(cart.([(0, 0), (1, 0), (1, 1)]))
  msh = simplexify(bez)
  @test eltype(msh) <: Segment
  @test nvertices(msh) == nelements(msh) + 1

  box = Box(cart(0, 0), cart(1, 1))
  ngon = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  poly = readpoly(T, joinpath(datadir, "taubin.line"))
  for geom in [box, ngon, poly]
    bound = boundary(geom)
    mesh = simplexify(geom)
    @test Set(vertices(bound)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2
  end

  # triangulation of multi geometries
  box1 = Box(cart(0, 0), cart(1, 1))
  box2 = Box(cart(1, 1), cart(2, 2))
  multi = Multi([box1, box2])
  mesh = simplexify(multi)
  @test nvertices(mesh) == 8
  @test nelements(mesh) == 4

  # triangulation of spheres
  sphere = Sphere(cart(0, 0, 0), T(1))
  mesh = simplexify(sphere)
  @test eltype(mesh) <: Triangle
  xs = to.(vertices(mesh))
  @test all(x -> norm(x) ≈ oneunit(ℳ), xs)

  # triangulation of cylinder surfaces
  cylsurf = CylinderSurface(T(1))
  mesh = simplexify(cylsurf)
  @test eltype(mesh) <: Triangle
  xs = to.(vertices(mesh))
  @test all(x -> -oneunit(ℳ) ≤ x[1] ≤ oneunit(ℳ), xs)
  @test all(x -> -oneunit(ℳ) ≤ x[2] ≤ oneunit(ℳ), xs)
  @test all(x -> zero(ℳ) ≤ x[3] ≤ oneunit(ℳ), xs)

  # triangulation of balls
  ball = Ball(cart(0, 0), T(1))
  mesh = simplexify(ball)
  @test eltype(mesh) <: Triangle
  xs = to.(vertices(mesh))
  @test all(x -> norm(x) ≤ oneunit(ℳ) + eps(T) * u"m", xs)

  # triangulation of meshes
  grid = cartgrid(3, 3)
  mesh = simplexify(grid)
  gpts = vertices(grid)
  mpts = vertices(mesh)
  @test nvertices(mesh) == 16
  @test nelements(mesh) == 18
  @test collect(mpts) == collect(gpts)
  @test eltype(mesh) <: Triangle
  @test measure(mesh) == measure(grid)

  # https://github.com/JuliaGeometry/Meshes.jl/issues/499
  quad = Quadrangle(cart(0, 1, -1), cart(0, 1, 1), cart(0, -1, 1), cart(0, -1, -1))
  mesh = simplexify(quad)
  @test vertices(mesh) == pointify(quad)

  if visualtests
    grid = cartgrid(3, 3)
    mesh = simplexify(grid)
    fig = Mke.Figure(size=(600, 300))
    viz(fig[1, 1], grid, showsegments=true)
    viz(fig[1, 2], mesh, showsegments=true)
    @test_reference "data/triangulate-$T.png" fig
  end

  # tetrahedralization
  box = Box(cart(0, 0, 0), cart(1, 1, 1))
  hex = Hexahedron(pointify(box)...)
  bmesh = simplexify(box)
  hmesh = simplexify(hex)
  @test bmesh == hmesh
  @test nvertices(bmesh) == 8
  @test nelements(bmesh) == 5
end
