@testset "Discretization" begin
  @testset "FanTriangulation" begin
    pts = P2[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.75, 1.5), (0.25, 1.5), (0.0, 1.0)]
    tris = [Triangle([pts[1], pts[i], pts[i + 1]]) for i in 2:(length(pts) - 1)]
    hex = Hexagon(pts)
    mesh = discretize(hex, FanTriangulation())
    @test nvertices(mesh) == 6
    @test nelements(mesh) == 4
    @test eltype(mesh) <: Triangle
    @test vertices(mesh) == pts
    @test collect(elements(mesh)) == tris
  end

  @testset "RegularDiscretization" begin
    # fix import conflict with Plots
    BezierCurve = Meshes.BezierCurve

    bezier = BezierCurve([P2(0, 0), P2(1, 0), P2(1, 1)])
    mesh = discretize(bezier, RegularDiscretization(10))
    @test nvertices(mesh) == 11
    @test nelements(mesh) == 10
    @test eltype(mesh) <: Segment
    @test nvertices.(mesh) ⊆ [2]

    sphere = Sphere(P2(0, 0), T(1))
    mesh = discretize(sphere, RegularDiscretization(10))
    @test nvertices(mesh) == 10
    @test nelements(mesh) == 10
    @test eltype(mesh) <: Segment
    @test nvertices.(mesh) ⊆ [2]

    sphere = Sphere(P3(0, 0, 0), T(1))
    mesh = discretize(sphere, RegularDiscretization(10))
    @test nvertices(mesh) == 10 * 10 + 2
    @test nelements(mesh) == (10) * (10 - 1) + 2 * (10)
    @test eltype(mesh) <: Ngon
    @test nvertices.(mesh) ⊆ [3, 4]

    ball = Ball(P2(0, 0), T(1))
    mesh = discretize(ball, RegularDiscretization(10))
    @test nvertices(mesh) == 10 * 10 + 1
    @test nelements(mesh) == (10) * (10 - 1) + 10
    @test eltype(mesh) <: Ngon
    @test nvertices.(mesh) ⊆ [3, 4]

    cylsurf = CylinderSurface(
      Plane(P3(0, 0, 0), V3(0, 0, 1)),
      Plane(P3(1, 1, 1), V3(0, 0, 1)),
      T(1)
    )
    mesh = discretize(cylsurf, RegularDiscretization(10))
    @test nvertices(mesh) == 10 * 10 + 2
    @test nelements(mesh) == 10 * (10 - 1) + 2 * 10
    @test eltype(mesh) <: Ngon
    @test nvertices.(mesh) ⊆ [3, 4]
  end

  @testset "Dehn1899" begin
    octa = Octagon(
      P2[
        (0.0, 0.0),
        (0.5, -0.5),
        (1.0, 0.0),
        (1.5, 0.5),
        (1.0, 1.0),
        (0.5, 1.5),
        (0.0, 1.0),
        (-0.5, 0.5)
      ]
    )
    mesh = discretize(octa, Dehn1899())
    @test nvertices(mesh) == 8
    @test nelements(mesh) == 6
    @test eltype(mesh) <: Triangle

    octa = Octagon(
      P3[
        (0.0, 0.0, 0.0),
        (0.5, -0.5, 0.0),
        (1.0, 0.0, 0.0),
        (1.5, 0.5, 0.0),
        (1.0, 1.0, 0.0),
        (0.5, 1.5, 0.0),
        (0.0, 1.0, 0.0),
        (-0.5, 0.5, 0.0)
      ]
    )
    mesh = discretize(octa, Dehn1899())
    @test nvertices(mesh) == 8
    @test nelements(mesh) == 6
    @test eltype(mesh) <: Triangle
  end

  @testset "FIST" begin
    𝒫 = Chain(P2[(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2), (0, 0)])
    @test Meshes.ears(𝒫) == [2, 4, 5]

    𝒫 = Chain(P2[(0, 0), (1, 0), (1, 1), (2, 1), (1, 2), (0, 0)])
    @test Meshes.ears(𝒫) == [2, 4]

    𝒫 = Chain(P2[(0, 0), (1, 0), (1, 1), (1, 2), (0, 0)])
    @test Meshes.ears(𝒫) == [2, 4]

    𝒫 = Chain(P2[(0, 0), (1, 1), (1, 2), (0, 0)])
    @test Meshes.ears(𝒫) == []

    𝒫 = Chain(
      P2[
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
        (0.385557753911967, 0.322338556632868),
        (0.443339268495331, 0.283757618605357)
      ]
    )
    @test Meshes.ears(𝒫) == [1, 3, 5, 6, 8, 10, 12, 14]

    points = P2[(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2), (0, 0)]
    connec = connect.([(4, 5, 6), (3, 4, 6), (3, 6, 1), (1, 2, 3)], Triangle)
    target = SimpleMesh(points[1:(end - 1)], connec)
    poly = PolyArea(points)
    mesh = discretize(poly, FIST(shuffle=false))
    @test mesh == target
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2
  end

  @testset "Miscellaneous" begin
    for method in [FIST(), Dehn1899()]
      triangle = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
      mesh = discretize(triangle, method)
      @test vertices(mesh) == [P2(0, 0), P2(1, 0), P2(0, 1)]
      @test collect(elements(mesh)) == [triangle]

      quadrangle = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
      mesh = discretize(quadrangle, method)
      elms = collect(elements(mesh))
      @test vertices(mesh) == [P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1)]
      @test eltype(elms) <: Triangle
      @test length(elms) == 2

      q = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
      t = Triangle(P2(1, 0), P2(2, 1), P2(1, 1))
      m = Multi([q, t])
      mesh = discretize(m, method)
      elms = collect(elements(mesh))
      @test vertices(mesh) == [vertices(q); vertices(t)]
      @test vertices(elms[1]) ⊆ vertices(q)
      @test vertices(elms[2]) ⊆ vertices(q)
      @test vertices(elms[3]) ⊆ vertices(t)
      @test eltype(elms) <: Triangle
      @test length(elms) == 3

      outer = P2[(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
      hole1 = P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4), (0.2, 0.2)]
      hole2 = P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4), (0.6, 0.2)]
      poly = PolyArea(outer, [hole1, hole2])
      chain, _ = bridge(poly, width=T(0.01))
      mesh = discretizewithin(chain, method)
      @test nvertices(mesh) == 16
      @test nelements(mesh) == 14
      @test all(t -> area(t) > zero(T), mesh)

      # 3D chains
      chain = Chain(P3[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 1), (0, 0, 0)])
      mesh = discretizewithin(chain, method)
      @test vertices(mesh) == vertices(chain)
      @test eltype(mesh) <: Triangle
      @test nelements(mesh) == 2
    end
  end

  @testset "Difficult examples" begin
    for method in [FIST(), Dehn1899()]
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
  end

  @testset "Discretize" begin
    ball = Ball(P2(0, 0), T(1))
    mesh = discretize(ball)
    @test !(eltype(mesh) <: Triangle)
    @test !(eltype(mesh) <: Quadrangle)
    @test nelements(mesh) == 2500

    sphere = Sphere(P3(0, 0, 0), T(1))
    mesh = discretize(ball)
    @test !(eltype(mesh) <: Triangle)
    @test !(eltype(mesh) <: Quadrangle)
    @test nelements(mesh) == 2500

    cylsurf = CylinderSurface(T(1))
    mesh = discretize(cylsurf)
    @test !(eltype(mesh) <: Triangle)
    @test !(eltype(mesh) <: Quadrangle)
    @test nelements(mesh) == 150

    grid = CartesianGrid(10)
    @test discretize(grid) == grid

    mesh = SimpleMesh(rand(P2, 3), connect.([(1, 2, 3)]))
    @test discretize(mesh) == mesh
  end

  @testset "Simplexify" begin
    # fix import conflict with Plots
    BezierCurve = Meshes.BezierCurve

    # simplexify is a helper function that calls an
    # appropriate discretization method depending on
    # the geometry type that is given to it
    box = Box(P1(0), P1(1))
    msh = simplexify(box)
    @test eltype(msh) <: Segment
    @test topology(msh) == GridTopology(1)
    @test nvertices(msh) == 2
    @test nelements(msh) == 1
    @test msh[1] == Segment(P1(0), P1(1))

    seg = Segment(P1(0), P1(1))
    msh = simplexify(seg)
    @test eltype(msh) <: Segment
    @test topology(msh) == GridTopology(1)
    @test nvertices(msh) == 2
    @test nelements(msh) == 1
    @test msh[1] == Segment(P1(0), P1(1))

    chn = Chain(P2[(0, 0), (1, 0), (1, 1)])
    msh = simplexify(chn)
    @test eltype(msh) <: Segment
    @test nvertices(msh) == 3
    @test nelements(msh) == 2
    @test msh[1] == Segment(P2(0, 0), P2(1, 0))
    @test msh[2] == Segment(P2(1, 0), P2(1, 1))
    chn = Chain(P2[(0, 0), (1, 0), (1, 1), (0, 0)])
    msh = simplexify(chn)
    @test eltype(msh) <: Segment
    @test nvertices(msh) == 3
    @test nelements(msh) == 3
    @test msh[1] == Segment(P2(0, 0), P2(1, 0))
    @test msh[2] == Segment(P2(1, 0), P2(1, 1))
    @test msh[3] == Segment(P2(1, 1), P2(0, 0))

    sph = Sphere(P2(0, 0), T(1))
    msh = simplexify(sph)
    @test eltype(msh) <: Segment
    @test nvertices(msh) == nelements(msh)

    bez = BezierCurve(P2[(0, 0), (1, 0), (1, 1)])
    msh = simplexify(bez)
    @test eltype(msh) <: Segment
    @test nvertices(msh) == nelements(msh) + 1

    box = Box(P2(0, 0), P2(1, 1))
    ngon = Quadrangle(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    poly = readpoly(T, joinpath(datadir, "taubin.line"))
    for geom in [box, ngon, poly]
      bound = boundary(geom)
      mesh = simplexify(geom)
      @test Set(vertices(bound)) == Set(vertices(mesh))
      @test nelements(mesh) == length(vertices(mesh)) - 2
    end

    # triangulation of multi geometries
    box1 = Box(P2(0, 0), P2(1, 1))
    box2 = Box(P2(1, 1), P2(2, 2))
    multi = Multi([box1, box2])
    mesh = simplexify(multi)
    @test nvertices(mesh) == 8
    @test nelements(mesh) == 4

    # triangulation of spheres
    sphere = Sphere(P3(0, 0, 0), T(1))
    mesh = simplexify(sphere)
    @test eltype(mesh) <: Triangle
    xs = coordinates.(vertices(mesh))
    @test all(x -> norm(x) ≈ T(1), xs)

    # triangulation of cylinder surfaces
    cylsurf = CylinderSurface(T(1))
    mesh = simplexify(cylsurf)
    @test eltype(mesh) <: Triangle
    xs = coordinates.(vertices(mesh))
    @test all(x -> T(-1) ≤ x[1] ≤ T(1), xs)
    @test all(x -> T(-1) ≤ x[2] ≤ T(1), xs)
    @test all(x -> T(0) ≤ x[3] ≤ T(1), xs)

    # triangulation of balls
    ball = Ball(P2(0, 0), T(1))
    mesh = simplexify(ball)
    @test eltype(mesh) <: Triangle
    xs = coordinates.(vertices(mesh))
    @test all(x -> norm(x) ≤ T(1) + eps(T), xs)

    # triangulation of meshes
    grid = CartesianGrid{T}(3, 3)
    mesh = simplexify(grid)
    gpts = vertices(grid)
    mpts = vertices(mesh)
    @test nvertices(mesh) == 16
    @test nelements(mesh) == 18
    @test collect(mpts) == collect(gpts)
    @test eltype(mesh) <: Triangle
    @test measure(mesh) == measure(grid)

    if visualtests
      fig = Mke.Figure(resolution=(600, 300))
      viz(fig[1, 1], grid, showfacets=true)
      viz(fig[1, 2], mesh, showfacets=true)
      @test_reference "data/triangulate-$T.png" fig
    end
  end
end
