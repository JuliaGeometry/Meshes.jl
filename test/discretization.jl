using Base: datatype_haspadding
@testset "Discretization" begin
  @testset "FIST" begin
    𝒫 = Chain(P2[(0,0),(1,0),(1,1),(2,1),(2,2),(1,2),(0,0)])
    @test Meshes.ears(𝒫) == [2,4,5]

    𝒫 = Chain(P2[(0,0),(1,0),(1,1),(2,1),(1,2),(0,0)])
    @test Meshes.ears(𝒫) == [2,4]

    𝒫 = Chain(P2[(0,0),(1,0),(1,1),(1,2),(0,0)])
    @test Meshes.ears(𝒫) == [2,4]

    𝒫 = Chain(P2[(0,0),(1,1),(1,2),(0,0)])
    @test Meshes.ears(𝒫) == []

    𝒫 = Chain(P2[(0.443339268495331, 0.283757618605357),
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
                 (0.443339268495331, 0.283757618605357)])
    @test Meshes.ears(𝒫) == [1,3,5,6,8,10,12,14]

    points = P2[(0,0),(1,0),(1,1),(2,1),(2,2),(1,2),(0,0)]
    connec = connect.([(4,5,6),(3,4,6),(3,6,1),(1,2,3)], Triangle)
    target = SimpleMesh(points[1:end-1], connec)
    poly = PolyArea(points)
    mesh = discretize(poly, FIST(shuffle=false))
    @test mesh == target
    @test Set(vertices(poly)) == Set(vertices(mesh))
    @test nelements(mesh) == length(vertices(mesh)) - 2
  end

  @testset "Miscellaneous" begin
    for method in [FIST(), Dehn1899()]
      triangle = Triangle(P2(0,0), P2(1,0), P2(0,1))
      mesh = discretize(triangle, method)
      @test vertices(mesh) == [P2(0,0), P2(1,0), P2(0,1)]
      @test collect(elements(mesh)) == [triangle]

      quadrangle = Quadrangle(P2(0,0), P2(1,0), P2(1,1), P2(0,1))
      mesh = discretize(quadrangle, method)
      elms = collect(elements(mesh))
      @test vertices(mesh) == [P2(0,0), P2(1,0), P2(1,1), P2(0,1)]
      @test eltype(elms) <: Triangle
      @test length(elms) == 2

      q = Quadrangle(P2(0,0), P2(1,0), P2(1,1), P2(0,1))
      t = Triangle(P2(1,0), P2(2,1), P2(1,1))
      m = Multi([q, t])
      mesh = discretize(m, method)
      elms = collect(elements(mesh))
      @test vertices(mesh) == [vertices(q); vertices(t)]
      @test vertices(elms[1]) ⊆ vertices(q)
      @test vertices(elms[2]) ⊆ vertices(q)
      @test vertices(elms[3]) ⊆ vertices(t)
      @test eltype(elms) <: Triangle
      @test length(elms) == 3

      outer = P2[(0,0),(1,0),(1,1),(0,1),(0,0)]
      hole1 = P2[(0.2,0.2),(0.4,0.2),(0.4,0.4),(0.2,0.4),(0.2,0.2)]
      hole2 = P2[(0.6,0.2),(0.8,0.2),(0.8,0.4),(0.6,0.4),(0.6,0.2)]
      poly  = PolyArea(outer, [hole1, hole2])
      chain, _ = bridge(poly, width=0.01)
      mesh  = discretize(chain, method)
      @test nvertices(mesh) == 16
      @test nelements(mesh) == 14
      @test all(t -> area(t) > zero(T), mesh)
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

  @testset "Triangulate" begin
    # triangulate is a helper function that calls an
    # appropriate discretization method depending on
    # the geometry type that is given to it
    box  = Box(P2(0,0), P2(1,1))
    ngon = Quadrangle(P2[(0,0),(1,0),(1,1),(0,1)])
    poly = readpoly(T, joinpath(datadir, "taubin.line"))
    for geom in [box, ngon, poly]
      mesh = triangulate(geom)
      @test Set(vertices(geom)) == Set(vertices(mesh))
      @test nelements(mesh) == length(vertices(mesh)) - 2
    end

    # triangulation of multi geometries
    box1  = Box(P2(0,0), P2(1,1))
    box2  = Box(P2(1,1), P2(2,2))
    multi = Multi([box1, box2])
    mesh  = triangulate(multi)
    @test nvertices(mesh) == 8
    @test nelements(mesh) == 4
  end
end
