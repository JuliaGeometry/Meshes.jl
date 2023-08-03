@testset "Viewing" begin
  @testset "Domain" begin
    g = CartesianGrid{T}(10, 10)
    v = view(g, 1:3)
    @test unview(v) == (g, 1:3)
    @test unview(g) == (g, 1:100)

    g = CartesianGrid{T}(10, 10)
    b = Box(P2(1, 1), P2(5, 5))
    v = view(g, b)
    @test v == CartesianGrid(P2(1, 1), P2(5, 5), dims=(4, 4))

    p = PointSet(collect(vertices(g)))
    v = view(p, b)
    @test centroid(v, 1) == P2(1, 1)
    @test centroid(v, nelements(v)) == P2(5, 5)

    g = CartesianGrid{T}(10, 10)
    p = PointSet(collect(vertices(g)))
    b = Ball(P2(0, 0), T(2))
    v = view(g, b)
    @test nelements(v) == 1
    @test v[1] == g[1]
    v = view(p, b)
    @test nelements(v) == 6
    @test coordinates.(v) == V2[(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (0, 2)]

    # convex polygons
    t = Triangle(P2(5, 7), P2(10, 12), P2(15, 7))
    q = Pentagon(P2(6, 1), P2(2, 10), P2(10, 16), P2(18, 10), P2(14, 1))
    p = PolyArea(pointify(q), [pointify(t)])

    g = CartesianGrid{T}(20, 20)
    vt = view(g, t)
    vq = view(g, q)
    vp = view(g, p)
    @test nelements(vt) == 36
    @test nelements(vq) == 162
    @test nelements(vp) == 178

    g = CartesianGrid((10, 10), T.((0, 0)), T.((2, 2)))
    vt = view(g, t)
    vq = view(g, q)
    vp = view(g, p)
    @test nelements(vt) == 11
    @test nelements(vq) == 45
    @test nelements(vp) == 44

    g = CartesianGrid(P2(-2, -2), P2(20, 20), T.((0.5, 1.5)))
    vt = view(g, t)
    vq = view(g, q)
    vp = view(g, p)
    @test nelements(vt) == 61
    @test nelements(vq) == 241
    @test nelements(vp) == 218
  end

  @testset "Data" begin
    dummydata(domain, table) = DummyData(domain, Dict(paramdim(domain) => table))
    dummymeta(domain, table) = meshdata(domain, Dict(paramdim(domain) => table))

    for dummy in [dummydata, dummymeta]
      g = CartesianGrid{T}(10, 10)
      t = (a=1:100, b=1:100)
      d = dummy(g, t)
      v = view(d, 1:3)
      @test unview(v) == (d, 1:3)
      @test unview(d) == (d, 1:100)

      g = CartesianGrid{T}(10, 10)
      t = (a=1:100, b=1:100)
      d = dummy(g, t)
      b = Box(P2(1, 1), P2(5, 5))
      v = view(d, b)
      @test domain(v) == CartesianGrid(P2(1, 1), P2(5, 5), dims=(4, 4))
      @test Tables.columntable(values(v)) == (
        a=[12, 13, 14, 15, 22, 23, 24, 25, 32, 33, 34, 35, 42, 43, 44, 45],
        b=[12, 13, 14, 15, 22, 23, 24, 25, 32, 33, 34, 35, 42, 43, 44, 45]
      )

      p = PointSet(collect(vertices(g)))
      d = dummy(p, t)
      v = view(d, b)
      dd = domain(v)
      @test centroid(dd, 1) == P2(1, 1)
      @test centroid(dd, nelements(dd)) == P2(5, 5)
      tt = Tables.columntable(values(v))
      @test tt == (
        a=[13, 14, 15, 16, 17, 24, 25, 26, 27, 28, 35, 36, 37, 38, 39, 46, 47, 48, 49, 50, 57, 58, 59, 60, 61],
        b=[13, 14, 15, 16, 17, 24, 25, 26, 27, 28, 35, 36, 37, 38, 39, 46, 47, 48, 49, 50, 57, 58, 59, 60, 61]
      )

      g = CartesianGrid{T}(250, 250)
      t = (a=rand(250 * 250), b=rand(250 * 250))
      d = dummy(g, t)
      s1 = slice(d, T(50.5):T(100.2), T(41.7):T(81.3))
      d1 = domain(s1)
      pts1 = [centroid(d1, i) for i in 1:nelements(d1)]
      X1 = reduce(hcat, coordinates.(pts1))
      @test all(T[50.5, 41.7] .≤ minimum(X1, dims=2))
      @test all(maximum(X1, dims=2) .≤ T[100.2, 81.3])

      p = sample(d, 100)
      s2 = slice(p, T(50.5):T(150.7), T(175.2):T(250.3))
      d2 = domain(s2)
      pts2 = [centroid(d2, i) for i in 1:nelements(d2)]
      X2 = reduce(hcat, coordinates.(pts2))
      @test all(T[50.5, 175.2] .≤ minimum(X2, dims=2))
      @test all(maximum(X2, dims=2) .≤ T[150.7, 250.3])
    end
  end
end
