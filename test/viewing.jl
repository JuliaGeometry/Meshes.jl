@testset "Viewing" begin
  @testset "Domain" begin
    g = CartesianGrid{T}(10, 10)
    v = view(g, 1:3)
    @test unview(v) == (g, 1:3)
    @test unview(g) == (g, 1:100)

    g = CartesianGrid{T}(10, 10)
    b = Box(P2(1, 1), P2(5, 5))
    v = view(g, b)
    @test v == CartesianGrid(P2(0, 0), P2(6, 6), dims=(6, 6))

    p = PointSet(collect(vertices(g)))
    v = view(p, b)
    @test centroid(v, 1) == P2(1, 1)
    @test centroid(v, nelements(v)) == P2(5, 5)

    g = CartesianGrid{T}(10, 10)
    p = PointSet(collect(vertices(g)))
    b = Ball(P2(0, 0), T(2))
    v = view(g, b)
    @test nelements(v) == 4
    @test v[1] == g[1]
    v = view(p, b)
    @test nelements(v) == 6
    @test coordinates.(v) == V2[(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (0, 2)]
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
      @test domain(v) == CartesianGrid(P2(0, 0), P2(6, 6), dims=(6, 6))
      x = [collect(1:6); collect(11:16); collect(21:26); collect(31:36); collect(41:46); collect(51:56)]
      @test Tables.columntable(values(v)) == (a=x, b=x)

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
    end
  end
end
