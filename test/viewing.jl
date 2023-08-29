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
end
