@testset "Bounding boxes" begin
  function bounding_box_tests()
    p = P2(0, 0)
    @test boundingbox(p) == Box(p, p)
    @testnoallocations boundingbox(p)

    b = Box(P2(0, 0), P2(1, 1))
    @test boundingbox(b) == b
    @testnoallocations boundingbox(b)

    s = Sphere(P2(0, 0), T(1))
    @test boundingbox(s) == Box(P2(-1, -1), P2(1, 1))
    @testnoallocations boundingbox(s)
    s = Sphere(P2(1, 1), T(1))
    @test boundingbox(s) == Box(P2(0, 0), P2(2, 2))
    @testnoallocations boundingbox(s)

    b = Ball(P2(0, 0), T(1))
    @test boundingbox(b) == Box(P2(-1, -1), P2(1, 1))
    @testnoallocations boundingbox(b)
    b = Ball(P2(1, 1), T(1))
    @test boundingbox(b) == Box(P2(0, 0), P2(2, 2))
    @testnoallocations boundingbox(b)

    b = Box(P2(-3, -1), P2(0.5, 0.5))
    s = Sphere(P2(0, 0), T(2))
    m = Multi([b, s])
    d = GeometrySet([b, s])
    @test boundingbox(m) == Box(P2(-3, -2), P2(2, 2))
    @test boundingbox(d) == Box(P2(-3, -2), P2(2, 2))
    @test @allocated(boundingbox(m)) < 2500
    @test @allocated(boundingbox(d)) < 2500

    b1 = Box(P2(0, 0), P2(1, 1))
    b2 = Box(P2(-1, -1), P2(0.5, 0.5))
    m = Multi([b1, b2])
    d = GeometrySet([b1, b2])
    @test boundingbox(m) == Box(P2(-1, -1), P2(1, 1))
    @test boundingbox(d) == Box(P2(-1, -1), P2(1, 1))
    @testnoallocations boundingbox(m)
    @testnoallocations boundingbox(d)

    d = PointSet(T[0 1 2; 0 2 1])
    @test boundingbox(d) == Box(P2(0, 0), P2(2, 2))
    @test @allocated(boundingbox(d)) < 50
    d = PointSet(T[1 2; 2 1])
    @test boundingbox(d) == Box(P2(1, 1), P2(2, 2))
    @test @allocated(boundingbox(d)) < 50

    d = CartesianGrid{T}(10, 10)
    @test boundingbox(d) == Box(P2(0, 0), P2(10, 10))
    @testnoallocations boundingbox(d)
    d = CartesianGrid{T}(100, 200)
    @test boundingbox(d) == Box(P2(0, 0), P2(100, 200))
    @testnoallocations boundingbox(d)
    d = CartesianGrid((10, 10), T.((1, 1)), T.((1, 1)))
    @test boundingbox(d) == Box(P2(1, 1), P2(11, 11))
    @testnoallocations boundingbox(d)

    d = PointSet(T[0 1 2; 0 2 1])
    v = view(d, 1:2)
    @test boundingbox(v) == Box(P2(0, 0), P2(1, 2))
    @test @allocated(boundingbox(v)) < 50

    d = CartesianGrid{T}(10, 10)
    v = view(d, 1:2)
    @test boundingbox(v) == Box(P2(0, 0), P2(2, 1))
    @test @allocated(boundingbox(v)) < 9000

    g = CartesianGrid{T}(10, 10)
    m = convert(SimpleMesh, g)
    @test boundingbox(m) == Box(P2(0, 0), P2(10, 10))
    @testnoallocations boundingbox(m)
  end

  bounding_box_tests()
end
