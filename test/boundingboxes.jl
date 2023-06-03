@testset "Bounding boxes" begin
  @test boundingbox(P2(0, 0)) == Box(P2(0, 0), P2(0, 0))

  @test boundingbox(Box(P2(0, 0), P2(1, 1))) == Box(P2(0, 0), P2(1, 1))

  @test boundingbox(Sphere(P2(0, 0), T(1))) == Box(P2(-1, -1), P2(1, 1))
  @test boundingbox(Sphere(P2(1, 1), T(1))) == Box(P2(0, 0), P2(2, 2))

  @test boundingbox(Ball(P2(0, 0), T(1))) == Box(P2(-1, -1), P2(1, 1))
  @test boundingbox(Ball(P2(1, 1), T(1))) == Box(P2(0, 0), P2(2, 2))

  b = Box(P2(-3, -1), P2(0.5, 0.5))
  s = Sphere(P2(0, 0), T(2))
  @test boundingbox(Multi([b, s])) == Box(P2(-3, -2), P2(2, 2))
  @test boundingbox(Collection([b, s])) == Box(P2(-3, -2), P2(2, 2))

  b1 = Box(P2(0, 0), P2(1, 1))
  b2 = Box(P2(-1, -1), P2(0.5, 0.5))
  @test boundingbox(Multi([b1, b2])) == Box(P2(-1, -1), P2(1, 1))
  @test boundingbox(Collection([b1, b2])) == Box(P2(-1, -1), P2(1, 1))

  @test boundingbox(PointSet(T[0 1 2; 0 2 1])) == Box(P2(0, 0), P2(2, 2))
  @test boundingbox(PointSet(T[1 2; 2 1])) == Box(P2(1, 1), P2(2, 2))

  @test boundingbox(CartesianGrid{T}(10, 10)) == Box(P2(0, 0), P2(10, 10))
  @test boundingbox(CartesianGrid{T}(100, 200)) == Box(P2(0, 0), P2(100, 200))
  @test boundingbox(CartesianGrid((10, 10), T.((1, 1)), T.((1, 1)))) == Box(P2(1, 1), P2(11, 11))

  d = PointSet(T[0 1 2; 0 2 1])
  v = view(d, 1:2)
  @test boundingbox(v) == Box(P2(0, 0), P2(1, 2))

  d = CartesianGrid{T}(10, 10)
  v = view(d, 1:2)
  @test boundingbox(v) == Box(P2(0, 0), P2(2, 1))
end
