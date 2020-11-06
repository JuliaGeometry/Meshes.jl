@testset "Points" begin
  @test embeddim(Point(1, 2)) == 2
  @test embeddim(Point(1, 2, 3)) == 3
  @test coordtype(Point(1, 1)) == Int
  @test coordtype(Point(1.,1.)) == Float64
  @test coordtype(Point(1f0, 1f0)) == Float32
  @test coordtype(Point2(1, 1)) == Float64
  @test coordtype(Point3(1, 1, 1)) == Float64
  @test coordtype(Point2f(1, 1)) == Float32
  @test coordtype(Point3f(1, 1, 1)) == Float32

  @test coordtype(Point{2,T}((1, 1))) == T
  @test coordtype(Point{2,T}(1, 1)) == T

  @test coordinates(P2(1, 2)) == [1.0, 2.0]
  @test coordinates(P3(1, 2, 3)) == [1.0, 2.0, 3.0]

  @test P2(1, 2) - P2(1, 1) == T[0, 1]
  @test P3(1, 2, 3) - P3(1, 1, 1) == T[0, 1, 2]
  @test_throws DimensionMismatch P2(1, 2) - P3(1, 2, 3)

  @test P2(1, 2) + V2(0, 0) == P2(1, 2)
  @test P2(2, 3) + V2(2, 1) == P2(4, 4)
  @test P3(1, 2, 3) + V3(0, 0, 0) == P3(1, 2, 3)
  @test P3(2, 3, 4) + V3(2, 1, 0) == P3(4, 4, 4)
  @test_throws DimensionMismatch P2(1, 2) + V3(1, 2, 3)

  @test P2(1, 2) - V2(0, 0) == P2(1, 2)
  @test P2(2, 3) - V2(2, 1) == P2(0, 2)
  @test P3(1, 2, 3) - V3(0, 0, 0) == P3(1, 2, 3)
  @test P3(2, 3, 4) - V3(2, 1, 0) == P3(0, 2, 4)

  @test embeddim(rand(P2)) == 2
  @test embeddim(rand(P3)) == 3
  @test coordtype(rand(P2)) == T
  @test coordtype(rand(P3)) == T

  @test eltype(rand(P2, 3)) == P2
  @test eltype(rand(P3, 3)) == P3

  @test P2(1, 2) ≈ P2(1 + eps(T), T(2))
  @test P3(1, 2, 3) ≈ P3(1 + eps(T), T(2), T(3))

  @test embeddim(Point([1,2])) == 2
  @test coordtype(Point([1,2])) == Int
  @test coordtype(Point([1.,2.])) == Float64

  @test embeddim(Point([1,2,3])) == 3
  @test coordtype(Point([1,2,3])) == Int
  @test coordtype(Point([1.,2.,3.])) == Float64
end
