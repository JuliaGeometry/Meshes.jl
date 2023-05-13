@testset "Points" begin
  @test embeddim(Point(1)) == 1
  @test embeddim(Point(1, 2)) == 2
  @test embeddim(Point(1, 2, 3)) == 3
  @test coordtype(Point(1, 1)) == Float64
  @test coordtype(Point(1.0, 1.0)) == Float64
  @test coordtype(Point(1.0f0, 1.0f0)) == Float32
  @test coordtype(Point1(1)) == Float64
  @test coordtype(Point2(1, 1)) == Float64
  @test coordtype(Point3(1, 1, 1)) == Float64
  @test coordtype(Point1f(1)) == Float32
  @test coordtype(Point2f(1, 1)) == Float32
  @test coordtype(Point3f(1, 1, 1)) == Float32

  @test coordtype(Point{2,T}((1, 1))) == T
  @test coordtype(Point{2,T}(1, 1)) == T

  @test coordinates(P1(1)) == T[1]
  @test coordinates(P2(1, 2)) == T[1, 2]
  @test coordinates(P3(1, 2, 3)) == T[1, 2, 3]

  @test P1(1) - P1(1) == T[0]
  @test P2(1, 2) - P2(1, 1) == T[0, 1]
  @test P3(1, 2, 3) - P3(1, 1, 1) == T[0, 1, 2]
  @test_throws DimensionMismatch P2(1, 2) - P3(1, 2, 3)

  @test P1(1) + V1(0) == P1(1)
  @test P1(2) + V1(2) == P1(4)
  @test P2(1, 2) + V2(0, 0) == P2(1, 2)
  @test P2(2, 3) + V2(2, 1) == P2(4, 4)
  @test P3(1, 2, 3) + V3(0, 0, 0) == P3(1, 2, 3)
  @test P3(2, 3, 4) + V3(2, 1, 0) == P3(4, 4, 4)
  @test_throws DimensionMismatch P2(1, 2) + V3(1, 2, 3)

  @test P1(1) - V1(0) == P1(1)
  @test P1(2) - V1(2) == P1(0)
  @test P2(1, 2) - V2(0, 0) == P2(1, 2)
  @test P2(2, 3) - V2(2, 1) == P2(0, 2)
  @test P3(1, 2, 3) - V3(0, 0, 0) == P3(1, 2, 3)
  @test P3(2, 3, 4) - V3(2, 1, 0) == P3(0, 2, 4)

  @test embeddim(rand(P1)) == 1
  @test embeddim(rand(P2)) == 2
  @test embeddim(rand(P3)) == 3
  @test coordtype(rand(P1)) == T
  @test coordtype(rand(P2)) == T
  @test coordtype(rand(P3)) == T

  @test eltype(rand(P1, 3)) == P1
  @test eltype(rand(P2, 3)) == P2
  @test eltype(rand(P3, 3)) == P3

  @test P1(1) ≈ P1(1 + eps(T))
  @test P2(1, 2) ≈ P2(1 + eps(T), T(2))
  @test P3(1, 2, 3) ≈ P3(1 + eps(T), T(2), T(3))

  @test embeddim(Point((1,))) == 1
  @test coordtype(Point((1,))) == Float64
  @test coordtype(Point((1.0,))) == Float64

  @test embeddim(Point((1, 2))) == 2
  @test coordtype(Point((1, 2))) == Float64
  @test coordtype(Point((1.0, 2.0))) == Float64

  @test embeddim(Point((1, 2, 3))) == 3
  @test coordtype(Point((1, 2, 3))) == Float64
  @test coordtype(Point((1.0, 2.0, 3.0))) == Float64

  # check all 1D Point constructors, because those tend to make trouble
  @test Point(1) == Point((1,))
  @test Point{1,T}(-2) == Point{1,T}((-2,))
  @test Point{1,T}(0) == Point{1,T}((0,))

  @test_throws DimensionMismatch Point{2,T}(1)
  @test_throws DimensionMismatch Point{3,T}((2, 3))
  @test_throws DimensionMismatch Point{-3,T}((4, 5, 6))

  # There are 2 cases that throw a MethodError instead of a DimensionMismatch:
  # `Point{1,T}((2,3))` because it tries to take the tuple as a whole and convert to T and:
  # `Point{1,T}(2,3)` which does about the same.
  # I don't think this can reasonably be fixed here without hurting performance

  # check that input of mixed coordinate types is allowed and works as expected
  @test Point(1, 0.2) == Point{2,Float64}(1.0, 0.2)
  @test Point((3.0, 4)) == Point{2,Float64}(3.0, 4.0)
  @test Point((5.0, 6.0, 7)) == Point{3,Float64}(5.0, 6.0, 7.0)
  @test Point{2,T}(8, 9.0) == Point{2,T}((8.0, 9.0))
  @test Point{2,T}((-1.0, -2)) == Point{2,T}((-1, -2))
  @test Point{4,T}((0, -1.0, +2, -4.0)) == Point{4,T}((0.0f0, -1.0f0, +2.0f0, -4.0f0))

  # Integer coordinates converted to Float64
  @test coordtype(Point(1)) == Float64
  @test coordtype(Point(1, 2)) == Float64
  @test coordtype(Point(1, 2, 3)) == Float64

  # generalized inequality
  @test P2(1, 1) ⪯ P2(1, 1)
  @test !(P2(1, 1) ≺ P2(1, 1))
  @test P2(1, 2) ⪯ P2(3, 4)
  @test P2(1, 2) ≺ P2(3, 4)
  @test P2(1, 1) ⪰ P2(1, 1)
  @test !(P2(1, 1) ≻ P2(1, 1))
  @test P2(3, 4) ⪰ P2(1, 2)
  @test P2(3, 4) ≻ P2(1, 2)

  # center and centroid
  @test Meshes.center(P2(1, 1)) == P2(1, 1)
  @test centroid(P2(1, 1)) == P2(1, 1)

  # measure of points is zero
  @test measure(P2(1, 2)) == zero(T)
  @test measure(P3(1, 2, 3)) == zero(T)

  # boundary of points is nothing
  @test isnothing(boundary(rand(P1)))
  @test isnothing(boundary(rand(P2)))
  @test isnothing(boundary(rand(P3)))

  # check broadcasting works as expected
  @test P2(2, 2) .- [P2(2, 3), P2(3, 1)] == [[0.0, -1.0], [-1.0, 1.0]]
  @test P3(2, 2, 2) .- [P3(2, 3, 1), P3(3, 1, 4)] == [[0.0, -1.0, 1.0], [-1.0, 1.0, -2.0]]
end
