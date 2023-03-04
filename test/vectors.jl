@testset "Vectors" begin
  # vararg constructors
  @test eltype(Vec(1, 1)) == Float64
  @test eltype(Vec(1.0, 1.0)) == Float64
  @test eltype(Vec(1.0f0, 1.0f0)) == Float32
  @test eltype(Vec1(1)) == Float64
  @test eltype(Vec2(1, 1)) == Float64
  @test eltype(Vec3(1, 1, 1)) == Float64
  @test eltype(Vec1f(1)) == Float32
  @test eltype(Vec2f(1, 1)) == Float32
  @test eltype(Vec3f(1, 1, 1)) == Float32

  # tuple constructors
  @test eltype(Vec((1, 1))) == Float64
  @test eltype(Vec((1.0, 1.0))) == Float64
  @test eltype(Vec((1.0f0, 1.0f0))) == Float32
  @test eltype(Vec1((1,))) == Float64
  @test eltype(Vec2((1, 1))) == Float64
  @test eltype(Vec3((1, 1, 1))) == Float64
  @test eltype(Vec1f((1,))) == Float32
  @test eltype(Vec2f((1, 1))) == Float32
  @test eltype(Vec3f((1, 1, 1))) == Float32

  # vector constructors
  @test eltype(Vec([1, 1])) == Float64
  @test eltype(Vec([1.0, 1.0])) == Float64
  @test eltype(Vec([1.0f0, 1.0f0])) == Float32
  @test eltype(Vec1([1])) == Float64
  @test eltype(Vec2([1, 1])) == Float64
  @test eltype(Vec3([1, 1, 1])) == Float64
  @test eltype(Vec1f([1])) == Float32
  @test eltype(Vec2f([1, 1])) == Float32
  @test eltype(Vec3f([1, 1, 1])) == Float32

  # parametric constructors 
  @test eltype(Vec{2,T}(1, 1)) == T
  @test eltype(Vec{2,T}((1, 1))) == T
  @test eltype(Vec{2,T}([1, 1])) == T

  # check all 1D Vec constructors, because those tend to make trouble
  @test Vec(1) == Vec((1,)) == Vec([1])
  @test Vec{1,T}(0) == Vec{1,T}((0,)) == Vec{1,T}([0])
  @test Vec{1,T}(-2) == Vec{1,T}((-2,)) == Vec{1,T}([-2])

  # check that input of mixed coordinate types is allowed and works as expected
  @test Vec(1, 0.2) == Vec{2,Float64}(1.0, 0.2)
  @test Vec((3.0, 4)) == Vec{2,Float64}(3.0, 4.0)
  @test Vec([5.0, 6.0, 7]) == Vec{3,Float64}(5.0, 6.0, 7.0)
  @test Vec{2,T}(8, 9.0) == Vec{2,T}((8.0, 9.0))
  @test Vec{2,T}((-1.0, -2)) == Vec{2,T}((-1, -2.0))
  @test Vec{4,T}([0, -1.0, +2, -4.0]) == Vec{4,T}((0.0f0, -1.0f0, +2.0f0, -4.0f0))

  # Integer coordinates converted to Float64
  @test eltype(Vec(1)) == Float64
  @test eltype(Vec(1, 2)) == Float64
  @test eltype(Vec(1, 2, 3)) == Float64
  @test Tuple(Vec(1)) == (1.0,)
  @test Tuple(Vec(1, 2)) == (1.0, 2.0)
  @test Tuple(Vec(1, 2, 3)) == (1.0, 2.0, 3.0)

  # throws
  @test_throws DimensionMismatch Vec{2,T}(1)
  @test_throws DimensionMismatch Vec{3,T}((2, 3))
  @test_throws DimensionMismatch Vec{-3,T}([4, 5, 6])
end
