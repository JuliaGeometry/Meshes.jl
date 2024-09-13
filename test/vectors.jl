@testitem "Vectors" setup = [Setup] begin
  # vararg constructors
  @test eltype(Vec(1, 1)) == Float64
  @test eltype(Vec(1.0, 1.0)) == Float64
  @test eltype(Vec(1.0f0, 1.0f0)) == Float32

  # tuple constructors
  @test eltype(Vec((1, 1))) == Float64
  @test eltype(Vec((1.0, 1.0))) == Float64
  @test eltype(Vec((1.0f0, 1.0f0))) == Float32

  # check all 1D Vec constructors, because those tend to make trouble
  @test Vec(T(1)) == Vec((T(1),))
  @test Vec(T(0)) == Vec((T(0),))
  @test Vec(T(-2)) == Vec((T(-2),))

  # check that input of mixed coordinate types is allowed and works as expected
  @test Vec(1, 0.2) == Vec(1.0, 0.2)
  @test Vec((3.0, 4)) == Vec(3.0, 4.0)
  @test Vec((5.0, 6.0, 7)) == Vec(5.0, 6.0, 7.0)
  @test Vec(8, T(9.0)) == Vec((T(8.0), T(9.0)))
  @test Vec((T(-1.0), -2)) == Vec((T(-1.0), T(-2.0)))
  @test Vec((0, T(-1.0), +2, T(-4.0))) == Vec((T(0.0), T(-1.0), T(+2.0), T(-4.0)))

  # integer coordinates are converted to float
  @test eltype(Vec(1)) == Float64
  @test eltype(Vec(1, 2)) == Float64
  @test eltype(Vec(1, 2, 3)) == Float64
  @test Tuple(Vec(1)) == (1.0,)
  @test Tuple(Vec(1, 2)) == (1.0, 2.0)
  @test Tuple(Vec(1, 2, 3)) == (1.0, 2.0, 3.0)

  # Unitful coordinates
  v = Vec(1u"m", 1u"m")
  @test unit(eltype(v)) == u"m"
  @test Unitful.numtype(eltype(v)) === Float64
  v = Vec(1.0u"m", 1.0u"m")
  @test unit(eltype(v)) == u"m"
  @test Unitful.numtype(eltype(v)) === Float64
  v = Vec(1.0f0u"m", 1.0f0u"m")
  @test unit(eltype(v)) == u"m"
  @test Unitful.numtype(eltype(v)) === Float32

  # angles between 2D vectors
  @test ∠(vector(1, 0), vector(0, 1)) ≈ T(π / 2)
  @test ∠(vector(1, 0), vector(0, -1)) ≈ T(-π / 2)
  @test ∠(vector(1, 0), vector(-1, 0)) ≈ T(π)
  @test ∠(vector(0, 1), vector(-1, 0)) ≈ T(π / 2)
  @test ∠(vector(0, 1), vector(0, -1)) ≈ T(π)
  @test ∠(vector(0, 1), vector(1, 1)) ≈ T(-π / 4)
  @test ∠(vector(0, -1), vector(1, 1)) ≈ T(π * 3 / 4)
  @test ∠(vector(-1, -1), vector(1, 1)) ≈ T(π)
  @test ∠(vector(-2, 0), vector(2, 0)) ≈ T(π)

  # angles between 3D vectors
  @test ∠(vector(0, 0, 1), vector(1, 1, 0)) ≈ T(π / 2)
  @test ∠(vector(1, 0, 1), vector(1, 1, 0)) ≈ T(π / 3)
  @test ∠(vector(-1, -1, 0), vector(1, 1, 0)) ≈ T(π)
  @test ∠(vector(0, -1, -1), vector(0, 1, 1)) ≈ T(π)
  @test ∠(vector(0, -1, -1), vector(0, 1, 0)) ≈ T(π * 3 / 4)
  @test ∠(vector(0, 1, 1), vector(1, 1, 0)) ≈ T(π / 3)

  v = vector(0, 1)
  @test sprint(show, v, context=:compact => true) == "(0.0, 1.0)"
  if T === Float32
    @test sprint(show, v) == "Vec(0.0f0, 1.0f0)"
  else
    @test sprint(show, v) == "Vec(0.0, 1.0)"
  end
  @test sprint(show, MIME("text/plain"), v) == sprint(show, v)
end
