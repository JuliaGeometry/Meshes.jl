@testset "Vectors" begin
  Met{T} = Quantity{T,u"𝐋",typeof(u"m")}

  # vararg constructors
  @test eltype(Vec(1, 1)) == Met{Float64}
  @test eltype(Vec(1.0, 1.0)) == Met{Float64}
  @test eltype(Vec(1.0f0, 1.0f0)) == Met{Float32}

  # tuple constructors
  @test eltype(Vec((1, 1))) == Met{Float64}
  @test eltype(Vec((1.0, 1.0))) == Met{Float64}
  @test eltype(Vec((1.0f0, 1.0f0))) == Met{Float32}

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
  @test eltype(Vec(1)) == Met{Float64}
  @test eltype(Vec(1, 2)) == Met{Float64}
  @test eltype(Vec(1, 2, 3)) == Met{Float64}
  @test Tuple(Vec(1)) == (1.0u"m",)
  @test Tuple(Vec(1, 2)) == (1.0u"m", 2.0u"m")
  @test Tuple(Vec(1, 2, 3)) == (1.0u"m", 2.0u"m", 3.0u"m")

  # Unitful coordinates
  vector = Vec(1u"m", 1u"m")
  @test unit(eltype(vector)) == u"m"
  @test Unitful.numtype(eltype(vector)) === Float64
  vector = Vec(1.0u"m", 1.0u"m")
  @test unit(eltype(vector)) == u"m"
  @test Unitful.numtype(eltype(vector)) === Float64
  vector = Vec(1.0f0u"m", 1.0f0u"m")
  @test unit(eltype(vector)) == u"m"
  @test Unitful.numtype(eltype(vector)) === Float32

  # angles between 2D vectors
  @test ∠(vec(1, 0), vec(0, 1)) ≈ T(π / 2)
  @test ∠(vec(1, 0), vec(0, -1)) ≈ T(-π / 2)
  @test ∠(vec(1, 0), vec(-1, 0)) ≈ T(π)
  @test ∠(vec(0, 1), vec(-1, 0)) ≈ T(π / 2)
  @test ∠(vec(0, 1), vec(0, -1)) ≈ T(π)
  @test ∠(vec(0, 1), vec(1, 1)) ≈ T(-π / 4)
  @test ∠(vec(0, -1), vec(1, 1)) ≈ T(π * 3 / 4)
  @test ∠(vec(-1, -1), vec(1, 1)) ≈ T(π)
  @test ∠(vec(-2, 0), vec(2, 0)) ≈ T(π)

  # angles between 3D vectors
  @test ∠(vec(0, 0, 1), vec(1, 1, 0)) ≈ T(π / 2)
  @test ∠(vec(1, 0, 1), vec(1, 1, 0)) ≈ T(π / 3)
  @test ∠(vec(-1, -1, 0), vec(1, 1, 0)) ≈ T(π)
  @test ∠(vec(0, -1, -1), vec(0, 1, 1)) ≈ T(π)
  @test ∠(vec(0, -1, -1), vec(0, 1, 0)) ≈ T(π * 3 / 4)
  @test ∠(vec(0, 1, 1), vec(1, 1, 0)) ≈ T(π / 3)

  v = vec(0, 1)
  @test sprint(show, v, context=:compact => true) == "(0.0 m, 1.0 m)"
  if T === Float32
    @test sprint(show, v) == "Vec(0.0f0 m, 1.0f0 m)"
  else
    @test sprint(show, v) == "Vec(0.0 m, 1.0 m)"
  end
  @test sprint(show, MIME("text/plain"), v) == sprint(show, v)
end
