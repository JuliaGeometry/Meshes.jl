@testitem "Tolerances" setup = [Setup] begin
  # absolute tolerance
  ℒ = ℳ
  𝒜 = typeof(zero(ℳ)^2)
  𝒱 = typeof(zero(ℳ)^3)
  @test Meshes.atol(T) == eps(T)^(3 // 4)
  @test Meshes.atol(ℒ) == Meshes.atol(T) * u"m"
  @test Meshes.atol(𝒜) == Meshes.atol(T) * u"m^2"
  @test Meshes.atol(𝒱) == Meshes.atol(T) * u"m^3"
  @test Meshes.atol(zero(T)) == Meshes.atol(T)
  @test Meshes.atol(zero(ℒ)) == Meshes.atol(ℒ)
  @test Meshes.atol(zero(𝒜)) == Meshes.atol(𝒜)
  @test Meshes.atol(zero(𝒱)) == Meshes.atol(𝒱)
  @inferred Meshes.atol(T)
  @inferred Meshes.atol(ℒ)
  @inferred Meshes.atol(𝒜)
  @inferred Meshes.atol(𝒱)

  # relative tolerance
  @test Meshes.rtol(T) == eps(T)^(1 // 3)
  @test Meshes.rtol(ℒ) == Meshes.rtol(T)
  @test Meshes.rtol(𝒜) == Meshes.rtol(T)
  @test Meshes.rtol(𝒱) == Meshes.rtol(T)
  @test Meshes.rtol(zero(T)) == Meshes.rtol(T)
  @test Meshes.rtol(zero(ℒ)) == Meshes.rtol(ℒ)
  @test Meshes.rtol(zero(𝒜)) == Meshes.rtol(𝒜)
  @test Meshes.rtol(zero(𝒱)) == Meshes.rtol(𝒱)
  @inferred Meshes.rtol(T)
  @inferred Meshes.rtol(ℒ)
  @inferred Meshes.rtol(𝒜)
  @inferred Meshes.rtol(𝒱)

  # generic AbstractFloat fallback (e.g. BigFloat)
  @test Meshes.atol(BigFloat) == eps(BigFloat)^(3 // 4)
  @test Meshes.rtol(BigFloat) == eps(BigFloat)^(1 // 3)
  @test Meshes.atol(zero(BigFloat)) == Meshes.atol(BigFloat)
  @test Meshes.rtol(zero(BigFloat)) == Meshes.rtol(BigFloat)

  # maximum length for discretization
  @test Meshes.maxlen() == 500u"km"
end
