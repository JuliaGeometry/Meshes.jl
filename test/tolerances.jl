@testitem "Tolerances" setup = [Setup] begin
  ℒ = ℳ
  𝒜 = typeof(zero(ℳ)^2)
  𝒱 = typeof(zero(ℳ)^3)

  @test Meshes.atol(T) == eps(T)^(3 // 4)
  @test Meshes.atol(ℒ) == Meshes.atol(T) * u"m"
  @test Meshes.atol(𝒜) == Meshes.atol(T)^2 * u"m^2"
  @test Meshes.atol(𝒱) == Meshes.atol(T)^3 * u"m^3"

  @test Meshes.atol(zero(T)) == Meshes.atol(T)
  @test Meshes.atol(zero(ℒ)) == Meshes.atol(ℒ)
  @test Meshes.atol(zero(𝒜)) == Meshes.atol(𝒜)
  @test Meshes.atol(zero(𝒱)) == Meshes.atol(𝒱)

  @inferred Meshes.atol(T)
  @inferred Meshes.atol(ℒ)
  @inferred Meshes.atol(𝒜)
  @inferred Meshes.atol(𝒱)
end
