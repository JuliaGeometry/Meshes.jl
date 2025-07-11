@testitem "Tolerances" setup = [Setup] begin
  ℒ = ℳ
  𝒜 = typeof(zero(ℳ)^2)
  𝒱 = typeof(zero(ℳ)^3)

  @test atol(T) == eps(T) ^ (3 // 4)
  @test atol(ℒ) == atol(T) * u"m"
  @test atol(𝒜) == atol(T)^2 * u"m^2"
  @test atol(𝒱) == atol(T)^3 * u"m^3"

  @test atol(zero(T)) == atol(T)
  @test atol(zero(ℒ)) == atol(ℒ)
  @test atol(zero(𝒜)) == atol(𝒜)
  @test atol(zero(𝒱)) == atol(𝒱)

  @inferred atol(T)
  @inferred atol(ℒ)
  @inferred atol(𝒜)
  @inferred atol(𝒱)
end
