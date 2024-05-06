@testset "tolerances" begin
  ℒ = ℳ
  𝒜 = typeof(zero(ℳ)^2)
  𝒱 = typeof(zero(ℳ)^3)
  if T === Float32
    @test atol(T) == 1.0f-5
    @test atol(ℒ) == 1.0f-5 * u"m"
    @test atol(𝒜) == 1.0f-5^1.3f0 * u"m^2"
    @test atol(𝒱) == 1.0f-5^3 * u"m^3"
  else
    @test atol(T) == 1e-10
    @test atol(ℒ) == 1e-10 * u"m"
    @test atol(𝒜) == 1e-10^1.3 * u"m^2"
    @test atol(𝒱) == 1e-10^3 * u"m^3"
  end
  @test atol(zero(T)) == atol(T)
  @test atol(zero(ℒ)) == atol(ℒ)
  @test atol(zero(𝒜)) == atol(𝒜)
  @test atol(zero(𝒱)) == atol(𝒱)
end
