@testset "tolerances" begin
  Q = typeof(zero(T) * u"m")
  if T === Float32
    @test atol(T) == 1.0f-5
    @test atol(Q) == 1.0f-5 * u"m"
  else
    @test atol(T) == 1e-10
    @test atol(Q) == 1e-10 * u"m"
  end
end
