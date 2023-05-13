@testset "Support function" begin
  t = Triangle(P2[(0,0), (1,0), (0,1)])
  @test supportfun(t, V2( 1,  0)) == P2(1, 0)
  @test supportfun(t, V2( 0,  1)) == P2(0, 1)
  @test supportfun(t, V2(-1, -1)) == P2(0, 0)
  @test supportfun(t, V2(-1,  1)) == P2(0, 1)
  
  b = Ball(P2(0,0), T(2))
  @test supportfun(b, V2( 1, 1)) ≈ P2(√2, √2)
  @test supportfun(b, V2( 1, 0)) ≈ P2(2, 0)
  @test supportfun(b, V2( 0, 1)) ≈ P2(0, 2)
  @test supportfun(b, V2(-1, 1)) ≈ P2(-√2, √2)
end
