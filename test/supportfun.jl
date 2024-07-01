@testset "Support function" begin
  t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
  @test supportfun(t, vector(1, 0)) == cart(1, 0)
  @test supportfun(t, vector(0, 1)) == cart(0, 1)
  @test supportfun(t, vector(-1, -1)) == cart(0, 0)
  @test supportfun(t, vector(-1, 1)) == cart(0, 1)

  b = Ball(cart(0, 0), T(2))
  @test supportfun(b, vector(1, 1)) ≈ cart(√2, √2)
  @test supportfun(b, vector(1, 0)) ≈ cart(2, 0)
  @test supportfun(b, vector(0, 1)) ≈ cart(0, 2)
  @test supportfun(b, vector(-1, 1)) ≈ cart(-√2, √2)

  b = Box(cart(0, 0), cart(1, 1))
  @test supportfun(b, vector(1, 1)) ≈ cart(1, 1)
  @test supportfun(b, vector(1, 0)) ≈ cart(1, 0)
  @test supportfun(b, vector(-1, 0)) ≈ cart(0, 0)
  @test supportfun(b, vector(-1, -1)) ≈ cart(0, 0)
  @test supportfun(b, vector(-1, 1)) ≈ cart(0, 1)
end
