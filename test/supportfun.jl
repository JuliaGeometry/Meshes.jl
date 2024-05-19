@testset "Support function" begin
  t = Triangle(point(0, 0), point(1, 0), point(0, 1))
  @test supportfun(t, vector(1, 0)) == point(1, 0)
  @test supportfun(t, vector(0, 1)) == point(0, 1)
  @test supportfun(t, vector(-1, -1)) == point(0, 0)
  @test supportfun(t, vector(-1, 1)) == point(0, 1)

  b = Ball(point(0, 0), T(2))
  @test supportfun(b, vector(1, 1)) ≈ point(√2, √2)
  @test supportfun(b, vector(1, 0)) ≈ point(2, 0)
  @test supportfun(b, vector(0, 1)) ≈ point(0, 2)
  @test supportfun(b, vector(-1, 1)) ≈ point(-√2, √2)

  b = Box(point(0, 0), point(1, 1))
  @test supportfun(b, vector(1, 1)) ≈ point(1, 1)
  @test supportfun(b, vector(1, 0)) ≈ point(1, 0)
  @test supportfun(b, vector(-1, 0)) ≈ point(0, 0)
  @test supportfun(b, vector(-1, -1)) ≈ point(0, 0)
  @test supportfun(b, vector(-1, 1)) ≈ point(0, 1)
end
