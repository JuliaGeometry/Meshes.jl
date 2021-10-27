@testset "Angles" begin
  # 2D points
  @test ∠(P2(0,1),P2(0,0),P2(1,0)) ≈ T(-π/2)
  @test ∠(P2(1,0),P2(0,0),P2(0,1)) ≈ T(π/2)
  @test ∠(P2(-1,0),P2(0,0),P2(0,1)) ≈ T(-π/2)
  @test ∠(P2(0,1),P2(0,0),P2(-1,0)) ≈ T(π/2)
  @test ∠(P2(0,-1),P2(0,0),P2(1,0)) ≈ T(π/2)
  @test ∠(P2(1,0),P2(0,0),P2(0,-1)) ≈ T(-π/2)
  @test ∠(P2(0,-1),P2(0,0),P2(-1,0)) ≈ T(-π/2)
  @test ∠(P2(-1,0),P2(0,0),P2(0,-1)) ≈ T(π/2)

  # 3D points
  @test ∠(P3(1,0,0),P3(0,0,0),P3(0,1,0)) ≈ T(π/2)
  @test ∠(P3(1,0,0),P3(0,0,0),P3(0,0,1)) ≈ T(π/2)
  @test ∠(P3(0,1,0),P3(0,0,0),P3(1,0,0)) ≈ T(π/2)
  @test ∠(P3(0,1,0),P3(0,0,0),P3(0,0,1)) ≈ T(π/2)
  @test ∠(P3(0,0,1),P3(0,0,0),P3(1,0,0)) ≈ T(π/2)
  @test ∠(P3(0,0,1),P3(0,0,0),P3(0,1,0)) ≈ T(π/2)
  
  # Ngon
  t = Triangle(Point2(0,0), Point2(1,0), Point2(0,1))
  @test all(isapprox.(rad2deg.(angles(t)), [90.0, 45.0, 45.0], atol=8*eps(45.0)))
  @test all(isapprox.(rad2deg.(innerangles(t)), [90.0, 45.0, 45.0], atol=8*eps(45.0)))     
end
