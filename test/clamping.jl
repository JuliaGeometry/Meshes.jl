@testset "Clamping" begin
  box = Box((zero(T), zero(T)), (one(T), one(T)))
  @test clamp(P2(0.5, 0.5), box) == P2(0.5, 0.5)
  @test clamp(P2(-1, 0.5), box) == P2(0, 0.5)
  @test clamp(P2(0.5, -1), box) == P2(0.5, 0)
  @test clamp(P2(2, 0.5), box) == P2(1, 0.5)
  @test clamp(P2(0.5, 2), box) == P2(0.5, 1)
  @test clamp(P2(2, 2), box) == P2(1, 1)
  @test clamp(P2(-1, -1), box) == P2(0, 0)
end