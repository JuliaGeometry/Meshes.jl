@testset "Clamping" begin
  box = Box((zero(T), zero(T)), (one(T), one(T)))
  @test clamp(cart(0.5, 0.5), box) == cart(0.5, 0.5)
  @test clamp(cart(-1, 0.5), box) == cart(0, 0.5)
  @test clamp(cart(0.5, -1), box) == cart(0.5, 0)
  @test clamp(cart(2, 0.5), box) == cart(1, 0.5)
  @test clamp(cart(0.5, 2), box) == cart(0.5, 1)
  @test clamp(cart(2, 2), box) == cart(1, 1)
  @test clamp(cart(-1, -1), box) == cart(0, 0)

  points = PointSet(cart(0.5, 0.5), cart(-1, 0.5), cart(0.5, 2))
  @test clamp(points, box) == PointSet(cart(0.5, 0.5), cart(0, 0.5), cart(0.5, 1))
end
