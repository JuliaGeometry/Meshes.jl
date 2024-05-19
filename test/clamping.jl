@testset "Clamping" begin
  box = Box((zero(T), zero(T)), (one(T), one(T)))
  @test clamp(point(0.5, 0.5), box) == point(0.5, 0.5)
  @test clamp(point(-1, 0.5), box) == point(0, 0.5)
  @test clamp(point(0.5, -1), box) == point(0.5, 0)
  @test clamp(point(2, 0.5), box) == point(1, 0.5)
  @test clamp(point(0.5, 2), box) == point(0.5, 1)
  @test clamp(point(2, 2), box) == point(1, 1)
  @test clamp(point(-1, -1), box) == point(0, 0)

  points = PointSet(point(0.5, 0.5), point(-1, 0.5), point(0.5, 2))
  @test clamp(points, box) == PointSet(point(0.5, 0.5), point(0, 0.5), point(0.5, 1))
end
