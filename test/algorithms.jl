@testset "Algorithms" begin
  @testset "Intersection" begin
    a = Segment(Point(0.0, 0.0), Point(4.0, 1.0))
    b = Segment(Point(0.0, 0.25), Point(3.0, 0.25))
    c = Segment(Point(0.0, 0.25), Point(0.5, 0.25))
    d = Segment(Point(0.0, 0.0), Point(0.0, 4.0))
    e = Segment(Point(1.0, 0.0), Point(0.0, 4.0))
    f = Segment(Point(5.0, 0.0), Point(6.0, 0.0))
    @test a ∩ b === (true, Point(1.0, 0.25))
    @test a ∩ c === (false, Point(0.0, 0.0))
    @test d ∩ d === (false, Point(0.0, 0.0))
    found, point = d ∩ e
    @test found && coordinates(point) ≈ [0.0, 4.0]
    @test a ∩ f === (false, Point(0.0, 0.0))
  end
end
