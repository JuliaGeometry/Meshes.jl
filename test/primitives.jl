@testset "Primitives" begin
  @testset "Boxes" begin
    box = Box((0.0, 0.0), (1.0, 1.0))
    @test embeddim(box) == 2
    @test coordtype(box) == Float64
    @test minimum(box) == Point(0.0,0.0)
    @test maximum(box) == Point(1.0,1.0)
    @test extrema(box) == (Point(0.0,0.0), Point(1.0,1.0))
    @test volume(box) == 1.0
    @test Point(1.0, 1.0) ∈ box
  end

  @testset "Ball" begin
    ball = Ball((1,2,3), 5)
    @test embeddim(ball) == 3
    @test coordtype(ball) == Int
    @test center(ball) == Point(1,2,3)
    @test radius(ball) == 5
  end

  @testset "Sphere" begin
    sphere = Sphere(Point3f(0,0,0), 1.0f0)
    @test embeddim(sphere) == 3
    @test coordtype(sphere) == Float32
    @test center(sphere) == Point3f(0,0,0)
    @test radius(sphere) == 1.0f0
  end

  @testset "Cylinder" begin
    c = Cylinder(Point3(1,2,3), Point3(4,5,6), 5.0)
    @test embeddim(c) == 3
    @test coordtype(c) == Float64
    @test radius(c) == 5.0
    @test height(c) == √27
    @test volume(c) == π*5.0^2*√27
  end
end
