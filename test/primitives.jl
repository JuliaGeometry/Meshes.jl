@testset "Primitives" begin
  @testset "Boxes" begin
    b = Box((0.0, 0.0), (1.0, 1.0))
    @test embeddim(b) == 2
    @test coordtype(b) == Float64
    @test minimum(b) == Point(0.0,0.0)
    @test maximum(b) == Point(1.0,1.0)
    @test extrema(b) == (Point(0.0,0.0), Point(1.0,1.0))

    b = Box((0.0, 0.0), (1.0, 1.0))
    @test measure(b) == 1.0
    @test Point(1.0, 1.0) ∈ b
  end

  @testset "Ball" begin
    b = Ball((1,2,3), 5)
    @test embeddim(b) == 3
    @test coordtype(b) == Int
    @test center(b) == Point(1,2,3)
    @test radius(b) == 5

    b = Ball((0,0), 2)
    @test measure(b) ≈ π*(2^2)
    b = Ball((0,0,0), 2)
    @test measure(b) ≈ (4/3)*π*(2^3)
  end

  @testset "Sphere" begin
    s = Sphere(Point3(0,0,0), 1.0)
    @test embeddim(s) == 3
    @test coordtype(s) == Float64
    @test center(s) == Point3(0, 0, 0)
    @test radius(s) == 1.0

    s = Sphere((0,0), 2)
    @test measure(s) == 2π*2
    s = Sphere((0,0,0), 2)
    @test measure(s) == 4π*(2^2)
  end

  @testset "Cylinder" begin
    c = Cylinder(Point3(1,2,3), Point3(4,5,6), 5.0)
    @test embeddim(c) == 3
    @test coordtype(c) == Float64
    @test radius(c) == 5.0
    @test height(c) == √27

    c = Cylinder(Point3(1,2,3), Point3(4,5,6), 5.0)
    @test measure(c) == π*5.0^2*√27
  end
end
