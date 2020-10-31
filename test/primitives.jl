for T in (Float32, Float64)
P2 = Point{2,T}
P3 = Point{3,T}
@testset "Primitives ($T)" begin
  @testset "Boxes" begin
    b = Box(P2(0, 0), P2(1, 1))
    @test embeddim(b) == 2
    @test coordtype(b) == T
    @test minimum(b) == P2(0, 0)
    @test maximum(b) == P2(1, 1)
    @test extrema(b) == (P2(0, 0), P2(1, 1))

    @test measure(b) == T(1)
    @test P2(1, 1) ∈ b
  end

  @testset "Ball" begin
    b = Ball(P3(1,2,3), T(5))
    @test embeddim(b) == 3
    @test coordtype(b) == T
    @test center(b) == P3(1,2,3)
    @test radius(b) == T(5)

    b = Ball(P2(0,0), T(2))
    @test measure(b) ≈ π*(2^2)
    b = Ball(P3(0,0,0), T(2))
    @test measure(b) ≈ (4/3)*π*(2^3)
  end

  @testset "Sphere" begin
    s = Sphere(P3(0,0,0), T(1))
    @test embeddim(s) == 3
    @test coordtype(s) == T
    @test center(s) == P3(0, 0, 0)
    @test radius(s) == T(1)

    s = Sphere(P2(0,0), T(2))
    @test measure(s) ≈ 2π*2
    s = Sphere(P3(0,0,0), T(2))
    @test measure(s) ≈ 4π*(2^2)
  end

  @testset "Cylinder" begin
    c = Cylinder(P3(1,2,3), P3(4,5,6), T(5))
    @test embeddim(c) == 3
    @test coordtype(c) == T
    @test radius(c) == T(5)
    @test height(c) ≈ √27

    @test measure(c) ≈ π*5.0^2*√27
  end
end
end
