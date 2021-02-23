@testset "Primitives" begin
  @testset "Boxes" begin
    b = Box(P2(0,0), P2(1,1))
    @test embeddim(b) == 2
    @test paramdim(b) == 2
    @test coordtype(b) == T
    @test minimum(b) == P2(0,0)
    @test maximum(b) == P2(1,1)
    @test extrema(b) == (P2(0,0), P2(1,1))
    @test measure(b) == T(1)
    @test P2(1,1) ∈ b

    b = Box(P2(1,1), P2(2,2))
    @test sides(b) == T[1,1]
    @test Meshes.center(b) == P2(1.5,1.5)
    @test diagonal(b) == √T(2)

    # intersection of boxes
    b1 = Box(P2(0,0), P2(1,1))
    b2 = Box(P2(0.5,0.5), P2(2,2))
    b3 = Box(P2(2,2), P2(3,3))
    b4 = Box(P2(1,1), P2(2,2))
    b5 = Box(P2(1.0,0.5), P2(2,2))
    @test intersecttype(b1, b2) isa OverlappingBoxes
    @test b1 ∩ b2 == Box(P2(0.5,0.5), P2(1,1))
    @test intersecttype(b1, b3) isa NonIntersectingBoxes
    @test b1 ∩ b3 === nothing
    @test intersecttype(b1, b4) isa CornerTouchingBoxes
    @test b1 ∩ b4 == P2(1,1)
    @test intersecttype(b1, b5) isa FaceTouchingBoxes
    @test b1 ∩ b5 == Box(P2(1.0,0.5), P2(1,1))
  end

  @testset "Ball" begin
    b = Ball(P3(1,2,3), T(5))
    @test embeddim(b) == 3
    @test paramdim(b) == 3
    @test coordtype(b) == T
    @test Meshes.center(b) == P3(1,2,3)
    @test radius(b) == T(5)

    b = Ball(P2(0,0), T(2))
    @test measure(b) ≈ π*(2^2)
    b = Ball(P3(0,0,0), T(2))
    @test measure(b) ≈ (4/3)*π*(2^3)

    b = Ball(P2(0,0), T(2))
    @test P2(1,0) ∈ b
    @test P2(0,1) ∈ b
    @test P2(2,0) ∈ b
    @test P2(0,2) ∈ b
    @test P2(3,5) ∉ b

    b = Ball(P3(0,0,0), T(2))
    @test P3(1,0,0) ∈ b
    @test P3(0,0,1) ∈ b
    @test P3(2,0,0) ∈ b
    @test P3(0,0,2) ∈ b
    @test P3(3,5,2) ∉ b
  end

  @testset "Sphere" begin
    s = Sphere(P3(0,0,0), T(1))
    @test embeddim(s) == 3
    @test paramdim(s) == 2
    @test coordtype(s) == T
    @test Meshes.center(s) == P3(0, 0, 0)
    @test radius(s) == T(1)

    s = Sphere(P2(0,0), T(2))
    @test measure(s) ≈ 2π*2
    s = Sphere(P3(0,0,0), T(2))
    @test measure(s) ≈ 4π*(2^2)

    s = Sphere(P2(0,0), T(2))
    @test P2(1,0) ∉ s
    @test P2(0,1) ∉ s
    @test P2(2,0) ∈ s
    @test P2(0,2) ∈ s
    @test P2(3,5) ∉ s

    s = Sphere(P3(0,0,0), T(2))
    @test P3(1,0,0) ∉ s
    @test P3(0,0,1) ∉ s
    @test P3(2,0,0) ∈ s
    @test P3(0,0,2) ∈ s
    @test P3(3,5,2) ∉ s
  end

  @testset "Cylinder" begin
    c = Cylinder(P3(1,2,3), P3(4,5,6), T(5))
    @test embeddim(c) == 3
    @test paramdim(c) == 3
    @test coordtype(c) == T
    @test radius(c) == T(5)
    @test height(c) ≈ √27

    @test measure(c) ≈ π*5.0^2*√27
  end
end
