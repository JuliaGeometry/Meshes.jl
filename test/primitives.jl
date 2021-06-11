@testset "Primitives" begin
  @testset "Lines" begin
    l = Line(P2(0,0), P2(1,1))
    @test paramdim(l) == 1
    @test isconvex(l)

    l = Line(P2(0,0), P2(1,1))
    @test (l(0), l(1)) == (P2(0,0), P2(1,1))
  end

  @testset "Rays" begin
    r = Ray(P2(0,0), V2(1,1))
    @test paramdim(r) == 1
    @test isconvex(r)

    r = Ray(P2(0,0), V2(1,1))
    @test r(T(0.)) == P2(0,0)
    @test r(T(1.)) == P2(1,1)
    @test r(T(Inf)) == P2(Inf,Inf)
    @test_throws DomainError(T(-1), "r(t) is not defined for t < 0.") r(T(-1))
  end

  @testset "Planes" begin
    p = Plane(P3(0, 0, 0), V3(1, 0, 0), V3(0, 1, 0))
    @test p(T(1.0), T(0.0)) == P3(1, 0, 0)
    @test paramdim(p) == 2
    @test embeddim(p) == 3
    @test isconvex(p)
  end

  @testset "Bezier curves" begin
    # fix import conflict with Plots
    BezierCurve = Meshes.BezierCurve

    b = BezierCurve(P2(0,0),P2(0.5,1),P2(1,0))
    for method in [DeCasteljau(), Horner()]
      @test b(T(0), method) == P2(0,0)
      @test b(T(1), method) == P2(1,0)
      @test b(T(0.5), method) == P2(0.5,0.5)
      @test b(T(0.5), method) == P2(0.5,0.5)
      @test_throws DomainError(T(-0.1), "b(t) is not defined for t outside [0, 1].") b(T(-0.1), method)
      @test_throws DomainError(T(1.2), "b(t) is not defined for t outside [0, 1].") b(T(1.2), method)
    end

    b = BezierCurve(P2.(randn(100), randn(100)))
    t1 = @timed b(T(0.2))
    t2 = @timed b(T(0.2), Horner())
    @test t1.time > t2.time
    @test t2.bytes < 100
  end

  @testset "Boxes" begin
    b = Box(P2(0,0), P2(1,1))
    @test embeddim(b) == 2
    @test paramdim(b) == 2
    @test coordtype(b) == T
    @test minimum(b) == P2(0,0)
    @test maximum(b) == P2(1,1)
    @test extrema(b) == (P2(0,0), P2(1,1))
    @test isconvex(b)

    b = Box(P2(0,0), P2(1,1))
    @test measure(b) == T(1)
    @test P2(1,1) ∈ b

    b = Box(P2(1,1), P2(2,2))
    @test sides(b) == T[1,1]
    @test Meshes.center(b) == P2(1.5,1.5)
    @test diagonal(b) == √T(2)

    b = Box(P2(1,2), P2(3,4))
    @test vertices(b) == P2[(1,2),(3,2),(3,4),(1,4)]

    b = Box(P3(1,2,3), P3(4,5,6))
    @test vertices(b) == P3[(1,2,3),(4,2,3),(4,5,3),(1,5,3),(1,2,6),(4,2,6),(4,5,6),(1,5,6)]

    # subsetting with boxes
    b1 = Box(P2(0,0), P2(0.5,0.5))
    b2 = Box(P2(0.1,0.1), P2(0.5,0.5))
    b3 = Box(P2(0,0), P2(1,1))
    @test b1 ⊆ b3
    @test b2 ⊆ b3
    @test !(b1 ⊆ b2)
    @test !(b3 ⊆ b1)
    @test !(b3 ⊆ b1)
  end

  @testset "Ball" begin
    b = Ball(P3(1,2,3), T(5))
    @test embeddim(b) == 3
    @test paramdim(b) == 3
    @test coordtype(b) == T
    @test Meshes.center(b) == P3(1,2,3)
    @test radius(b) == T(5)
    @test isconvex(b)

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
    @test !isconvex(s)

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

    # 2D sphere passing through 3 points
    s = Sphere(P2(0,0), P2(0.5,0), P2(1,1))
    @test Meshes.center(s) == P2(0.25, 0.75)
    @test radius(s) == T(0.7905694150420949)
    s = Sphere(P2(0,0), P2(1,0), P2(0,1))
    @test Meshes.center(s) == P2(0.5, 0.5)
    @test radius(s) == T(0.7071067811865476)
    s = Sphere(P2(0,0), P2(1,0), P2(1,1))
    @test Meshes.center(s) == P2(0.5, 0.5)
    @test radius(s) == T(0.7071067811865476)
  end

  @testset "Cylinder" begin
    c = Cylinder(P3(1,2,3), P3(4,5,6), T(5))
    @test embeddim(c) == 3
    @test paramdim(c) == 3
    @test coordtype(c) == T
    @test radius(c) == T(5)
    @test height(c) ≈ √27
    @test isconvex(c)

    @test measure(c) ≈ π*5.0^2*√27
  end
end
