@testset "Primitives" begin
  @testset "Lines" begin
    l = Line(P2(0,0), P2(1,1))
    @test paramdim(l) == 1
    @test isconvex(l)
    @test measure(l) == T(Inf)
    @test length(l) == T(Inf)
    @test isnothing(boundary(l))
    @test perimeter(l) == zero(T)

    l = Line(P2(0,0), P2(1,1))
    @test (l(0), l(1)) == (P2(0,0), P2(1,1))
  end

  @testset "Rays" begin
    r = Ray(P2(0,0), V2(1,1))
    @test paramdim(r) == 1
    @test isconvex(r)
    @test measure(r) == T(Inf)
    @test length(r) == T(Inf)
    @test origin(r) == P2(0,0)
    @test direction(r) == V2(1,1)
    @test boundary(r) == P2(0,0)
    @test perimeter(r) == zero(T)

    r = Ray(P2(0,0), V2(1,1))
    @test r(T(0.)) == P2(0,0)
    @test r(T(1.)) == P2(1,1)
    @test r(T(Inf)) == P2(Inf,Inf)
    @test_throws DomainError(T(-1), "r(t) is not defined for t < 0.") r(T(-1))

    p₁ = P3(3,3,3)
    p₂ = P3(-3,-3,-3)
    p₃ = P3(1,0,0)
    r = Ray(P3(0,0,0), V3(1,1,1))
    @test p₁ ∈ r
    @test p₂ ∉ r
    @test p₃ ∉ r

    r1 = Ray(P3(0,0,0), V3(1,0,0))
    r2 = Ray(P3(1,1,1), V3(1,2,1))
    @test r1 != r2

    r1 = Ray(P3(0,0,0), V3(1,0,0))
    r2 = Ray(P3(1,0,0), V3(-1,0,0))
    @test r1 != r2

    r1 = Ray(P3(0,0,0), V3(1,0,0))
    r2 = Ray(P3(1,0,0), V3(1,0,0))
    @test r1 != r2

    r1 = Ray(P3(0,0,0), V3(2,0,0))
    r2 = Ray(P3(0,0,0), V3(1,0,0))
    @test r1 == r2
  end

  @testset "Planes" begin
    p = Plane(P3(0,0,0), V3(1,0,0), V3(0,1,0))
    @test p(T(1), T(0)) == P3(1,0,0)
    @test paramdim(p) == 2
    @test embeddim(p) == 3
    @test isconvex(p)
    @test measure(p) == T(Inf)
    @test area(p) == T(Inf)
    @test p(T(0), T(0)) == P3(0,0,0)
    @test normal(p) == Vec(0,0,1)
    @test isnothing(boundary(p))
    @test perimeter(p) == zero(T)

    p = Plane(P3(0,0,0), V3(0,0,1))
    @test p(T(1), T(0)) == P3(1,0,0)
    @test p(T(0), T(1)) == P3(0,1,0)

    p₁ = Plane(P3(0,0,0), V3(1,0,0), V3(0,1,0))
    p₂ = Plane(P3(0,0,0), V3(0,1,0), V3(1,0,0))
    @test p₁ == p₂
    p₁ = Plane(P3(0,0,0),  V3(1,1,0))
    p₂ = Plane(P3(0,0,0), -V3(1,1,0))
    @test p₁ == p₂

    # normal to plane has norm one regardless of basis
    p = Plane(P3(0,0,0), V3(2,0,0), V3(0,3,0))
    n = normal(p)
    @test isapprox(norm(n), T(1), atol=atol(T))

    # plane passing through three points
    p₁ = P3(0,0,0)
    p₂ = P3(1,2,3)
    p₃ = P3(3,2,1)
    p = Plane(p₁, p₂, p₃)
    @test p₁ ∈ p
    @test p₂ ∈ p
    @test p₃ ∈ p
  end

  @testset "BezierCurve" begin
    # fix import conflict with Plots
    BezierCurve = Meshes.BezierCurve

    b = BezierCurve(P2(0,0),P2(0.5,1),P2(1,0))
    @test embeddim(b) == 2
    @test paramdim(b) == 1

    b = BezierCurve(P2(0,0),P2(0.5,1),P2(1,0))
    for method in [DeCasteljau(), Horner()]
      @test b(T(0), method) == P2(0,0)
      @test b(T(1), method) == P2(1,0)
      @test b(T(0.5), method) == P2(0.5,0.5)
      @test b(T(0.5), method) == P2(0.5,0.5)
      @test_throws DomainError(T(-0.1), "b(t) is not defined for t outside [0, 1].") b(T(-0.1), method)
      @test_throws DomainError(T(1.2), "b(t) is not defined for t outside [0, 1].") b(T(1.2), method)
    end

    @test boundary(b) == PointSet(P2(0,0), P2(1,0))
    b = BezierCurve(P2(0,0), P2(1,1))
    @test boundary(b) == PointSet([P2(0,0), P2(1,1)])
    @test perimeter(b) == zero(T)

    b = BezierCurve(P2.(randn(100), randn(100)))
    t1 = @timed b(T(0.2))
    t2 = @timed b(T(0.2), Horner())
    @test t1.time > t2.time
    @test t2.bytes < 100
  end

  @testset "Boxes" begin
    b = Box(P1(0), P1(1))
    @test embeddim(b) == 1
    @test paramdim(b) == 1
    @test coordtype(b) == T
    @test minimum(b) == P1(0)
    @test maximum(b) == P1(1)
    @test extrema(b) == (P1(0), P1(1))
    @test isconvex(b)

    b = Box(P1(0), P1(1))
    @test vertices(b) == [P1(0), P1(1)]
    @test measure(b) == T(1)
    @test P1(0) ∈ b
    @test P1(1) ∈ b
    @test P1(0.5) ∈ b
    @test P1(-0.5) ∉ b
    @test P1(1.5) ∉ b

    b = Box(P2(0,0), P2(1,1))
    @test embeddim(b) == 2
    @test paramdim(b) == 2
    @test coordtype(b) == T
    @test minimum(b) == P2(0,0)
    @test maximum(b) == P2(1,1)
    @test extrema(b) == (P2(0,0), P2(1,1))
    @test isconvex(b)

    b = Box(P2(0,0), P2(1,1))
    @test measure(b) == area(b) == T(1)
    @test P2(1,1) ∈ b
    @test perimeter(b) ≈ T(4)

    b = Box(P2(1,1), P2(2,2))
    @test sides(b) == T.((1,1))
    @test Meshes.center(b) == P2(1.5,1.5)
    @test diagonal(b) == √T(2)

    b = Box(P2(1,2), P2(3,4))
    @test vertices(b) == P2[(1,2),(3,2),(3,4),(1,4)]

    b = Box(P3(1,2,3), P3(4,5,6))
    @test vertices(b) == P3[(1,2,3),(4,2,3),(4,5,3),(1,5,3),(1,2,6),(4,2,6),(4,5,6),(1,5,6)]

    b = Box(P2(0,0), P2(1,1))
    @test boundary(b) == Chain(P2[(0,0),(1,0),(1,1),(0,1),(0,0)])

    b = Box(P3(0,0,0), P3(1,1,1))
    m = boundary(b)
    @test m isa Mesh
    @test nvertices(m) == 8
    @test nelements(m) == 6

    # subsetting with boxes
    b1 = Box(P2(0,0), P2(0.5,0.5))
    b2 = Box(P2(0.1,0.1), P2(0.5,0.5))
    b3 = Box(P2(0,0), P2(1,1))
    @test b1 ⊆ b3
    @test b2 ⊆ b3
    @test !(b1 ⊆ b2)
    @test !(b3 ⊆ b1)
    @test !(b3 ⊆ b1)

    b = Box(P2(0,0), P2(10,20))
    @test b(T(0.0), T(0.0)) == P2(0,0)
    @test b(T(0.5), T(0.0)) == P2(5,0)
    @test b(T(1.0), T(0.0)) == P2(10,0)
    @test b(T(0.0), T(0.5)) == P2(0,10)
    @test b(T(0.0), T(1.0)) == P2(0,20)

    b = Box(P3(0,0,0), P3(10,20,30))
    @test b(T(0.0), T(0.0), T(0.0)) == P3(0,0,0)
    @test b(T(1.0), T(1.0), T(1.0)) == P3(10,20,30)
  end

  @testset "Ball" begin
    b = Ball(P3(1,2,3), T(5))
    @test embeddim(b) == 3
    @test paramdim(b) == 3
    @test coordtype(b) == T
    @test Meshes.center(b) == P3(1,2,3)
    @test radius(b) == T(5)
    @test isconvex(b)

    b = Ball(P3(1,2,3), 4)
    @test coordtype(b) == T

    b = Ball(P2(0,0), T(2))
    @test measure(b) ≈ T(π)*(T(2)^2)
    b = Ball(P3(0,0,0), T(2))
    @test measure(b) ≈ T(4/3)*T(π)*(T(2)^3)

    b = Ball(P2(0,0), T(2))
    @test P2(1,0) ∈ b
    @test P2(0,1) ∈ b
    @test P2(2,0) ∈ b
    @test P2(0,2) ∈ b
    @test P2(3,5) ∉ b
    @test perimeter(b) ≈ T(4π)

    b = Ball(P3(0,0,0), T(2))
    @test P3(1,0,0) ∈ b
    @test P3(0,0,1) ∈ b
    @test P3(2,0,0) ∈ b
    @test P3(0,0,2) ∈ b
    @test P3(3,5,2) ∉ b

    b = Ball(P2(0,0), T(2))
    @test b(T(0), T(0)) ≈ P2(0, 0)
    @test b(T(1), T(0)) ≈ P2(2, 0)

    b = Ball(P3(0,0,0), T(2))
    @test b(T(0), T(0), T(0)) ≈ P3(0, 0, 0)
    @test b(T(1), T(0), T(0)) ≈ P3(0, 0, 2)
  end

  @testset "Sphere" begin
    s = Sphere(P3(0,0,0), T(1))
    @test embeddim(s) == 3
    @test paramdim(s) == 2
    @test coordtype(s) == T
    @test Meshes.center(s) == P3(0, 0, 0)
    @test radius(s) == T(1)
    @test extrema(s) == (P3(-1,-1,-1), P3(1,1,1))
    @test !isconvex(s)
    @test isnothing(boundary(s))
    @test isperiodic(s) == (true, true)
    @test perimeter(s) == zero(T)

    s = Sphere(P3(1,2,3), 4)
    @test coordtype(s) == T

    s = Sphere(P2(0,0), T(1))
    @test embeddim(s) == 2
    @test paramdim(s) == 1
    @test coordtype(s) == T
    @test Meshes.center(s) == P2(0, 0)
    @test radius(s) == T(1)
    @test extrema(s) == (P2(-1,-1), P2(1,1))
    @test !isconvex(s)
    @test isnothing(boundary(s))
    @test isperiodic(s) == (true,)

    s = Sphere(P2(0,0), T(2))
    @test measure(s) ≈ 2π*2
    @test length(s) ≈ 2π*2
    @test extrema(s) == (P2(-2,-2), P2(2,2))
    s = Sphere(P3(0,0,0), T(2))
    @test measure(s) ≈ 4π*(2^2)
    @test area(s) ≈ 4π*(2^2)

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

    # 3D sphere passing through 4 points
    s = Sphere(P3(0,0,0), P3(5,0,1), P3(1,1,1), P3(3,2,1))
    @test P3(0,0,0) ∈ s
    @test P3(5,0,1) ∈ s
    @test P3(1,1,1) ∈ s
    @test P3(3,2,1) ∈ s
    O = Meshes.center(s)
    r = radius(s)
    @test isapprox(r, norm(P3(0,0,0) - O))    

    s = Sphere(P2(0,0), T(2))
    @test s(T(0)) ≈ P2(2, 0)
    @test s(T(0.5)) ≈ P2(-2, 0)

    s = Sphere(P3(0,0,0), T(2))
    @test s(T(0), T(0)) ≈ P3(0, 0, 2)
    @test s(T(0.5), T(0.5)) ≈ P3(0, 0, -2)
  end

  @testset "Disk" begin
    p = Plane(P3(0,0,0), V3(0,0,1))
    d = Disk(p, T(2))
    @test embeddim(d) == 3
    @test paramdim(d) == 2
    @test coordtype(d) == T
    @test Meshes.center(d) == P3(0,0,0)
    @test radius(d) == T(2)
    @test isconvex(d)
    @test measure(d) == T(π)*T(2)^2
    @test area(d) == measure(d)
    @test P3(0,0,0) ∈ d
    @test P3(0,0,1) ∉ d
    @test boundary(d) == Circle(p, T(2))
  end

  @testset "Circle" begin
    p = Plane(P3(0,0,0), V3(0,0,1))
    c = Circle(p, T(2))
    @test embeddim(c) == 3
    @test paramdim(c) == 1
    @test coordtype(c) == T
    @test Meshes.center(c) == P3(0,0,0)
    @test radius(c) == T(2)
    @test !isconvex(c)
    @test measure(c) == 2*T(π)*T(2)
    @test length(c) == measure(c)
    @test P3(2,0,0) ∈ c
    @test P3(0,2,0) ∈ c
    @test P3(0,0,0) ∉ c
    @test isnothing(boundary(c))

    # 3D circumcircle
    p1 = P3(0,4,0)
    p2 = P3(0,-4,0)
    p3 = P3(0,0,4)
    c = Circle(p1, p2, p3)
    @test p1 ∈ c
    @test p2 ∈ c
    @test p3 ∈ c
  end

  @testset "Cylinder" begin
    c = Cylinder(T(5),
                 Plane(P3(1,2,3), V3(0,0,1)),
                 Plane(P3(4,5,6), V3(0,0,1)))
    @test embeddim(c) == 3
    @test paramdim(c) == 3
    @test coordtype(c) == T
    @test radius(c) == T(5)
    @test bottom(c) == Plane(P3(1,2,3), V3(0,0,1))
    @test top(c) == Plane(P3(4,5,6), V3(0,0,1))
    @test axis(c) == Line(P3(1,2,3), P3(4,5,6))
    @test isconvex(c)
    @test !isright(c)
    @test measure(c) == volume(c) ≈ T(5)^2 * pi * T(3)*sqrt(T(3))
    @test P3(1,2,3) ∈ c
    @test P3(4,5,6) ∈ c
    @test P3(0.99,1.99,2.99) ∉ c
    @test P3(4.01,5.01,6.01) ∉ c

    c = Cylinder(1.0)
    @test coordtype(c) == Float64
    c = Cylinder(1f0)
    @test coordtype(c) == Float32
    c = Cylinder(1)
    @test coordtype(c) == Float64

    c = Cylinder(T(1), Segment(P3(0,0,0), P3(0,0,1)))
    @test radius(c) == T(1)
    @test bottom(c) == Plane(P3(0,0,0), V3(0,0,1))
    @test top(c) == Plane(P3(0,0,1), V3(0,0,1))
    @test center(c) == P3(0.0,0.0,0.5)
    @test centroid(c) == P3(0.0,0.0,0.5)
    @test axis(c) == Line(P3(0,0,0), P3(0,0,1))
    @test isright(c)
    @test boundary(c) == CylinderSurface(T(1), Segment(P3(0,0,0), P3(0,0,1)))
    @test measure(c) == volume(c) ≈ pi 
    @test P3(0,0,0) ∈ c
    @test P3(0,0,1) ∈ c
    @test P3(1,0,0) ∈ c
    @test P3(0,1,0) ∈ c
    @test P3(cosd(60),sind(60),0.5) ∈ c
    @test P3(0,0,-0.001) ∉ c
    @test P3(0,0,1.001) ∉ c
    @test P3(1,1,1) ∉ c
  end

  @testset "CylinderSurface" begin
    c = CylinderSurface(T(2))
    @test embeddim(c) == 3
    @test paramdim(c) == 2
    @test coordtype(c) == T
    @test radius(c) == T(2)
    @test bottom(c) == Plane(P3(0,0,0), V3(0,0,1))
    @test top(c) == Plane(P3(0,0,1), V3(0,0,1))
    @test center(c) == P3(0.0,0.0,0.5)
    @test centroid(c) == P3(0.0,0.0,0.5)
    @test axis(c) == Line(P3(0,0,0), P3(0,0,1))
    @test isconvex(c)
    @test isright(c)
    @test isnothing(boundary(c))
    @test measure(c) == area(c) ≈ 2 * T(2)^2 * pi + 2 * T(2) * pi 

    c = CylinderSurface(T(5),
                 Plane(P3(1,2,3), V3(0,0,1)),
                 Plane(P3(4,5,6), V3(0,0,1)))
    @test measure(c) == area(c) ≈ 2 * T(5)^2 * pi + 2 * T(5) * pi * sqrt(3*T(3)^2)

    c = CylinderSurface(1.0)
    @test coordtype(c) == Float64
    c = CylinderSurface(1f0)
    @test coordtype(c) == Float32
    c = CylinderSurface(1)
    @test coordtype(c) == Float64
  end

  @testset "Cone" begin
    p = Plane(P3(0,0,0), V3(0,0,1))
    d = Disk(p, T(2))
    a = P3(0,0,1)
    c = Cone(d, a)
    @test embeddim(c) == 3
    @test paramdim(c) == 3
    @test coordtype(c) == T
    @test isconvex(c)
    @test boundary(c) == ConeSurface(d, a)
  end

  @testset "ConeSurface" begin
    p = Plane(P3(0,0,0), V3(0,0,1))
    d = Disk(p, T(2))
    a = P3(0,0,1)
    s = ConeSurface(d, a)
    @test embeddim(s) == 3
    @test paramdim(s) == 2
    @test coordtype(s) == T
    @test !isconvex(s)
    @test isnothing(boundary(s))
  end

  @testset "Torus" begin
    t = Torus(T.((1,1,1)), T.((1,0,0)), 2, 1)
    @test P3(1,1,-1) ∈ t
    @test P3(1,1,1) ∉ t
    @test paramdim(t) == 2
    @test !isconvex(t)
    @test isperiodic(t) == (true, true)
    @test Meshes.center(t) == P3(1,1,1)
    @test normal(t) == V3(1,0,0)
    @test radii(t) == (T(2), T(1))
    @test axis(t) == Line(P3(1,1,1), P3(2,1,1))

    # torus passing through three points
    p₁ = P3(0,0,0)
    p₂ = P3(1,2,3)
    p₃ = P3(3,2,1)
    t = Torus(p₁, p₂, p₃, 1)
    c = center(t)
    R, r = radii(t)
    @test r == 1
    @test norm(p₁ - c) ≈ R
    @test norm(p₂ - c) ≈ R
    @test norm(p₃ - c) ≈ R
    @test p₁ ∈ t
    @test p₂ ∈ t
    @test p₃ ∈ t
    
    # constructor with tuples
    c₁ = T.((0,0,0))
    c₂ = T.((1,2,3))
    c₃ = T.((3,2,1))
    q = Torus(c₁, c₂, c₃, 1)
    @test q == t
  end
end
