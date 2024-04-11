@testset "Primitives" begin
  @testset "Point" begin
    @test embeddim(Point(1)) == 1
    @test embeddim(Point(1, 2)) == 2
    @test embeddim(Point(1, 2, 3)) == 3
    @test coordtype(Point(1, 1)) == Float64
    @test coordtype(Point(1.0, 1.0)) == Float64
    @test coordtype(Point(1.0f0, 1.0f0)) == Float32
    @test coordtype(Point1(1)) == Float64
    @test coordtype(Point2(1, 1)) == Float64
    @test coordtype(Point3(1, 1, 1)) == Float64
    @test coordtype(Point1f(1)) == Float32
    @test coordtype(Point2f(1, 1)) == Float32
    @test coordtype(Point3f(1, 1, 1)) == Float32

    @test coordtype(Point{2,T}((1, 1))) == T
    @test coordtype(Point{2,T}(1, 1)) == T

    @test coordinates(P1(1)) == T[1]
    @test coordinates(P2(1, 2)) == T[1, 2]
    @test coordinates(P3(1, 2, 3)) == T[1, 2, 3]

    @test P1(1) - P1(1) == T[0]
    @test P2(1, 2) - P2(1, 1) == T[0, 1]
    @test P3(1, 2, 3) - P3(1, 1, 1) == T[0, 1, 2]
    @test_throws DimensionMismatch P2(1, 2) - P3(1, 2, 3)

    @test P1(1) + V1(0) == P1(1)
    @test P1(2) + V1(2) == P1(4)
    @test P2(1, 2) + V2(0, 0) == P2(1, 2)
    @test P2(2, 3) + V2(2, 1) == P2(4, 4)
    @test P3(1, 2, 3) + V3(0, 0, 0) == P3(1, 2, 3)
    @test P3(2, 3, 4) + V3(2, 1, 0) == P3(4, 4, 4)
    @test_throws DimensionMismatch P2(1, 2) + V3(1, 2, 3)

    @test P1(1) - V1(0) == P1(1)
    @test P1(2) - V1(2) == P1(0)
    @test P2(1, 2) - V2(0, 0) == P2(1, 2)
    @test P2(2, 3) - V2(2, 1) == P2(0, 2)
    @test P3(1, 2, 3) - V3(0, 0, 0) == P3(1, 2, 3)
    @test P3(2, 3, 4) - V3(2, 1, 0) == P3(0, 2, 4)

    @test embeddim(rand(P1)) == 1
    @test embeddim(rand(P2)) == 2
    @test embeddim(rand(P3)) == 3
    @test coordtype(rand(P1)) == T
    @test coordtype(rand(P2)) == T
    @test coordtype(rand(P3)) == T

    @test eltype(rand(P1, 3)) == P1
    @test eltype(rand(P2, 3)) == P2
    @test eltype(rand(P3, 3)) == P3

    @test P1(1) ≈ P1(1 + eps(T))
    @test P2(1, 2) ≈ P2(1 + eps(T), T(2))
    @test P3(1, 2, 3) ≈ P3(1 + eps(T), T(2), T(3))

    @test embeddim(Point((1,))) == 1
    @test coordtype(Point((1,))) == Float64
    @test coordtype(Point((1.0,))) == Float64

    @test embeddim(Point((1, 2))) == 2
    @test coordtype(Point((1, 2))) == Float64
    @test coordtype(Point((1.0, 2.0))) == Float64

    @test embeddim(Point((1, 2, 3))) == 3
    @test coordtype(Point((1, 2, 3))) == Float64
    @test coordtype(Point((1.0, 2.0, 3.0))) == Float64

    # check all 1D Point constructors, because those tend to make trouble
    @test Point(1) == Point((1,))
    @test Point{1,T}(-2) == Point{1,T}((-2,))
    @test Point{1,T}(0) == Point{1,T}((0,))

    @test_throws DimensionMismatch Point{2,T}(1)
    @test_throws DimensionMismatch Point{3,T}((2, 3))
    @test_throws DimensionMismatch Point{-3,T}((4, 5, 6))

    # There are 2 cases that throw a MethodError instead of a DimensionMismatch:
    # `Point{1,T}((2,3))` because it tries to take the tuple as a whole and convert to T and:
    # `Point{1,T}(2,3)` which does about the same.
    # I don't think this can reasonably be fixed here without hurting performance

    # check that input of mixed coordinate types is allowed and works as expected
    @test Point(1, 0.2) == Point{2,Float64}(1.0, 0.2)
    @test Point((3.0, 4)) == Point{2,Float64}(3.0, 4.0)
    @test Point((5.0, 6.0, 7)) == Point{3,Float64}(5.0, 6.0, 7.0)
    @test Point{2,T}(8, 9.0) == Point{2,T}((8.0, 9.0))
    @test Point{2,T}((-1.0, -2)) == Point{2,T}((-1, -2))
    @test Point{4,T}((0, -1.0, +2, -4.0)) == Point{4,T}((0.0f0, -1.0f0, +2.0f0, -4.0f0))

    # Integer coordinates converted to Float64
    @test coordtype(Point(1)) == Float64
    @test coordtype(Point(1, 2)) == Float64
    @test coordtype(Point(1, 2, 3)) == Float64

    # Unitful coordinates
    point = Point(1u"m", 1u"m")
    @test unit(coordtype(point)) == u"m"
    @test Unitful.numtype(coordtype(point)) === Float64
    point = Point(1.0u"m", 1.0u"m")
    @test unit(coordtype(point)) == u"m"
    @test Unitful.numtype(coordtype(point)) === Float64
    point = Point(1.0f0u"m", 1.0f0u"m")
    @test unit(coordtype(point)) == u"m"
    @test Unitful.numtype(coordtype(point)) === Float32

    # generalized inequality
    @test P2(1, 1) ⪯ P2(1, 1)
    @test !(P2(1, 1) ≺ P2(1, 1))
    @test P2(1, 2) ⪯ P2(3, 4)
    @test P2(1, 2) ≺ P2(3, 4)
    @test P2(1, 1) ⪰ P2(1, 1)
    @test !(P2(1, 1) ≻ P2(1, 1))
    @test P2(3, 4) ⪰ P2(1, 2)
    @test P2(3, 4) ≻ P2(1, 2)

    # center and centroid
    @test Meshes.center(P2(1, 1)) == P2(1, 1)
    @test centroid(P2(1, 1)) == P2(1, 1)

    # measure of points is zero
    @test measure(P2(1, 2)) == zero(T)
    @test measure(P3(1, 2, 3)) == zero(T)

    # boundary of points is nothing
    @test isnothing(boundary(rand(P1)))
    @test isnothing(boundary(rand(P2)))
    @test isnothing(boundary(rand(P3)))

    # check broadcasting works as expected
    @test P2(2, 2) .- [P2(2, 3), P2(3, 1)] == [[0.0, -1.0], [-1.0, 1.0]]
    @test P3(2, 2, 2) .- [P3(2, 3, 1), P3(3, 1, 4)] == [[0.0, -1.0, 1.0], [-1.0, 1.0, -2.0]]

    # angles between 2D points
    @test ∠(P2(0, 1), P2(0, 0), P2(1, 0)) ≈ T(-π / 2)
    @test ∠(P2(1, 0), P2(0, 0), P2(0, 1)) ≈ T(π / 2)
    @test ∠(P2(-1, 0), P2(0, 0), P2(0, 1)) ≈ T(-π / 2)
    @test ∠(P2(0, 1), P2(0, 0), P2(-1, 0)) ≈ T(π / 2)
    @test ∠(P2(0, -1), P2(0, 0), P2(1, 0)) ≈ T(π / 2)
    @test ∠(P2(1, 0), P2(0, 0), P2(0, -1)) ≈ T(-π / 2)
    @test ∠(P2(0, -1), P2(0, 0), P2(-1, 0)) ≈ T(-π / 2)
    @test ∠(P2(-1, 0), P2(0, 0), P2(0, -1)) ≈ T(π / 2)

    # angles between 3D points
    @test ∠(P3(1, 0, 0), P3(0, 0, 0), P3(0, 1, 0)) ≈ T(π / 2)
    @test ∠(P3(1, 0, 0), P3(0, 0, 0), P3(0, 0, 1)) ≈ T(π / 2)
    @test ∠(P3(0, 1, 0), P3(0, 0, 0), P3(1, 0, 0)) ≈ T(π / 2)
    @test ∠(P3(0, 1, 0), P3(0, 0, 0), P3(0, 0, 1)) ≈ T(π / 2)
    @test ∠(P3(0, 0, 1), P3(0, 0, 0), P3(1, 0, 0)) ≈ T(π / 2)
    @test ∠(P3(0, 0, 1), P3(0, 0, 0), P3(0, 1, 0)) ≈ T(π / 2)

    # a point pertains to itself
    p = P2(0, 0)
    q = P2(1, 1)
    @test p ∈ p
    @test q ∈ q
    @test p ∉ q
    @test q ∉ p
    p = P3(0, 0, 0)
    q = P3(1, 1, 1)
    @test p ∈ p
    @test q ∈ q
    @test p ∉ q
    @test q ∉ p

    p = P2(0, 1)
    @test sprint(show, p, context=:compact => true) == "(0.0, 1.0)"
    if T === Float32
      @test sprint(show, p) == "Point(0.0f0, 1.0f0)"
    else
      @test sprint(show, p) == "Point(0.0, 1.0)"
    end
  end

  @testset "Ray" begin
    r = Ray(P2(0, 0), V2(1, 1))
    @test paramdim(r) == 1
    @test measure(r) == T(Inf)
    @test length(r) == T(Inf)
    @test boundary(r) == P2(0, 0)
    @test perimeter(r) == zero(T)

    r = Ray(P2(0, 0), V2(1, 1))
    @test r(T(0.0)) == P2(0, 0)
    @test r(T(1.0)) == P2(1, 1)
    @test r(T(Inf)) == P2(Inf, Inf)
    @test r(T(1.0)) - r(T(0.0)) == V2(1, 1)
    @test_throws DomainError(T(-1), "r(t) is not defined for t < 0.") r(T(-1))

    p₁ = P3(3, 3, 3)
    p₂ = P3(-3, -3, -3)
    p₃ = P3(1, 0, 0)
    r = Ray(P3(0, 0, 0), V3(1, 1, 1))
    @test p₁ ∈ r
    @test p₂ ∉ r
    @test p₃ ∉ r

    r1 = Ray(P3(0, 0, 0), V3(1, 0, 0))
    r2 = Ray(P3(1, 1, 1), V3(1, 2, 1))
    @test r1 != r2

    r1 = Ray(P3(0, 0, 0), V3(1, 0, 0))
    r2 = Ray(P3(1, 0, 0), V3(-1, 0, 0))
    @test r1 != r2

    r1 = Ray(P3(0, 0, 0), V3(1, 0, 0))
    r2 = Ray(P3(1, 0, 0), V3(1, 0, 0))
    @test r1 != r2

    r1 = Ray(P3(0, 0, 0), V3(2, 0, 0))
    r2 = Ray(P3(0, 0, 0), V3(1, 0, 0))
    @test r1 == r2

    r2 = rand(Ray{2,T})
    r3 = rand(Ray{3,T})
    @test r2 isa Ray
    @test r3 isa Ray
    @test embeddim(r2) == 2
    @test embeddim(r3) == 3

    r = Ray(P2(0, 0), V2(1, 1))
    @test sprint(show, r) == "Ray(p: (0.0, 0.0), v: (1.0, 1.0))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), r) == """
      Ray{2,Float32}
      ├─ p: Point(0.0f0, 0.0f0)
      └─ v: Vec(1.0f0, 1.0f0)"""
    else
      @test sprint(show, MIME("text/plain"), r) == """
      Ray{2,Float64}
      ├─ p: Point(0.0, 0.0)
      └─ v: Vec(1.0, 1.0)"""
    end
  end

  @testset "Line" begin
    l = Line(P2(0, 0), P2(1, 1))
    @test paramdim(l) == 1
    @test measure(l) == T(Inf)
    @test length(l) == T(Inf)
    @test isnothing(boundary(l))
    @test perimeter(l) == zero(T)

    l = Line(P2(0, 0), P2(1, 1))
    @test (l(0), l(1)) == (P2(0, 0), P2(1, 1))

    l2 = rand(Line{2,T})
    l3 = rand(Line{3,T})
    @test l2 isa Line
    @test l3 isa Line
    @test embeddim(l2) == 2
    @test embeddim(l3) == 3

    l = Line(P2(0, 0), P2(1, 1))
    @test sprint(show, l) == "Line(a: (0.0, 0.0), b: (1.0, 1.0))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), l) == """
      Line{2,Float32}
      ├─ a: Point(0.0f0, 0.0f0)
      └─ b: Point(1.0f0, 1.0f0)"""
    else
      @test sprint(show, MIME("text/plain"), l) == """
      Line{2,Float64}
      ├─ a: Point(0.0, 0.0)
      └─ b: Point(1.0, 1.0)"""
    end
  end

  @testset "Plane" begin
    p = Plane(P3(0, 0, 0), V3(1, 0, 0), V3(0, 1, 0))
    @test p(T(1), T(0)) == P3(1, 0, 0)
    @test paramdim(p) == 2
    @test embeddim(p) == 3
    @test measure(p) == T(Inf)
    @test area(p) == T(Inf)
    @test p(T(0), T(0)) == P3(0, 0, 0)
    @test normal(p) == Vec(0, 0, 1)
    @test isnothing(boundary(p))
    @test perimeter(p) == zero(T)

    p = Plane(P3(0, 0, 0), V3(0, 0, 1))
    @test p(T(1), T(0)) == P3(1, 0, 0)
    @test p(T(0), T(1)) == P3(0, 1, 0)

    p₁ = Plane(P3(0, 0, 0), V3(1, 0, 0), V3(0, 1, 0))
    p₂ = Plane(P3(0, 0, 0), V3(0, 1, 0), V3(1, 0, 0))
    @test p₁ ≈ p₂
    p₁ = Plane(P3(0, 0, 0), V3(1, 1, 0))
    p₂ = Plane(P3(0, 0, 0), -V3(1, 1, 0))
    @test p₁ ≈ p₂

    # https://github.com/JuliaGeometry/Meshes.jl/issues/624
    p₁ = Plane(P3(0, 0, 0), V3(0, 0, 1))
    p₂ = Plane(P3(0, 0, 10), V3(0, 0, 1))
    @test !(p₁ ≈ p₂)

    # normal to plane has norm one regardless of basis
    p = Plane(P3(0, 0, 0), V3(2, 0, 0), V3(0, 3, 0))
    n = normal(p)
    @test isapprox(norm(n), T(1), atol=atol(T))

    # plane passing through three points
    p₁ = P3(0, 0, 0)
    p₂ = P3(1, 2, 3)
    p₃ = P3(3, 2, 1)
    p = Plane(p₁, p₂, p₃)
    @test p₁ ∈ p
    @test p₂ ∈ p
    @test p₃ ∈ p

    p = rand(Plane{T})
    @test p isa Plane
    @test embeddim(p) == 3

    p = Plane(P3(0, 0, 0), V3(1, 0, 0), V3(0, 1, 0))
    @test sprint(show, p) == "Plane(p: (0.0, 0.0, 0.0), u: (1.0, 0.0, 0.0), v: (0.0, 1.0, 0.0))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), p) == """
      Plane{3,Float32}
      ├─ p: Point(0.0f0, 0.0f0, 0.0f0)
      ├─ u: Vec(1.0f0, 0.0f0, 0.0f0)
      └─ v: Vec(0.0f0, 1.0f0, 0.0f0)"""
    else
      @test sprint(show, MIME("text/plain"), p) == """
      Plane{3,Float64}
      ├─ p: Point(0.0, 0.0, 0.0)
      ├─ u: Vec(1.0, 0.0, 0.0)
      └─ v: Vec(0.0, 1.0, 0.0)"""
    end
  end

  @testset "BezierCurve" begin
    b = BezierCurve(P2(0, 0), P2(0.5, 1), P2(1, 0))
    @test embeddim(b) == 2
    @test paramdim(b) == 1

    b = BezierCurve(P2(0, 0), P2(0.5, 1), P2(1, 0))
    for method in [DeCasteljau(), Horner()]
      @test b(T(0), method) == P2(0, 0)
      @test b(T(1), method) == P2(1, 0)
      @test b(T(0.5), method) == P2(0.5, 0.5)
      @test b(T(0.5), method) == P2(0.5, 0.5)
      @test_throws DomainError(T(-0.1), "b(t) is not defined for t outside [0, 1].") b(T(-0.1), method)
      @test_throws DomainError(T(1.2), "b(t) is not defined for t outside [0, 1].") b(T(1.2), method)
    end

    @test boundary(b) == Multi([P2(0, 0), P2(1, 0)])
    b = BezierCurve(P2(0, 0), P2(1, 1))
    @test boundary(b) == Multi([P2(0, 0), P2(1, 1)])
    @test perimeter(b) == zero(T)

    b = BezierCurve(P2.(randn(100), randn(100)))
    t1 = @timed b(T(0.2))
    t2 = @timed b(T(0.2), Horner())
    @test t1.time > t2.time
    @test t2.bytes < 100

    b2 = rand(BezierCurve{2,T})
    b3 = rand(BezierCurve{3,T})
    @test b2 isa BezierCurve
    @test b3 isa BezierCurve
    @test embeddim(b2) == 2
    @test embeddim(b3) == 3

    b = BezierCurve(P2(0, 0), P2(0.5, 1), P2(1, 0))
    if T === Float32
      @test sprint(show, b) == "BezierCurve(controls: Point2f[(0.0, 0.0), (0.5, 1.0), (1.0, 0.0)])"
      @test sprint(show, MIME("text/plain"), b) == """
      BezierCurve{2,Float32}
      └─ controls: Point2f[Point(0.0f0, 0.0f0), Point(0.5f0, 1.0f0), Point(1.0f0, 0.0f0)]"""
    else
      @test sprint(show, b) == "BezierCurve(controls: Point2[(0.0, 0.0), (0.5, 1.0), (1.0, 0.0)])"
      @test sprint(show, MIME("text/plain"), b) == """
      BezierCurve{2,Float64}
      └─ controls: Point2[Point(0.0, 0.0), Point(0.5, 1.0), Point(1.0, 0.0)]"""
    end
  end

  @testset "Box" begin
    b = Box(P1(0), P1(1))
    @test embeddim(b) == 1
    @test paramdim(b) == 1
    @test coordtype(b) == T
    @test minimum(b) == P1(0)
    @test maximum(b) == P1(1)
    @test extrema(b) == (P1(0), P1(1))

    b = Box(P2(0, 0), P2(1, 1))
    @test embeddim(b) == 2
    @test paramdim(b) == 2
    @test coordtype(b) == T
    @test minimum(b) == P2(0, 0)
    @test maximum(b) == P2(1, 1)
    @test extrema(b) == (P2(0, 0), P2(1, 1))

    b = Box(P3(0, 0, 0), P3(1, 1, 1))
    @test embeddim(b) == 3
    @test paramdim(b) == 3
    @test coordtype(b) == T
    @test minimum(b) == P3(0, 0, 0)
    @test maximum(b) == P3(1, 1, 1)
    @test extrema(b) == (P3(0, 0, 0), P3(1, 1, 1))

    b = Box(P1(0), P1(1))
    @test boundary(b) == Multi([P1(0), P1(1)])
    @test measure(b) == T(1)
    @test P1(0) ∈ b
    @test P1(1) ∈ b
    @test P1(0.5) ∈ b
    @test P1(-0.5) ∉ b
    @test P1(1.5) ∉ b

    b = Box(P2(0, 0), P2(1, 1))
    @test measure(b) == area(b) == T(1)
    @test P2(1, 1) ∈ b
    @test perimeter(b) ≈ T(4)

    b = Box(P2(1, 1), P2(2, 2))
    @test sides(b) == T.((1, 1))
    @test Meshes.center(b) == P2(1.5, 1.5)
    @test diagonal(b) == √T(2)

    b = Box(P2(1, 2), P2(3, 4))
    v = P2[(1, 2), (3, 2), (3, 4), (1, 4)]
    @test boundary(b) == Ring(v)

    b = Box(P3(1, 2, 3), P3(4, 5, 6))
    v = P3[(1, 2, 3), (4, 2, 3), (4, 5, 3), (1, 5, 3), (1, 2, 6), (4, 2, 6), (4, 5, 6), (1, 5, 6)]
    c = connect.([(4, 3, 2, 1), (6, 5, 1, 2), (3, 7, 6, 2), (4, 8, 7, 3), (1, 5, 8, 4), (6, 7, 8, 5)])
    @test boundary(b) == SimpleMesh(v, c)

    b = Box(P2(0, 0), P2(1, 1))
    @test boundary(b) == Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])

    b = Box(P3(0, 0, 0), P3(1, 1, 1))
    m = boundary(b)
    @test m isa Mesh
    @test nvertices(m) == 8
    @test nelements(m) == 6

    # subsetting with boxes
    b1 = Box(P2(0, 0), P2(0.5, 0.5))
    b2 = Box(P2(0.1, 0.1), P2(0.5, 0.5))
    b3 = Box(P2(0, 0), P2(1, 1))
    @test b1 ⊆ b3
    @test b2 ⊆ b3
    @test !(b1 ⊆ b2)
    @test !(b3 ⊆ b1)
    @test !(b3 ⊆ b1)

    b = Box(P2(0, 0), P2(10, 20))
    @test b(T(0.0), T(0.0)) == P2(0, 0)
    @test b(T(0.5), T(0.0)) == P2(5, 0)
    @test b(T(1.0), T(0.0)) == P2(10, 0)
    @test b(T(0.0), T(0.5)) == P2(0, 10)
    @test b(T(0.0), T(1.0)) == P2(0, 20)

    b = Box(P3(0, 0, 0), P3(10, 20, 30))
    @test b(T(0.0), T(0.0), T(0.0)) == P3(0, 0, 0)
    @test b(T(1.0), T(1.0), T(1.0)) == P3(10, 20, 30)

    b1 = rand(Box{1,T})
    b2 = rand(Box{2,T})
    b3 = rand(Box{3,T})
    @test b1 isa Box
    @test b2 isa Box
    @test b3 isa Box
    @test embeddim(b1) == 1
    @test embeddim(b2) == 2
    @test embeddim(b3) == 3

    @test_throws AssertionError Box(P1(1), P1(0))
    @test_throws AssertionError Box(P2(1, 1), P2(0, 0))
    @test_throws AssertionError Box(P3(1, 1, 1), P3(0, 0, 0))

    b = Box(P2(0, 0), P2(1, 1))
    q = convert(Quadrangle, b)
    @test q isa Quadrangle
    @test q == Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))

    b = Box(P3(0, 0, 0), P3(1, 1, 1))
    h = convert(Hexahedron, b)
    @test h isa Hexahedron
    @test h == Hexahedron(
      P3(0, 0, 0),
      P3(1, 0, 0),
      P3(1, 1, 0),
      P3(0, 1, 0),
      P3(0, 0, 1),
      P3(1, 0, 1),
      P3(1, 1, 1),
      P3(0, 1, 1)
    )

    b = Box(P2(0, 0), P2(1, 1))
    @test sprint(show, b) == "Box(min: (0.0, 0.0), max: (1.0, 1.0))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), b) == """
      Box{2,Float32}
      ├─ min: Point(0.0f0, 0.0f0)
      └─ max: Point(1.0f0, 1.0f0)"""
    else
      @test sprint(show, MIME("text/plain"), b) == """
      Box{2,Float64}
      ├─ min: Point(0.0, 0.0)
      └─ max: Point(1.0, 1.0)"""
    end
  end

  @testset "Ball" begin
    b = Ball(P3(1, 2, 3), T(5))
    @test embeddim(b) == 3
    @test paramdim(b) == 3
    @test coordtype(b) == T
    @test Meshes.center(b) == P3(1, 2, 3)
    @test radius(b) == T(5)

    b = Ball(P3(1, 2, 3), 4)
    @test coordtype(b) == T

    b1 = Ball(P2(0, 0), T(1))
    b2 = Ball(P2(0, 0))
    b3 = Ball(T.((0, 0)))
    @test b1 == b2 == b3

    b = Ball(P2(0, 0), T(2))
    @test measure(b) ≈ T(π) * (T(2)^2)
    b = Ball(P3(0, 0, 0), T(2))
    @test measure(b) ≈ T(4 / 3) * T(π) * (T(2)^3)
    @test_throws ArgumentError length(b)
    @test_throws ArgumentError area(b)

    b = Ball(P2(0, 0), T(2))
    @test P2(1, 0) ∈ b
    @test P2(0, 1) ∈ b
    @test P2(2, 0) ∈ b
    @test P2(0, 2) ∈ b
    @test P2(3, 5) ∉ b
    @test perimeter(b) ≈ T(4π)

    b = Ball(P3(0, 0, 0), T(2))
    @test P3(1, 0, 0) ∈ b
    @test P3(0, 0, 1) ∈ b
    @test P3(2, 0, 0) ∈ b
    @test P3(0, 0, 2) ∈ b
    @test P3(3, 5, 2) ∉ b

    b = Ball(P2(0, 0), T(2))
    @test b(T(0), T(0)) ≈ P2(0, 0)
    @test b(T(1), T(0)) ≈ P2(2, 0)

    b = Ball(P2(7, 7), T(1.5))
    ps = b.(1, rand(T, 100))
    all(∈(b), ps)

    b = Ball(P3(0, 0, 0), T(2))
    @test b(T(0), T(0), T(0)) ≈ P3(0, 0, 0)
    @test b(T(1), T(0), T(0)) ≈ P3(0, 0, 2)

    b = Ball(P3(7, 7, 7), T(1.5))
    ps = b.(1, rand(T, 100), rand(T, 100))
    all(∈(b), ps)

    b1 = rand(Ball{1,T})
    b2 = rand(Ball{2,T})
    b3 = rand(Ball{3,T})
    @test b1 isa Ball
    @test b2 isa Ball
    @test b3 isa Ball
    @test embeddim(b1) == 1
    @test embeddim(b2) == 2
    @test embeddim(b3) == 3

    b = Ball(P2(0, 0), T(1))
    @test sprint(show, b) == "Ball(center: (0.0, 0.0), radius: 1.0)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), b) == """
      Ball{2,Float32}
      ├─ center: Point(0.0f0, 0.0f0)
      └─ radius: 1.0f0"""
    else
      @test sprint(show, MIME("text/plain"), b) == """
      Ball{2,Float64}
      ├─ center: Point(0.0, 0.0)
      └─ radius: 1.0"""
    end
  end

  @testset "Sphere" begin
    s = Sphere(P3(0, 0, 0), T(1))
    @test embeddim(s) == 3
    @test paramdim(s) == 2
    @test coordtype(s) == T
    @test Meshes.center(s) == P3(0, 0, 0)
    @test radius(s) == T(1)
    @test extrema(s) == (P3(-1, -1, -1), P3(1, 1, 1))
    @test isnothing(boundary(s))
    @test perimeter(s) == zero(T)

    s = Sphere(P3(1, 2, 3), 4)
    @test coordtype(s) == T

    s = Sphere(P2(0, 0), T(1))
    @test embeddim(s) == 2
    @test paramdim(s) == 1
    @test coordtype(s) == T
    @test Meshes.center(s) == P2(0, 0)
    @test radius(s) == T(1)
    @test extrema(s) == (P2(-1, -1), P2(1, 1))
    @test isnothing(boundary(s))

    s1 = Sphere(P2(0, 0), T(1))
    s2 = Sphere(P2(0, 0))
    s3 = Sphere(T.((0, 0)))
    @test s1 == s2 == s3

    s = Sphere(P2(0, 0), T(2))
    @test measure(s) ≈ 2π * 2
    @test length(s) ≈ 2π * 2
    @test extrema(s) == (P2(-2, -2), P2(2, 2))
    s = Sphere(P3(0, 0, 0), T(2))
    @test measure(s) ≈ 4π * (2^2)
    @test area(s) ≈ 4π * (2^2)

    s = Sphere(P2(0, 0), T(2))
    @test P2(1, 0) ∉ s
    @test P2(0, 1) ∉ s
    @test P2(2, 0) ∈ s
    @test P2(0, 2) ∈ s
    @test P2(3, 5) ∉ s

    s = Sphere(P3(0, 0, 0), T(2))
    @test P3(1, 0, 0) ∉ s
    @test P3(0, 0, 1) ∉ s
    @test P3(2, 0, 0) ∈ s
    @test P3(0, 0, 2) ∈ s
    @test P3(3, 5, 2) ∉ s

    # 2D sphere passing through 3 points
    s = Sphere(P2(0, 0), P2(0.5, 0), P2(1, 1))
    @test Meshes.center(s) == P2(0.25, 0.75)
    @test radius(s) == T(0.7905694150420949)
    s = Sphere(P2(0, 0), P2(1, 0), P2(0, 1))
    @test Meshes.center(s) == P2(0.5, 0.5)
    @test radius(s) == T(0.7071067811865476)
    s = Sphere(P2(0, 0), P2(1, 0), P2(1, 1))
    @test Meshes.center(s) == P2(0.5, 0.5)
    @test radius(s) == T(0.7071067811865476)

    # 3D sphere passing through 4 points
    s = Sphere(P3(0, 0, 0), P3(5, 0, 1), P3(1, 1, 1), P3(3, 2, 1))
    @test P3(0, 0, 0) ∈ s
    @test P3(5, 0, 1) ∈ s
    @test P3(1, 1, 1) ∈ s
    @test P3(3, 2, 1) ∈ s
    O = Meshes.center(s)
    r = radius(s)
    @test isapprox(r, norm(P3(0, 0, 0) - O))

    s = Sphere(P2(0, 0), T(2))
    @test s(T(0)) ≈ P2(2, 0)
    @test s(T(0.5)) ≈ P2(-2, 0)

    s = Sphere(P3(0, 0, 0), T(2))
    @test s(T(0), T(0)) ≈ P3(0, 0, 2)
    @test s(T(0.5), T(0.5)) ≈ P3(-2, 0, 0)

    s1 = rand(Sphere{1,T})
    s2 = rand(Sphere{2,T})
    s3 = rand(Sphere{3,T})
    @test s1 isa Sphere
    @test s2 isa Sphere
    @test s3 isa Sphere
    @test embeddim(s1) == 1
    @test embeddim(s2) == 2
    @test embeddim(s3) == 3

    s = Sphere(P3(0, 0, 0), T(1))
    @test sprint(show, s) == "Sphere(center: (0.0, 0.0, 0.0), radius: 1.0)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), s) == """
      Sphere{3,Float32}
      ├─ center: Point(0.0f0, 0.0f0, 0.0f0)
      └─ radius: 1.0f0"""
    else
      @test sprint(show, MIME("text/plain"), s) == """
      Sphere{3,Float64}
      ├─ center: Point(0.0, 0.0, 0.0)
      └─ radius: 1.0"""
    end
  end

  @testset "Ellipsoid" begin
    e = Ellipsoid((T(3), T(2), T(1)))
    @test embeddim(e) == 3
    @test paramdim(e) == 2
    @test coordtype(e) == T
    @test radii(e) == (T(3), T(2), T(1))
    @test center(e) == P3(0, 0, 0)
    @test isnothing(boundary(e))
    @test perimeter(e) == zero(T)

    e = Ellipsoid((T(3), T(2), T(1)))
    @test sprint(show, e) == "Ellipsoid(radii: (3.0, 2.0, 1.0), center: (0.0, 0.0, 0.0), rotation: UniformScaling{Bool}(true))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), e) == """
      Ellipsoid{3,Float32}
      ├─ radii: (3.0f0, 2.0f0, 1.0f0)
      ├─ center: Point(0.0f0, 0.0f0, 0.0f0)
      └─ rotation: UniformScaling{Bool}(true)"""
    else
      @test sprint(show, MIME("text/plain"), e) == """
      Ellipsoid{3,Float64}
      ├─ radii: (3.0, 2.0, 1.0)
      ├─ center: Point(0.0, 0.0, 0.0)
      └─ rotation: UniformScaling{Bool}(true)"""
    end
  end

  @testset "Disk" begin
    p = Plane(P3(0, 0, 0), V3(0, 0, 1))
    d = Disk(p, T(2))
    @test embeddim(d) == 3
    @test paramdim(d) == 2
    @test coordtype(d) == T
    @test plane(d) == p
    @test Meshes.center(d) == P3(0, 0, 0)
    @test radius(d) == T(2)
    @test normal(d) == V3(0, 0, 1)
    @test measure(d) == T(π) * T(2)^2
    @test area(d) == measure(d)
    @test P3(0, 0, 0) ∈ d
    @test P3(0, 0, 1) ∉ d
    @test boundary(d) == Circle(p, T(2))

    d = rand(Disk{T})
    @test d isa Disk
    @test embeddim(d) == 3

    p = Plane(P3(0, 0, 0), V3(0, 0, 1))
    d = Disk(p, T(2))
    @test sprint(show, d) ==
          "Disk(plane: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0)), radius: 2.0)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), d) == """
      Disk{3,Float32}
      ├─ plane: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0))
      └─ radius: 2.0f0"""
    else
      @test sprint(show, MIME("text/plain"), d) == """
      Disk{3,Float64}
      ├─ plane: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0))
      └─ radius: 2.0"""
    end
  end

  @testset "Circle" begin
    p = Plane(P3(0, 0, 0), V3(0, 0, 1))
    c = Circle(p, T(2))
    @test embeddim(c) == 3
    @test paramdim(c) == 1
    @test coordtype(c) == T
    @test plane(c) == p
    @test Meshes.center(c) == P3(0, 0, 0)
    @test radius(c) == T(2)
    @test measure(c) == 2 * T(π) * T(2)
    @test length(c) == measure(c)
    @test P3(2, 0, 0) ∈ c
    @test P3(0, 2, 0) ∈ c
    @test P3(0, 0, 0) ∉ c
    @test isnothing(boundary(c))

    # 3D circumcircle
    p1 = P3(0, 4, 0)
    p2 = P3(0, -4, 0)
    p3 = P3(0, 0, 4)
    c = Circle(p1, p2, p3)
    @test p1 ∈ c
    @test p2 ∈ c
    @test p3 ∈ c

    # circle parametrization
    p = Plane(P3(0, 0, 0), V3(0, 0, 1))
    c = Circle(p, T(2))
    @test c(T(0)) ≈ P3(2, 0, 0)
    @test c(T(0.25)) ≈ P3(0, 2, 0)
    @test c(T(0.5)) ≈ P3(-2, 0, 0)
    @test c(T(0.75)) ≈ P3(0, -2, 0)
    @test c(T(1)) ≈ P3(2, 0, 0)

    c = rand(Circle{T})
    @test c isa Circle
    @test embeddim(c) == 3

    p = Plane(P3(0, 0, 0), V3(0, 0, 1))
    c = Circle(p, T(2))
    @test sprint(show, c) ==
          "Circle(plane: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0)), radius: 2.0)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Circle{3,Float32}
      ├─ plane: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0))
      └─ radius: 2.0f0"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Circle{3,Float64}
      ├─ plane: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0))
      └─ radius: 2.0"""
    end
  end

  @testset "Cylinder" begin
    c = Cylinder(Plane(P3(1, 2, 3), V3(0, 0, 1)), Plane(P3(4, 5, 6), V3(0, 0, 1)), T(5))
    @test embeddim(c) == 3
    @test paramdim(c) == 3
    @test coordtype(c) == T
    @test radius(c) == T(5)
    @test bottom(c) == Plane(P3(1, 2, 3), V3(0, 0, 1))
    @test top(c) == Plane(P3(4, 5, 6), V3(0, 0, 1))
    @test axis(c) == Line(P3(1, 2, 3), P3(4, 5, 6))
    @test !isright(c)
    @test measure(c) == volume(c) ≈ T(5)^2 * pi * T(3) * sqrt(T(3))
    @test P3(1, 2, 3) ∈ c
    @test P3(4, 5, 6) ∈ c
    @test P3(0.99, 1.99, 2.99) ∉ c
    @test P3(4.01, 5.01, 6.01) ∉ c
    @test !Meshes.hasintersectingplanes(c)
    @test c(0, 0, 0) ≈ bottom(c)(0, 0)
    @test c(0, 0, 1) ≈ top(c)(0, 0)
    @test c(1, 0.25, 0.5) ≈ Point(T(4.330127018922193), T(10.330127018922191), T(4.5))
    @test_throws DomainError c(1.1, 0, 0)

    c = Cylinder(Plane(P3(0, 0, 0), V3(0, 0, 1)), Plane(P3(0, 0, 1), V3(1, 0, 1)), T(5))
    @test Meshes.hasintersectingplanes(c)

    c1 = Cylinder(P3(0, 0, 0), P3(0, 0, 1), T(1))
    c2 = Cylinder(P3(0, 0, 0), P3(0, 0, 1))
    c3 = Cylinder(T(1))
    @test c1 == c2 == c3
    @test c1 ≈ c2 ≈ c3

    c = Cylinder(T(1))
    @test coordtype(c) == T
    c = Cylinder(1)
    @test coordtype(c) == Float64

    c = Cylinder(P3(0, 0, 0), P3(0, 0, 1), T(1))
    @test radius(c) == T(1)
    @test bottom(c) == Plane(P3(0, 0, 0), V3(0, 0, 1))
    @test top(c) == Plane(P3(0, 0, 1), V3(0, 0, 1))
    @test center(c) == P3(0.0, 0.0, 0.5)
    @test centroid(c) == P3(0.0, 0.0, 0.5)
    @test axis(c) == Line(P3(0, 0, 0), P3(0, 0, 1))
    @test isright(c)
    @test boundary(c) == CylinderSurface(P3(0, 0, 0), P3(0, 0, 1), T(1))
    @test measure(c) == volume(c) ≈ pi
    @test P3(0, 0, 0) ∈ c
    @test P3(0, 0, 1) ∈ c
    @test P3(1, 0, 0) ∈ c
    @test P3(0, 1, 0) ∈ c
    @test P3(cosd(60), sind(60), 0.5) ∈ c
    @test P3(0, 0, -0.001) ∉ c
    @test P3(0, 0, 1.001) ∉ c
    @test P3(1, 1, 1) ∉ c

    c = rand(Cylinder{T})
    @test c isa Cylinder
    @test embeddim(c) == 3

    c = Cylinder(P3(0, 0, 0), P3(0, 0, 1), T(1))
    @test sprint(show, c) ==
          "Cylinder(bot: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0)), top: Plane(p: (0.0, 0.0, 1.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0)), radius: 1.0)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Cylinder{3,Float32}
      ├─ bot: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0))
      ├─ top: Plane(p: (0.0, 0.0, 1.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0))
      └─ radius: 1.0f0"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Cylinder{3,Float64}
      ├─ bot: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0))
      ├─ top: Plane(p: (0.0, 0.0, 1.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0))
      └─ radius: 1.0"""
    end
  end

  @testset "CylinderSurface" begin
    c = CylinderSurface(T(2))
    @test embeddim(c) == 3
    @test paramdim(c) == 2
    @test coordtype(c) == T
    @test radius(c) == T(2)
    @test bottom(c) == Plane(P3(0, 0, 0), V3(0, 0, 1))
    @test top(c) == Plane(P3(0, 0, 1), V3(0, 0, 1))
    @test center(c) == P3(0.0, 0.0, 0.5)
    @test centroid(c) == P3(0.0, 0.0, 0.5)
    @test axis(c) == Line(P3(0, 0, 0), P3(0, 0, 1))
    @test isright(c)
    @test isnothing(boundary(c))
    @test measure(c) == area(c) ≈ 2 * T(2)^2 * pi + 2 * T(2) * pi
    @test !Meshes.hasintersectingplanes(c)

    c = CylinderSurface(Plane(P3(0, 0, 0), V3(0, 0, 1)), Plane(P3(0, 0, 1), V3(1, 0, 1)), T(5))
    @test Meshes.hasintersectingplanes(c)

    c1 = CylinderSurface(P3(0, 0, 0), P3(0, 0, 1), T(1))
    c2 = CylinderSurface(P3(0, 0, 0), P3(0, 0, 1))
    c3 = CylinderSurface(T(1))
    @test c1 == c2 == c3
    @test c1 ≈ c2 ≈ c3

    c = CylinderSurface(Plane(P3(1, 2, 3), V3(0, 0, 1)), Plane(P3(4, 5, 6), V3(0, 0, 1)), T(5))
    @test measure(c) == area(c) ≈ 2 * T(5)^2 * pi + 2 * T(5) * pi * sqrt(3 * T(3)^2)

    c = CylinderSurface(T(1))
    @test c(T(0), T(0)) ≈ P3(1, 0, 0)
    @test c(T(0.5), T(0)) ≈ P3(-1, 0, 0)
    @test c(T(0), T(1)) ≈ P3(1, 0, 1)
    @test c(T(0.5), T(1)) ≈ P3(-1, 0, 1)

    c = CylinderSurface(1.0)
    @test coordtype(c) == Float64
    c = CylinderSurface(1.0f0)
    @test coordtype(c) == Float32
    c = CylinderSurface(1)
    @test coordtype(c) == Float64

    c = rand(CylinderSurface{T})
    @test c isa CylinderSurface
    @test embeddim(c) == 3

    c = CylinderSurface(T(1))
    @test sprint(show, c) ==
          "CylinderSurface(bot: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0)), top: Plane(p: (0.0, 0.0, 1.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0)), radius: 1.0)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      CylinderSurface{3,Float32}
      ├─ bot: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0))
      ├─ top: Plane(p: (0.0, 0.0, 1.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0))
      └─ radius: 1.0f0"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      CylinderSurface{3,Float64}
      ├─ bot: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0))
      ├─ top: Plane(p: (0.0, 0.0, 1.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0))
      └─ radius: 1.0"""
    end
  end

  @testset "ParaboloidSurface" begin
    p = ParaboloidSurface(P3(0, 0, 0), T(1), T(2))
    @test embeddim(p) == 3
    @test paramdim(p) == 2
    @test coordtype(p) == T
    @test focallength(p) == T(2)
    @test radius(p) == T(1)
    @test axis(p) == Line(P3(0, 0, 0), P3(0, 0, T(2)))
    @test measure(p) == area(p) ≈ T(32π / 3 * (17√17 / 64 - 1))

    p1 = ParaboloidSurface(P3(1, 2, 3), T(1), T(1))
    p2 = ParaboloidSurface(P3(1, 2, 3), T(1))
    p3 = ParaboloidSurface(P3(1, 2, 3))
    @test p1 == p2 == p3
    @test p1 ≈ p2 ≈ p3

    p1 = ParaboloidSurface((1, 2, 3), 1.0, 1.0)
    p2 = ParaboloidSurface((1, 2, 3), 1.0)
    p3 = ParaboloidSurface((1, 2, 3))
    @test p1 == p2 == p3
    @test p1 ≈ p2 ≈ p3

    p = ParaboloidSurface((1.0, 2.0, 3.0), 4.0, 5.0)
    @test coordtype(p) == Float64
    @test radius(p) == 4.0
    @test focallength(p) == 5.0

    p = ParaboloidSurface(P3(1, 5, 2), T(3), T(4))
    @test measure(p) == area(p) ≈ T(128π / 3 * (73√73 / 512 - 1))
    @test p(T(0), T(0)) ≈ P3(1, 5, 2)
    @test p(T(1), T(0)) ≈ P3(4, 5, 2 + 3^2 / (4 * 4))
    @test_throws DomainError p(T(-0.1), T(0))
    @test_throws DomainError p(T(1.1), T(0))

    p = ParaboloidSurface()
    @test coordtype(p) == Float64
    @test p(0.0, 0.0) ≈ Point3(0, 0, 0)
    @test p(0.5, 0.0) ≈ Point3(0.5, 0, 0.5^2 / 4)
    @test p(0.0, 0.5) ≈ Point3(0, 0, 0)
    @test p(0.5, 0.5) ≈ Point3(-0.5, 0, 0.5^2 / 4)

    p = ParaboloidSurface(Point3(0, 0, 0))
    @test coordtype(p) == Float64
    p = ParaboloidSurface(Point3f(0, 0, 0))
    @test coordtype(p) == Float32

    p = rand(ParaboloidSurface{T})
    @test p isa ParaboloidSurface
    @test embeddim(p) == 3

    p = ParaboloidSurface(P3(0, 0, 0))
    @test sprint(show, p) == "ParaboloidSurface(apex: (0.0, 0.0, 0.0), radius: 1.0, focallength: 1.0)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), p) == """
      ParaboloidSurface{3,Float32}
      ├─ apex: Point(0.0f0, 0.0f0, 0.0f0)
      ├─ radius: 1.0f0
      └─ focallength: 1.0f0"""
    else
      @test sprint(show, MIME("text/plain"), p) == """
      ParaboloidSurface{3,Float64}
      ├─ apex: Point(0.0, 0.0, 0.0)
      ├─ radius: 1.0
      └─ focallength: 1.0"""
    end
  end

  @testset "Cone" begin
    p = Plane(P3(0, 0, 0), V3(0, 0, 1))
    d = Disk(p, T(2))
    a = P3(0, 0, 1)
    c = Cone(d, a)
    @test embeddim(c) == 3
    @test paramdim(c) == 3
    @test coordtype(c) == T
    @test boundary(c) == ConeSurface(d, a)

    c = rand(Cone{T})
    @test c isa Cone
    @test embeddim(c) == 3

    p = Plane(P3(0, 0, 0), V3(0, 0, 1))
    d = Disk(p, T(2))
    a = P3(0, 0, 1)
    c = Cone(d, a)
    @test sprint(show, c) ==
          "Cone(base: Disk(plane: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0)), radius: 2.0), apex: (0.0, 0.0, 1.0))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Cone{3,Float32}
      ├─ base: Disk(plane: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0)), radius: 2.0)
      └─ apex: Point(0.0f0, 0.0f0, 1.0f0)"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Cone{3,Float64}
      ├─ base: Disk(plane: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0)), radius: 2.0)
      └─ apex: Point(0.0, 0.0, 1.0)"""
    end

    # cone: apex at (5,4,3); base center at (5,1,3)
    # halfangle: 30° -> radius: sqrt(3)
    # axis of the cone is parallel to y axis
    p = Plane(P3(5, 1, 3), V3(0, 1, 0))
    d = Disk(p, sqrt(T(3)))
    a = P3(5, 4, 3)
    c = Cone(d, a)

    @test rad2deg(Meshes.halfangle(c)) ≈ T(30)
    @test Meshes.height(c) ≈ T(3)

    @test P3(5, 1, 3) ∈ c
    @test P3(5, 4, 3) ∈ c
    @test P3(5, 1, 3 - sqrt(3)) ∈ c
    @test P3(5, 1, 3 + sqrt(3)) ∈ c
    @test P3(5 - sqrt(3), 1, 3) ∈ c
    @test P3(5 + sqrt(3), 1, 3) ∈ c
    @test P3(5, 2.5, 3) ∈ c
    @test P3(5 + sqrt(3) / 2, 2.5, 3) ∈ c
    @test P3(5 - sqrt(3) / 2, 2.5, 3) ∈ c

    @test P3(5, 0.9, 3) ∉ c
    @test P3(5, 4.1, 3) ∉ c
    @test P3(5, 1, 1) ∉ c
    @test P3(5 + sqrt(3) + 0.01, 1, 3) ∉ c
    @test P3(5 + sqrt(3) / 2 + 0.01, 2.5, 3) ∉ c
    @test P3(5 - sqrt(3) / 2 - 0.01, 2.5, 3) ∉ c
  end

  @testset "ConeSurface" begin
    p = Plane(P3(0, 0, 0), V3(0, 0, 1))
    d = Disk(p, T(2))
    a = P3(0, 0, 1)
    s = ConeSurface(d, a)
    @test embeddim(s) == 3
    @test paramdim(s) == 2
    @test coordtype(s) == T
    @test isnothing(boundary(s))

    c = rand(ConeSurface{T})
    @test c isa ConeSurface
    @test embeddim(c) == 3

    p = Plane(P3(0, 0, 0), V3(0, 0, 1))
    d = Disk(p, T(2))
    a = P3(0, 0, 1)
    s = ConeSurface(d, a)
    @test sprint(show, s) ==
          "ConeSurface(base: Disk(plane: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0)), radius: 2.0), apex: (0.0, 0.0, 1.0))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), s) == """
      ConeSurface{3,Float32}
      ├─ base: Disk(plane: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0)), radius: 2.0)
      └─ apex: Point(0.0f0, 0.0f0, 1.0f0)"""
    else
      @test sprint(show, MIME("text/plain"), s) == """
      ConeSurface{3,Float64}
      ├─ base: Disk(plane: Plane(p: (0.0, 0.0, 0.0), u: (1.0, -0.0, -0.0), v: (-0.0, 1.0, -0.0)), radius: 2.0)
      └─ apex: Point(0.0, 0.0, 1.0)"""
    end
  end

  @testset "Frustum" begin
    pb = Plane(P3(0, 0, 0), V3(0, 0, 1))
    db = Disk(pb, T(1))
    pt = Plane(P3(0, 0, 10), V3(0, 0, 1))
    dt = Disk(pt, T(2))
    f = Frustum(db, dt)
    @test embeddim(f) == 3
    @test coordtype(f) == T
    @test boundary(f) == FrustumSurface(db, dt)

    @test_throws AssertionError Frustum(db, db)

    f = rand(Frustum{T})
    @test f isa Frustum

    f = Frustum(db, dt)
    @test P3(0, 0, 0) ∈ f
    @test P3(0, 0, 10) ∈ f
    @test P3(1, 0, 0) ∈ f
    @test P3(2, 0, 10) ∈ f
    @test P3(1, 0, 5) ∈ f

    @test P3(1, 1, 0) ∉ f
    @test P3(2, 2, 10) ∉ f
    @test P3(0, 0, -0.01) ∉ f
    @test P3(0, 0, 10.01) ∉ f

    # reverse order, when top is larger than bottom
    # the frustum is the same geometry
    f = Frustum(dt, db)
    @test P3(0, 0, 0) ∈ f
    @test P3(0, 0, 10) ∈ f
    @test P3(1, 0, 0) ∈ f
    @test P3(2, 0, 10) ∈ f
    @test P3(1, 0, 5) ∈ f

    @test P3(1, 1, 0) ∉ f
    @test P3(2, 2, 10) ∉ f
    @test P3(0, 0, -0.01) ∉ f
    @test P3(0, 0, 10.01) ∉ f
  end

  @testset "FrustumSurface" begin
    pb = Plane(P3(0, 0, 0), V3(0, 0, 1))
    db = Disk(pb, T(1))
    pt = Plane(P3(0, 0, 10), V3(0, 0, 1))
    dt = Disk(pt, T(2))
    f = FrustumSurface(db, dt)
    @test embeddim(f) == 3
    @test paramdim(f) == 2
    @test coordtype(f) == T
    @test isnothing(boundary(f))

    @test_throws AssertionError FrustumSurface(db, db)

    f = rand(FrustumSurface{T})
    @test f isa FrustumSurface
  end

  @testset "Torus" begin
    t = Torus(T.((1, 1, 1)), T.((1, 0, 0)), 2, 1)
    @test P3(1, 1, -1) ∈ t
    @test P3(1, 1, 1) ∉ t
    @test paramdim(t) == 2
    @test Meshes.center(t) == P3(1, 1, 1)
    @test normal(t) == V3(1, 0, 0)
    @test radii(t) == (T(2), T(1))
    @test axis(t) == Line(P3(1, 1, 1), P3(2, 1, 1))
    @test measure(t) ≈ 8 * T(π)^2
    @test_throws ArgumentError length(t)
    @test_throws ArgumentError volume(t)

    # torus passing through three points
    p₁ = P3(0, 0, 0)
    p₂ = P3(1, 2, 3)
    p₃ = P3(3, 2, 1)
    t = Torus(p₁, p₂, p₃, T(1))
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
    c₁ = T.((0, 0, 0))
    c₂ = T.((1, 2, 3))
    c₃ = T.((3, 2, 1))
    q = Torus(c₁, c₂, c₃, 1)
    @test q == t

    t = rand(Torus{T})
    @test t isa Torus
    @test embeddim(t) == 3
    @test coordtype(t) == T
    @test isnothing(boundary(t))

    t = Torus(P3(1, 1, 1), V3(1, 0, 0), 2, 1)
    @test sprint(show, t) == "Torus(center: (1.0, 1.0, 1.0), normal: (1.0, 0.0, 0.0), major: 2.0, minor: 1.0)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), t) == """
      Torus{3,Float32}
      ├─ center: Point(1.0f0, 1.0f0, 1.0f0)
      ├─ normal: Vec(1.0f0, 0.0f0, 0.0f0)
      ├─ major: 2.0f0
      └─ minor: 1.0f0"""
    else
      @test sprint(show, MIME("text/plain"), t) == """
      Torus{3,Float64}
      ├─ center: Point(1.0, 1.0, 1.0)
      ├─ normal: Vec(1.0, 0.0, 0.0)
      ├─ major: 2.0
      └─ minor: 1.0"""
    end
  end
end
