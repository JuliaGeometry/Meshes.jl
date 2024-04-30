@testset "Primitives" begin
  @testset "Point" begin
    @test embeddim(Point(1)) == 1
    @test embeddim(Point(1, 2)) == 2
    @test embeddim(Point(1, 2, 3)) == 3
    @test Meshes.lentype(Point(1, 1)) == Meshes.Met{Float64}
    @test Meshes.lentype(Point(1.0, 1.0)) == Meshes.Met{Float64}
    @test Meshes.lentype(Point(1.0f0, 1.0f0)) == Meshes.Met{Float32}

    @test Meshes.lentype(Point((T(1), T(1)))) == ℳ
    @test Meshes.lentype(Point(T(1), T(1))) == ℳ

    @test coordinates(point(1)) == vec(1)
    @test coordinates(point(1, 2)) == vec(1, 2)
    @test coordinates(point(1, 2, 3)) == vec(1, 2, 3)

    @test point(1) - point(1) == vec(0)
    @test point(1, 2) - point(1, 1) == vec(0, 1)
    @test point(1, 2, 3) - point(1, 1, 1) == vec(0, 1, 2)
    @test_throws MethodError point(1, 2) - point(1, 2, 3)

    @test point(1) + vec(0) == point(1)
    @test point(2) + vec(2) == point(4)
    @test point(1, 2) + vec(0, 0) == point(1, 2)
    @test point(2, 3) + vec(2, 1) == point(4, 4)
    @test point(1, 2, 3) + vec(0, 0, 0) == point(1, 2, 3)
    @test point(2, 3, 4) + vec(2, 1, 0) == point(4, 4, 4)
    @test_throws MethodError point(1, 2) + vec(1, 2, 3)

    @test point(1) - vec(0) == point(1)
    @test point(2) - vec(2) == point(0)
    @test point(1, 2) - vec(0, 0) == point(1, 2)
    @test point(2, 3) - vec(2, 1) == point(0, 2)
    @test point(1, 2, 3) - vec(0, 0, 0) == point(1, 2, 3)
    @test point(2, 3, 4) - vec(2, 1, 0) == point(0, 2, 4)

    @test embeddim(rand(Point{1})) == 1
    @test embeddim(rand(Point{2})) == 2
    @test embeddim(rand(Point{3})) == 3
    @test Meshes.lentype(rand(Point{1})) == Meshes.Met{Float64}
    @test Meshes.lentype(rand(Point{2})) == Meshes.Met{Float64}
    @test Meshes.lentype(rand(Point{3})) == Meshes.Met{Float64}

    @test eltype(rand(Point{1}, 3)) == Point{1}
    @test eltype(rand(Point{2}, 3)) == Point{2}
    @test eltype(rand(Point{3}, 3)) == Point{3}

    @test point(1) ≈ point(1 + eps(T))
    @test point(1, 2) ≈ point(1 + eps(T), T(2))
    @test point(1, 2, 3) ≈ point(1 + eps(T), T(2), T(3))

    @test embeddim(Point((1,))) == 1
    @test Meshes.lentype(Point((1,))) == Meshes.Met{Float64}
    @test Meshes.lentype(Point((1.0,))) == Meshes.Met{Float64}

    @test embeddim(Point((1, 2))) == 2
    @test Meshes.lentype(Point((1, 2))) == Meshes.Met{Float64}
    @test Meshes.lentype(Point((1.0, 2.0))) == Meshes.Met{Float64}

    @test embeddim(Point((1, 2, 3))) == 3
    @test Meshes.lentype(Point((1, 2, 3))) == Meshes.Met{Float64}
    @test Meshes.lentype(Point((1.0, 2.0, 3.0))) == Meshes.Met{Float64}

    # check all 1D Point constructors, because those tend to make trouble
    @test Point(1) == Point((1,))
    @test Point(T(-2)) == Point((T(-2),))
    @test Point(T(0)) == Point((T(0),))

    # check that input of mixed coordinate types is allowed and works as expected
    @test Point(1, 0.2) == Point(1.0, 0.2)
    @test Point((3.0, 4)) == Point(3.0, 4.0)
    @test Point((5.0, 6.0, 7)) == Point(5.0, 6.0, 7.0)
    @test Point(8, T(9.0)) == Point((T(8.0), T(9.0)))
    @test Point((T(-1.0), -2)) == Point((T(-1.0), T(-2.0)))
    @test Point((0, T(-1.0), +2, T(-4.0))) == Point((T(0.0), T(-1.0), T(+2.0), T(-4.0)))

    # Integer coordinates converted to Float64
    @test Meshes.lentype(Point(1)) == Meshes.Met{Float64}
    @test Meshes.lentype(Point(1, 2)) == Meshes.Met{Float64}
    @test Meshes.lentype(Point(1, 2, 3)) == Meshes.Met{Float64}

    # Unitful coordinates
    p = Point(1u"m", 1u"m")
    @test unit(Meshes.lentype(p)) == u"m"
    @test Unitful.numtype(Meshes.lentype(p)) === Float64
    p = Point(1.0u"m", 1.0u"m")
    @test unit(Meshes.lentype(p)) == u"m"
    @test Unitful.numtype(Meshes.lentype(p)) === Float64
    p = Point(1.0f0u"m", 1.0f0u"m")
    @test unit(Meshes.lentype(p)) == u"m"
    @test Unitful.numtype(Meshes.lentype(p)) === Float32

    # generalized inequality
    @test point(1, 1) ⪯ point(1, 1)
    @test !(point(1, 1) ≺ point(1, 1))
    @test point(1, 2) ⪯ point(3, 4)
    @test point(1, 2) ≺ point(3, 4)
    @test point(1, 1) ⪰ point(1, 1)
    @test !(point(1, 1) ≻ point(1, 1))
    @test point(3, 4) ⪰ point(1, 2)
    @test point(3, 4) ≻ point(1, 2)

    # center and centroid
    @test Meshes.center(point(1, 1)) == point(1, 1)
    @test centroid(point(1, 1)) == point(1, 1)

    # measure of points is zero
    @test measure(point(1, 2)) == zero(ℳ)
    @test measure(point(1, 2, 3)) == zero(ℳ)

    # boundary of points is nothing
    @test isnothing(boundary(rand(Point{1})))
    @test isnothing(boundary(rand(Point{2})))
    @test isnothing(boundary(rand(Point{3})))

    # check broadcasting works as expected
    @test point(2, 2) .- [point(2, 3), point(3, 1)] == [vec(0.0, -1.0), vec(-1.0, 1.0)]
    @test point(2, 2, 2) .- [point(2, 3, 1), point(3, 1, 4)] == [vec(0.0, -1.0, 1.0), vec(-1.0, 1.0, -2.0)]

    # angles between 2D points
    @test ∠(point(0, 1), point(0, 0), point(1, 0)) ≈ T(-π / 2)
    @test ∠(point(1, 0), point(0, 0), point(0, 1)) ≈ T(π / 2)
    @test ∠(point(-1, 0), point(0, 0), point(0, 1)) ≈ T(-π / 2)
    @test ∠(point(0, 1), point(0, 0), point(-1, 0)) ≈ T(π / 2)
    @test ∠(point(0, -1), point(0, 0), point(1, 0)) ≈ T(π / 2)
    @test ∠(point(1, 0), point(0, 0), point(0, -1)) ≈ T(-π / 2)
    @test ∠(point(0, -1), point(0, 0), point(-1, 0)) ≈ T(-π / 2)
    @test ∠(point(-1, 0), point(0, 0), point(0, -1)) ≈ T(π / 2)

    # angles between 3D points
    @test ∠(point(1, 0, 0), point(0, 0, 0), point(0, 1, 0)) ≈ T(π / 2)
    @test ∠(point(1, 0, 0), point(0, 0, 0), point(0, 0, 1)) ≈ T(π / 2)
    @test ∠(point(0, 1, 0), point(0, 0, 0), point(1, 0, 0)) ≈ T(π / 2)
    @test ∠(point(0, 1, 0), point(0, 0, 0), point(0, 0, 1)) ≈ T(π / 2)
    @test ∠(point(0, 0, 1), point(0, 0, 0), point(1, 0, 0)) ≈ T(π / 2)
    @test ∠(point(0, 0, 1), point(0, 0, 0), point(0, 1, 0)) ≈ T(π / 2)

    # a point pertains to itself
    p = point(0, 0)
    q = point(1, 1)
    @test p ∈ p
    @test q ∈ q
    @test p ∉ q
    @test q ∉ p
    p = point(0, 0, 0)
    q = point(1, 1, 1)
    @test p ∈ p
    @test q ∈ q
    @test p ∉ q
    @test q ∉ p

    p = point(0, 1)
    @test sprint(show, p, context=:compact => true) == "(x: 0.0 m, y: 1.0 m)"
    @test sprint(show, p) == "Point(x: 0.0 m, y: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), p) == """
      Point with Cartesian{NoDatum} coordinates
      ├─ x: 0.0f0 m
      └─ y: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), p) == """
      Point with Cartesian{NoDatum} coordinates
      ├─ x: 0.0 m
      └─ y: 1.0 m"""
    end
  end

  @testset "Ray" begin
    r = Ray(point(0, 0), vec(1, 1))
    @test paramdim(r) == 1
    @test measure(r) == typemax(ℳ)
    @test length(r) == typemax(ℳ)
    @test boundary(r) == point(0, 0)
    @test perimeter(r) == zero(ℳ)

    r = Ray(point(0, 0), vec(1, 1))
    @test r(T(0.0)) == point(0, 0)
    @test r(T(1.0)) == point(1, 1)
    @test r(T(Inf)) == point(Inf, Inf)
    @test r(T(1.0)) - r(T(0.0)) == vec(1, 1)
    @test_throws DomainError(T(-1), "r(t) is not defined for t < 0.") r(T(-1))

    p₁ = point(3, 3, 3)
    p₂ = point(-3, -3, -3)
    p₃ = point(1, 0, 0)
    r = Ray(point(0, 0, 0), vec(1, 1, 1))
    @test p₁ ∈ r
    @test p₂ ∉ r
    @test p₃ ∉ r

    r1 = Ray(point(0, 0, 0), vec(1, 0, 0))
    r2 = Ray(point(1, 1, 1), vec(1, 2, 1))
    @test r1 != r2

    r1 = Ray(point(0, 0, 0), vec(1, 0, 0))
    r2 = Ray(point(1, 0, 0), vec(-1, 0, 0))
    @test r1 != r2

    r1 = Ray(point(0, 0, 0), vec(1, 0, 0))
    r2 = Ray(point(1, 0, 0), vec(1, 0, 0))
    @test r1 != r2

    r1 = Ray(point(0, 0, 0), vec(2, 0, 0))
    r2 = Ray(point(0, 0, 0), vec(1, 0, 0))
    @test r1 == r2

    r2 = rand(Ray{2})
    r3 = rand(Ray{3})
    @test r2 isa Ray
    @test r3 isa Ray
    @test embeddim(r2) == 2
    @test embeddim(r3) == 3

    r = Ray(point(0, 0), vec(1, 1))
    @test sprint(show, r) == "Ray(p: (x: 0.0 m, y: 0.0 m), v: (1.0 m, 1.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), r) == """
      Ray
      ├─ p: Point(x: 0.0 m, y: 0.0 m)
      └─ v: Vec(1.0f0 m, 1.0f0 m)"""
    else
      @test sprint(show, MIME("text/plain"), r) == """
      Ray
      ├─ p: Point(x: 0.0 m, y: 0.0 m)
      └─ v: Vec(1.0 m, 1.0 m)"""
    end
  end

  @testset "Line" begin
    l = Line(point(0, 0), point(1, 1))
    @test paramdim(l) == 1
    @test measure(l) == typemax(ℳ)
    @test length(l) == typemax(ℳ)
    @test isnothing(boundary(l))
    @test perimeter(l) == zero(ℳ)

    l = Line(point(0, 0), point(1, 1))
    @test (l(0), l(1)) == (point(0, 0), point(1, 1))

    l2 = rand(Line{2})
    l3 = rand(Line{3})
    @test l2 isa Line
    @test l3 isa Line
    @test embeddim(l2) == 2
    @test embeddim(l3) == 3

    l = Line(point(0, 0), point(1, 1))
    @test sprint(show, l) == "Line(a: (x: 0.0 m, y: 0.0 m), b: (x: 1.0 m, y: 1.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), l) == """
      Line
      ├─ a: Point(x: 0.0 m, y: 0.0 m)
      └─ b: Point(x: 1.0 m, y: 1.0 m)"""
    else
      @test sprint(show, MIME("text/plain"), l) == """
      Line
      ├─ a: Point(x: 0.0 m, y: 0.0 m)
      └─ b: Point(x: 1.0 m, y: 1.0 m)"""
    end
  end

  @testset "Plane" begin
    p = Plane(point(0, 0, 0), vec(1, 0, 0), vec(0, 1, 0))
    @test p(T(1), T(0)) == point(1, 0, 0)
    @test paramdim(p) == 2
    @test embeddim(p) == 3
    @test measure(p) == typemax(ℳ)^2
    @test area(p) == typemax(ℳ)^2
    @test p(T(0), T(0)) == point(0, 0, 0)
    @test normal(p) == Vec(0, 0, 1)
    @test isnothing(boundary(p))
    @test perimeter(p) == zero(ℳ)

    p = Plane(point(0, 0, 0), vec(0, 0, 1))
    @test p(T(1), T(0)) == point(1, 0, 0)
    @test p(T(0), T(1)) == point(0, 1, 0)

    p₁ = Plane(point(0, 0, 0), vec(1, 0, 0), vec(0, 1, 0))
    p₂ = Plane(point(0, 0, 0), vec(0, 1, 0), vec(1, 0, 0))
    @test p₁ ≈ p₂
    p₁ = Plane(point(0, 0, 0), vec(1, 1, 0))
    p₂ = Plane(point(0, 0, 0), -vec(1, 1, 0))
    @test p₁ ≈ p₂

    # https://github.com/JuliaGeometry/Meshes.jl/issues/624
    p₁ = Plane(point(0, 0, 0), vec(0, 0, 1))
    p₂ = Plane(point(0, 0, 10), vec(0, 0, 1))
    @test !(p₁ ≈ p₂)

    # normal to plane has norm one regardless of basis
    p = Plane(point(0, 0, 0), vec(2, 0, 0), vec(0, 3, 0))
    n = normal(p)
    @test isapprox(norm(n), oneunit(ℳ), atol=atol(ℳ))

    # plane passing through three points
    p₁ = point(0, 0, 0)
    p₂ = point(1, 2, 3)
    p₃ = point(3, 2, 1)
    p = Plane(p₁, p₂, p₃)
    @test p₁ ∈ p
    @test p₂ ∈ p
    @test p₃ ∈ p

    p = rand(Plane)
    @test p isa Plane
    @test embeddim(p) == 3

    p = Plane(point(0, 0, 0), vec(1, 0, 0), vec(0, 1, 0))
    @test sprint(show, p) == "Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, 0.0 m, 0.0 m), v: (0.0 m, 1.0 m, 0.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), p) == """
      Plane
      ├─ p: Point(x: 0.0 m, y: 0.0 m, z: 0.0 m)
      ├─ u: Vec(1.0f0 m, 0.0f0 m, 0.0f0 m)
      └─ v: Vec(0.0f0 m, 1.0f0 m, 0.0f0 m)"""
    else
      @test sprint(show, MIME("text/plain"), p) == """
      Plane
      ├─ p: Point(x: 0.0 m, y: 0.0 m, z: 0.0 m)
      ├─ u: Vec(1.0 m, 0.0 m, 0.0 m)
      └─ v: Vec(0.0 m, 1.0 m, 0.0 m)"""
    end
  end

  @testset "BezierCurve" begin
    b = BezierCurve(point(0, 0), point(0.5, 1), point(1, 0))
    @test embeddim(b) == 2
    @test paramdim(b) == 1

    b = BezierCurve(point(0, 0), point(0.5, 1), point(1, 0))
    for method in [DeCasteljau(), Horner()]
      @test b(T(0), method) == point(0, 0)
      @test b(T(1), method) == point(1, 0)
      @test b(T(0.5), method) == point(0.5, 0.5)
      @test b(T(0.5), method) == point(0.5, 0.5)
      @test_throws DomainError(T(-0.1), "b(t) is not defined for t outside [0, 1].") b(T(-0.1), method)
      @test_throws DomainError(T(1.2), "b(t) is not defined for t outside [0, 1].") b(T(1.2), method)
    end

    @test boundary(b) == Multi([point(0, 0), point(1, 0)])
    b = BezierCurve(point(0, 0), point(1, 1))
    @test boundary(b) == Multi([point(0, 0), point(1, 1)])
    @test perimeter(b) == zero(ℳ)

    b = BezierCurve(point.(randn(100), randn(100)))
    t1 = @timed b(T(0.2))
    t2 = @timed b(T(0.2), Horner())
    @test t1.time > t2.time
    @test t2.bytes < 100

    b2 = rand(BezierCurve{2})
    b3 = rand(BezierCurve{3})
    @test b2 isa BezierCurve
    @test b3 isa BezierCurve
    @test embeddim(b2) == 2
    @test embeddim(b3) == 3

    b = BezierCurve(point(0, 0), point(0.5, 1), point(1, 0))
    if T === Float32
      @test sprint(show, b) == "BezierCurve(controls: Point{2, Cartesian{NoDatum, 2, Quantity{Float32, 𝐋, FreeUnits{(m,), 𝐋, nothing}}}}[(x: 0.0 m, y: 0.0 m), (x: 0.5 m, y: 1.0 m), (x: 1.0 m, y: 0.0 m)])"
      @test sprint(show, MIME("text/plain"), b) == """
      BezierCurve
      └─ controls: Point{2, CoordRefSystems.Cartesian{CoordRefSystems.NoDatum, 2, Quantity{Float32, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}}[Point(x: 0.0 m, y: 0.0 m), Point(x: 0.5 m, y: 1.0 m), Point(x: 1.0 m, y: 0.0 m)]"""
    else
      @test sprint(show, b) == "BezierCurve(controls: Point{2, Cartesian{NoDatum, 2, Quantity{Float64, 𝐋, FreeUnits{(m,), 𝐋, nothing}}}}[(x: 0.0 m, y: 0.0 m), (x: 0.5 m, y: 1.0 m), (x: 1.0 m, y: 0.0 m)])"
      @test sprint(show, MIME("text/plain"), b) == """
      BezierCurve
      └─ controls: Point{2, CoordRefSystems.Cartesian{CoordRefSystems.NoDatum, 2, Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}}[Point(x: 0.0 m, y: 0.0 m), Point(x: 0.5 m, y: 1.0 m), Point(x: 1.0 m, y: 0.0 m)]"""
    end
  end

  @testset "Box" begin
    b = Box(point(0), point(1))
    @test embeddim(b) == 1
    @test paramdim(b) == 1
    @test Meshes.lentype(b) == ℳ
    @test minimum(b) == point(0)
    @test maximum(b) == point(1)
    @test extrema(b) == (point(0), point(1))

    b = Box(point(0, 0), point(1, 1))
    @test embeddim(b) == 2
    @test paramdim(b) == 2
    @test Meshes.lentype(b) == ℳ
    @test minimum(b) == point(0, 0)
    @test maximum(b) == point(1, 1)
    @test extrema(b) == (point(0, 0), point(1, 1))

    b = Box(point(0, 0, 0), point(1, 1, 1))
    @test embeddim(b) == 3
    @test paramdim(b) == 3
    @test Meshes.lentype(b) == ℳ
    @test minimum(b) == point(0, 0, 0)
    @test maximum(b) == point(1, 1, 1)
    @test extrema(b) == (point(0, 0, 0), point(1, 1, 1))

    b = Box(point(0), point(1))
    @test boundary(b) == Multi([point(0), point(1)])
    @test measure(b) == T(1) * u"m"
    @test point(0) ∈ b
    @test point(1) ∈ b
    @test point(0.5) ∈ b
    @test point(-0.5) ∉ b
    @test point(1.5) ∉ b

    b = Box(point(0, 0), point(1, 1))
    @test measure(b) == area(b) == T(1) * u"m^2"
    @test point(1, 1) ∈ b
    @test perimeter(b) ≈ T(4) * u"m"

    b = Box(point(1, 1), point(2, 2))
    @test sides(b) == (T(1) * u"m", T(1) * u"m")
    @test Meshes.center(b) == point(1.5, 1.5)
    @test diagonal(b) == √T(2) * u"m"

    b = Box(point(1, 2), point(3, 4))
    v = point.([(1, 2), (3, 2), (3, 4), (1, 4)])
    @test boundary(b) == Ring(v)

    b = Box(point(1, 2, 3), point(4, 5, 6))
    v = point.([(1, 2, 3), (4, 2, 3), (4, 5, 3), (1, 5, 3), (1, 2, 6), (4, 2, 6), (4, 5, 6), (1, 5, 6)])
    c = connect.([(4, 3, 2, 1), (6, 5, 1, 2), (3, 7, 6, 2), (4, 8, 7, 3), (1, 5, 8, 4), (6, 7, 8, 5)])
    @test boundary(b) == SimpleMesh(v, c)

    b = Box(point(0, 0), point(1, 1))
    @test boundary(b) == Ring(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))

    b = Box(point(0, 0, 0), point(1, 1, 1))
    m = boundary(b)
    @test m isa Mesh
    @test nvertices(m) == 8
    @test nelements(m) == 6

    # subsetting with boxes
    b1 = Box(point(0, 0), point(0.5, 0.5))
    b2 = Box(point(0.1, 0.1), point(0.5, 0.5))
    b3 = Box(point(0, 0), point(1, 1))
    @test b1 ⊆ b3
    @test b2 ⊆ b3
    @test !(b1 ⊆ b2)
    @test !(b3 ⊆ b1)
    @test !(b3 ⊆ b1)

    b = Box(point(0, 0), point(10, 20))
    @test b(T(0.0), T(0.0)) == point(0, 0)
    @test b(T(0.5), T(0.0)) == point(5, 0)
    @test b(T(1.0), T(0.0)) == point(10, 0)
    @test b(T(0.0), T(0.5)) == point(0, 10)
    @test b(T(0.0), T(1.0)) == point(0, 20)

    b = Box(point(0, 0, 0), point(10, 20, 30))
    @test b(T(0.0), T(0.0), T(0.0)) == point(0, 0, 0)
    @test b(T(1.0), T(1.0), T(1.0)) == point(10, 20, 30)

    b1 = rand(Box{1})
    b2 = rand(Box{2})
    b3 = rand(Box{3})
    @test b1 isa Box
    @test b2 isa Box
    @test b3 isa Box
    @test embeddim(b1) == 1
    @test embeddim(b2) == 2
    @test embeddim(b3) == 3

    @test_throws AssertionError Box(point(1), point(0))
    @test_throws AssertionError Box(point(1, 1), point(0, 0))
    @test_throws AssertionError Box(point(1, 1, 1), point(0, 0, 0))

    b = Box(point(0, 0), point(1, 1))
    q = convert(Quadrangle, b)
    @test q isa Quadrangle
    @test q == Quadrangle(point(0, 0), point(1, 0), point(1, 1), point(0, 1))

    b = Box(point(0, 0, 0), point(1, 1, 1))
    h = convert(Hexahedron, b)
    @test h isa Hexahedron
    @test h == Hexahedron(
      point(0, 0, 0),
      point(1, 0, 0),
      point(1, 1, 0),
      point(0, 1, 0),
      point(0, 0, 1),
      point(1, 0, 1),
      point(1, 1, 1),
      point(0, 1, 1)
    )

    b = Box(point(0, 0), point(1, 1))
    @test sprint(show, b) == "Box(min: (x: 0.0 m, y: 0.0 m), max: (x: 1.0 m, y: 1.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), b) == """
      Box
      ├─ min: Point(x: 0.0 m, y: 0.0 m)
      └─ max: Point(x: 1.0 m, y: 1.0 m)"""
    else
      @test sprint(show, MIME("text/plain"), b) == """
      Box
      ├─ min: Point(x: 0.0 m, y: 0.0 m)
      └─ max: Point(x: 1.0 m, y: 1.0 m)"""
    end
  end

  @testset "Ball" begin
    b = Ball(point(1, 2, 3), T(5))
    @test embeddim(b) == 3
    @test paramdim(b) == 3
    @test Meshes.lentype(b) == ℳ
    @test Meshes.center(b) == point(1, 2, 3)
    @test radius(b) == T(5) * u"m"

    b = Ball(point(1, 2, 3), 4)
    @test Meshes.lentype(b) == ℳ

    # TODO: radius must have the same lentype as point in the Ball(point) constructor?
    # b1 = Ball(point(0, 0), T(1))
    # b2 = Ball(point(0, 0))
    # b3 = Ball(T.((0, 0)))
    # @test b1 == b2 == b3

    b = Ball(point(0, 0), T(2))
    @test measure(b) ≈ T(π) * (T(2)^2) * u"m^2"
    b = Ball(point(0, 0, 0), T(2))
    @test measure(b) ≈ T(4 / 3) * T(π) * (T(2)^3) * u"m^3"
    @test_throws ArgumentError length(b)
    @test_throws ArgumentError area(b)

    b = Ball(point(0, 0), T(2))
    @test point(1, 0) ∈ b
    @test point(0, 1) ∈ b
    @test point(2, 0) ∈ b
    @test point(0, 2) ∈ b
    @test point(3, 5) ∉ b
    @test perimeter(b) ≈ T(4π) * u"m"

    b = Ball(point(0, 0, 0), T(2))
    @test point(1, 0, 0) ∈ b
    @test point(0, 0, 1) ∈ b
    @test point(2, 0, 0) ∈ b
    @test point(0, 0, 2) ∈ b
    @test point(3, 5, 2) ∉ b

    b = Ball(point(0, 0), T(2))
    @test b(T(0), T(0)) ≈ point(0, 0)
    @test b(T(1), T(0)) ≈ point(2, 0)

    b = Ball(point(7, 7), T(1.5))
    ps = b.(1, rand(T, 100))
    all(∈(b), ps)

    b = Ball(point(0, 0, 0), T(2))
    @test b(T(0), T(0), T(0)) ≈ point(0, 0, 0)
    @test b(T(1), T(0), T(0)) ≈ point(0, 0, 2)

    b = Ball(point(7, 7, 7), T(1.5))
    ps = b.(1, rand(T, 100), rand(T, 100))
    all(∈(b), ps)

    b1 = rand(Ball{1})
    b2 = rand(Ball{2})
    b3 = rand(Ball{3})
    @test b1 isa Ball
    @test b2 isa Ball
    @test b3 isa Ball
    @test embeddim(b1) == 1
    @test embeddim(b2) == 2
    @test embeddim(b3) == 3

    b = Ball(point(0, 0), T(1))
    @test sprint(show, b) == "Ball(center: (x: 0.0 m, y: 0.0 m), radius: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), b) == """
      Ball
      ├─ center: Point(x: 0.0 m, y: 0.0 m)
      └─ radius: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), b) == """
      Ball
      ├─ center: Point(x: 0.0 m, y: 0.0 m)
      └─ radius: 1.0 m"""
    end
  end

  @testset "Sphere" begin
    s = Sphere(point(0, 0, 0), T(1))
    @test embeddim(s) == 3
    @test paramdim(s) == 2
    @test Meshes.lentype(s) == ℳ
    @test Meshes.center(s) == point(0, 0, 0)
    @test radius(s) == T(1) * u"m"
    @test extrema(s) == (point(-1, -1, -1), point(1, 1, 1))
    @test isnothing(boundary(s))
    @test perimeter(s) == zero(ℳ)

    s = Sphere(point(1, 2, 3), 4)
    @test Meshes.lentype(s) == ℳ

    s = Sphere(point(0, 0), T(1))
    @test embeddim(s) == 2
    @test paramdim(s) == 1
    @test Meshes.lentype(s) == ℳ
    @test Meshes.center(s) == point(0, 0)
    @test radius(s) == T(1) * u"m"
    @test extrema(s) == (point(-1, -1), point(1, 1))
    @test isnothing(boundary(s))

    # TODO: radius must have the same lentype as point in the Sphere(point) constructor?
    # s1 = Sphere(point(0, 0), T(1))
    # s2 = Sphere(point(0, 0))
    # s3 = Sphere(T.((0, 0)))
    # @test s1 == s2 == s3

    s = Sphere(point(0, 0), T(2))
    @test measure(s) ≈ T(2π) * 2 * u"m"
    @test length(s) ≈ T(2π) * 2 * u"m"
    @test extrema(s) == (point(-2, -2), point(2, 2))
    s = Sphere(point(0, 0, 0), T(2))
    @test measure(s) ≈ T(4π) * (2^2) * u"m^2"
    @test area(s) ≈ T(4π) * (2^2) * u"m^2"

    s = Sphere(point(0, 0), T(2))
    @test point(1, 0) ∉ s
    @test point(0, 1) ∉ s
    @test point(2, 0) ∈ s
    @test point(0, 2) ∈ s
    @test point(3, 5) ∉ s

    s = Sphere(point(0, 0, 0), T(2))
    @test point(1, 0, 0) ∉ s
    @test point(0, 0, 1) ∉ s
    @test point(2, 0, 0) ∈ s
    @test point(0, 0, 2) ∈ s
    @test point(3, 5, 2) ∉ s

    # 2D sphere passing through 3 points
    s = Sphere(point(0, 0), point(0.5, 0), point(1, 1))
    @test Meshes.center(s) == point(0.25, 0.75)
    @test radius(s) == T(0.7905694150420949) * u"m"
    s = Sphere(point(0, 0), point(1, 0), point(0, 1))
    @test Meshes.center(s) == point(0.5, 0.5)
    @test radius(s) == T(0.7071067811865476) * u"m"
    s = Sphere(point(0, 0), point(1, 0), point(1, 1))
    @test Meshes.center(s) == point(0.5, 0.5)
    @test radius(s) == T(0.7071067811865476) * u"m"

    # 3D sphere passing through 4 points
    s = Sphere(point(0, 0, 0), point(5, 0, 1), point(1, 1, 1), point(3, 2, 1))
    @test point(0, 0, 0) ∈ s
    @test point(5, 0, 1) ∈ s
    @test point(1, 1, 1) ∈ s
    @test point(3, 2, 1) ∈ s
    O = Meshes.center(s)
    r = radius(s)
    @test isapprox(r, norm(point(0, 0, 0) - O))

    s = Sphere(point(0, 0), T(2))
    @test s(T(0)) ≈ point(2, 0)
    @test s(T(0.5)) ≈ point(-2, 0)

    s = Sphere(point(0, 0, 0), T(2))
    @test s(T(0), T(0)) ≈ point(0, 0, 2)
    @test s(T(0.5), T(0.5)) ≈ point(-2, 0, 0)

    s1 = rand(Sphere{1})
    s2 = rand(Sphere{2})
    s3 = rand(Sphere{3})
    @test s1 isa Sphere
    @test s2 isa Sphere
    @test s3 isa Sphere
    @test embeddim(s1) == 1
    @test embeddim(s2) == 2
    @test embeddim(s3) == 3

    s = Sphere(point(0, 0, 0), T(1))
    @test sprint(show, s) == "Sphere(center: (x: 0.0 m, y: 0.0 m, z: 0.0 m), radius: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), s) == """
      Sphere
      ├─ center: Point(x: 0.0 m, y: 0.0 m, z: 0.0 m)
      └─ radius: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), s) == """
      Sphere
      ├─ center: Point(x: 0.0 m, y: 0.0 m, z: 0.0 m)
      └─ radius: 1.0 m"""
    end
  end

  @testset "Ellipsoid" begin
    e = Ellipsoid((T(3), T(2), T(1)), point(0, 0, 0))
    @test embeddim(e) == 3
    @test paramdim(e) == 2
    @test Meshes.lentype(e) == ℳ
    @test radii(e) == (T(3) * u"m", T(2) * u"m", T(1) * u"m")
    @test center(e) == point(0, 0, 0)
    @test isnothing(boundary(e))
    @test perimeter(e) == zero(ℳ)

    e = Ellipsoid((T(3), T(2), T(1)), point(0, 0, 0))
    @test sprint(show, e) ==
          "Ellipsoid(radii: (3.0 m, 2.0 m, 1.0 m), center: (x: 0.0 m, y: 0.0 m, z: 0.0 m), rotation: UniformScaling{Bool}(true))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), e) == """
      Ellipsoid
      ├─ radii: (3.0f0 m, 2.0f0 m, 1.0f0 m)
      ├─ center: Point(x: 0.0 m, y: 0.0 m, z: 0.0 m)
      └─ rotation: UniformScaling{Bool}(true)"""
    else
      @test sprint(show, MIME("text/plain"), e) == """
      Ellipsoid
      ├─ radii: (3.0 m, 2.0 m, 1.0 m)
      ├─ center: Point(x: 0.0 m, y: 0.0 m, z: 0.0 m)
      └─ rotation: UniformScaling{Bool}(true)"""
    end
  end

  @testset "Disk" begin
    p = Plane(point(0, 0, 0), vec(0, 0, 1))
    d = Disk(p, T(2))
    @test embeddim(d) == 3
    @test paramdim(d) == 2
    @test Meshes.lentype(d) == ℳ
    @test plane(d) == p
    @test Meshes.center(d) == point(0, 0, 0)
    @test radius(d) == T(2) * u"m"
    @test normal(d) == vec(0, 0, 1)
    @test measure(d) == T(π) * T(2)^2 * u"m^2"
    @test area(d) == measure(d)
    @test point(0, 0, 0) ∈ d
    @test point(0, 0, 1) ∉ d
    @test boundary(d) == Circle(p, T(2))

    d = rand(Disk)
    @test d isa Disk
    @test embeddim(d) == 3

    p = Plane(point(0, 0, 0), vec(0, 0, 1))
    d = Disk(p, T(2))
    @test sprint(show, d) ==
          "Disk(plane: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m)), radius: 2.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), d) == """
      Disk
      ├─ plane: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m))
      └─ radius: 2.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), d) == """
      Disk
      ├─ plane: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m))
      └─ radius: 2.0 m"""
    end
  end

  @testset "Circle" begin
    p = Plane(point(0, 0, 0), vec(0, 0, 1))
    c = Circle(p, T(2))
    @test embeddim(c) == 3
    @test paramdim(c) == 1
    @test Meshes.lentype(c) == ℳ
    @test plane(c) == p
    @test Meshes.center(c) == point(0, 0, 0)
    @test radius(c) == T(2) * u"m"
    @test measure(c) == 2 * T(π) * T(2) * u"m"
    @test length(c) == measure(c)
    @test point(2, 0, 0) ∈ c
    @test point(0, 2, 0) ∈ c
    @test point(0, 0, 0) ∉ c
    @test isnothing(boundary(c))

    # 3D circumcircle
    p1 = point(0, 4, 0)
    p2 = point(0, -4, 0)
    p3 = point(0, 0, 4)
    c = Circle(p1, p2, p3)
    @test p1 ∈ c
    @test p2 ∈ c
    @test p3 ∈ c

    # circle parametrization
    p = Plane(point(0, 0, 0), vec(0, 0, 1))
    c = Circle(p, T(2))
    @test c(T(0)) ≈ point(2, 0, 0)
    @test c(T(0.25)) ≈ point(0, 2, 0)
    @test c(T(0.5)) ≈ point(-2, 0, 0)
    @test c(T(0.75)) ≈ point(0, -2, 0)
    @test c(T(1)) ≈ point(2, 0, 0)

    c = rand(Circle)
    @test c isa Circle
    @test embeddim(c) == 3

    p = Plane(point(0, 0, 0), vec(0, 0, 1))
    c = Circle(p, T(2))
    @test sprint(show, c) ==
          "Circle(plane: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m)), radius: 2.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Circle
      ├─ plane: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m))
      └─ radius: 2.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Circle
      ├─ plane: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m))
      └─ radius: 2.0 m"""
    end
  end

  @testset "Cylinder" begin
    c = Cylinder(Plane(point(1, 2, 3), vec(0, 0, 1)), Plane(point(4, 5, 6), vec(0, 0, 1)), T(5))
    @test embeddim(c) == 3
    @test paramdim(c) == 3
    @test Meshes.lentype(c) == ℳ
    @test radius(c) == T(5) * u"m"
    @test bottom(c) == Plane(point(1, 2, 3), vec(0, 0, 1))
    @test top(c) == Plane(point(4, 5, 6), vec(0, 0, 1))
    @test axis(c) == Line(point(1, 2, 3), point(4, 5, 6))
    @test !isright(c)
    @test measure(c) == volume(c) ≈ T(5)^2 * pi * T(3) * sqrt(T(3)) * u"m^3"
    @test point(1, 2, 3) ∈ c
    @test point(4, 5, 6) ∈ c
    @test point(0.99, 1.99, 2.99) ∉ c
    @test point(4.01, 5.01, 6.01) ∉ c
    # TODO: fix intersection(f, Plane, Plane)
    # @test !Meshes.hasintersectingplanes(c)
    @test c(0, 0, 0) ≈ bottom(c)(0, 0)
    @test c(0, 0, 1) ≈ top(c)(0, 0)
    @test c(1, 0.25, 0.5) ≈ Point(T(4.330127018922193), T(10.330127018922191), T(4.5))
    @test_throws DomainError c(1.1, 0, 0)

    # c = Cylinder(Plane(point(0, 0, 0), vec(0, 0, 1)), Plane(point(0, 0, 1), vec(1, 0, 1)), T(5))
    # @test Meshes.hasintersectingplanes(c)

    # TODO: radius must have the same lentype as point in the Cylinder(point, point, radius) constructor?
    # c1 = Cylinder(point(0, 0, 0), point(0, 0, 1), T(1))
    # c2 = Cylinder(point(0, 0, 0), point(0, 0, 1))
    # c3 = Cylinder(T(1))
    # @test c1 == c2 == c3
    # @test c1 ≈ c2 ≈ c3

    c = Cylinder(T(1))
    @test Meshes.lentype(c) == Meshes.Met{Float64}
    c = Cylinder(1)
    @test Meshes.lentype(c) == Meshes.Met{Float64}

    c = Cylinder(point(0, 0, 0), point(0, 0, 1), T(1))
    @test radius(c) == T(1) * u"m"
    @test bottom(c) == Plane(point(0, 0, 0), vec(0, 0, 1))
    @test top(c) == Plane(point(0, 0, 1), vec(0, 0, 1))
    @test center(c) == point(0.0, 0.0, 0.5)
    @test centroid(c) == point(0.0, 0.0, 0.5)
    @test axis(c) == Line(point(0, 0, 0), point(0, 0, 1))
    @test isright(c)
    @test boundary(c) == CylinderSurface(point(0, 0, 0), point(0, 0, 1), T(1))
    @test measure(c) == volume(c) ≈ T(π) * u"m^3"
    @test point(0, 0, 0) ∈ c
    @test point(0, 0, 1) ∈ c
    @test point(1, 0, 0) ∈ c
    @test point(0, 1, 0) ∈ c
    @test point(cosd(60), sind(60), 0.5) ∈ c
    @test point(0, 0, -0.001) ∉ c
    @test point(0, 0, 1.001) ∉ c
    @test point(1, 1, 1) ∉ c

    c = rand(Cylinder)
    @test c isa Cylinder
    @test embeddim(c) == 3

    c = Cylinder(point(0, 0, 0), point(0, 0, 1), T(1))
    @test sprint(show, c) ==
          "Cylinder(bot: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m)), top: Plane(p: (x: 0.0 m, y: 0.0 m, z: 1.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m)), radius: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Cylinder
      ├─ bot: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m))
      ├─ top: Plane(p: (x: 0.0 m, y: 0.0 m, z: 1.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m))
      └─ radius: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Cylinder
      ├─ bot: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m))
      ├─ top: Plane(p: (x: 0.0 m, y: 0.0 m, z: 1.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m))
      └─ radius: 1.0 m"""
    end
  end

  @testset "CylinderSurface" begin
    c = CylinderSurface(point(0, 0, 0), point(0, 0, 1), T(2))
    @test embeddim(c) == 3
    @test paramdim(c) == 2
    @test Meshes.lentype(c) == ℳ
    @test radius(c) == T(2) * u"m"
    @test bottom(c) == Plane(point(0, 0, 0), vec(0, 0, 1))
    @test top(c) == Plane(point(0, 0, 1), vec(0, 0, 1))
    @test center(c) == point(0.0, 0.0, 0.5)
    @test centroid(c) == point(0.0, 0.0, 0.5)
    @test axis(c) == Line(point(0, 0, 0), point(0, 0, 1))
    @test isright(c)
    @test isnothing(boundary(c))
    @test measure(c) == area(c) ≈ (2 * T(2)^2 * pi + 2 * T(2) * pi) * u"m^2"
    # TODO: fix intersection(f, Plane, Plane)
    # @test !Meshes.hasintersectingplanes(c)

    # c = CylinderSurface(Plane(point(0, 0, 0), vec(0, 0, 1)), Plane(point(0, 0, 1), vec(1, 0, 1)), T(5))
    # @test Meshes.hasintersectingplanes(c)

    # # TODO: radius must have the same lentype as point in the CylinderSurface(point, point, radius) constructor?
    # c1 = CylinderSurface(point(0, 0, 0), point(0, 0, 1), T(1))
    # c2 = CylinderSurface(point(0, 0, 0), point(0, 0, 1))
    # c3 = CylinderSurface(T(1))
    # @test c1 == c2 == c3
    # @test c1 ≈ c2 ≈ c3

    c = CylinderSurface(Plane(point(1, 2, 3), vec(0, 0, 1)), Plane(point(4, 5, 6), vec(0, 0, 1)), T(5))
    @test measure(c) == area(c) ≈ (2 * T(5)^2 * pi + 2 * T(5) * pi * sqrt(3 * T(3)^2)) * u"m^2"

    c = CylinderSurface(point(0, 0, 0), point(0, 0, 1), T(1))
    @test c(T(0), T(0)) ≈ point(1, 0, 0)
    @test c(T(0.5), T(0)) ≈ point(-1, 0, 0)
    @test c(T(0), T(1)) ≈ point(1, 0, 1)
    @test c(T(0.5), T(1)) ≈ point(-1, 0, 1)

    c = CylinderSurface(1.0)
    @test Meshes.lentype(c) == Meshes.Met{Float64}
    c = CylinderSurface(1.0f0)
    @test Meshes.lentype(c) == Meshes.Met{Float64}
    c = CylinderSurface(1)
    @test Meshes.lentype(c) == Meshes.Met{Float64}

    c = rand(CylinderSurface)
    @test c isa CylinderSurface
    @test embeddim(c) == 3

    c = CylinderSurface(T(1))
    @test sprint(show, c) ==
          "CylinderSurface(bot: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m)), top: Plane(p: (x: 0.0 m, y: 0.0 m, z: 1.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m)), radius: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      CylinderSurface
      ├─ bot: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m))
      ├─ top: Plane(p: (x: 0.0 m, y: 0.0 m, z: 1.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m))
      └─ radius: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      CylinderSurface
      ├─ bot: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m))
      ├─ top: Plane(p: (x: 0.0 m, y: 0.0 m, z: 1.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m))
      └─ radius: 1.0 m"""
    end
  end

  @testset "ParaboloidSurface" begin
    p = ParaboloidSurface(point(0, 0, 0), T(1), T(2))
    @test embeddim(p) == 3
    @test paramdim(p) == 2
    @test Meshes.lentype(p) == ℳ
    @test focallength(p) == T(2) * u"m"
    @test radius(p) == T(1) * u"m"
    @test axis(p) == Line(point(0, 0, 0), point(0, 0, T(2)))
    @test measure(p) == area(p) ≈ T(32π / 3 * (17√17 / 64 - 1)) * u"m^2"

    # TODO: radius must have the same lentype as apex in the ParaboloidSurface(apex, radius, focallength) constructor?
    # p1 = ParaboloidSurface(point(1, 2, 3), T(1), T(1))
    # p2 = ParaboloidSurface(point(1, 2, 3), T(1))
    # p3 = ParaboloidSurface(point(1, 2, 3))
    # @test p1 == p2 == p3
    # @test p1 ≈ p2 ≈ p3

    p1 = ParaboloidSurface((1, 2, 3), 1.0, 1.0)
    p2 = ParaboloidSurface((1, 2, 3), 1.0)
    p3 = ParaboloidSurface((1, 2, 3))
    @test p1 == p2 == p3
    @test p1 ≈ p2 ≈ p3

    p = ParaboloidSurface((1.0, 2.0, 3.0), 4.0, 5.0)
    @test Meshes.lentype(p) == Meshes.Met{Float64}
    @test radius(p) == 4.0 * u"m"
    @test focallength(p) == 5.0 * u"m"

    p = ParaboloidSurface(point(1, 5, 2), T(3), T(4))
    @test measure(p) == area(p) ≈ T(128π / 3 * (73√73 / 512 - 1)) * u"m^2"
    @test p(T(0), T(0)) ≈ point(1, 5, 2)
    @test p(T(1), T(0)) ≈ point(4, 5, 2 + 3^2 / (4 * 4))
    @test_throws DomainError p(T(-0.1), T(0))
    @test_throws DomainError p(T(1.1), T(0))

    p = ParaboloidSurface()
    @test Meshes.lentype(p) == Meshes.Met{Float64}
    @test p(0.0, 0.0) ≈ Point(0, 0, 0)
    @test p(0.5, 0.0) ≈ Point(0.5, 0, 0.5^2 / 4)
    @test p(0.0, 0.5) ≈ Point(0, 0, 0)
    @test p(0.5, 0.5) ≈ Point(-0.5, 0, 0.5^2 / 4)

    p = ParaboloidSurface(Point(0.0, 0.0, 0.0))
    @test Meshes.lentype(p) == Meshes.Met{Float64}
    p = ParaboloidSurface(Point(0f0, 0f0, 0f0))
    @test Meshes.lentype(p) == Meshes.Met{Float32}

    p = rand(ParaboloidSurface)
    @test p isa ParaboloidSurface
    @test embeddim(p) == 3

    p = ParaboloidSurface(point(0, 0, 0), T(1), T(1))
    @test sprint(show, p) == "ParaboloidSurface(apex: (x: 0.0 m, y: 0.0 m, z: 0.0 m), radius: 1.0 m, focallength: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), p) == """
      ParaboloidSurface
      ├─ apex: Point(x: 0.0 m, y: 0.0 m, z: 0.0 m)
      ├─ radius: 1.0f0 m
      └─ focallength: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), p) == """
      ParaboloidSurface
      ├─ apex: Point(x: 0.0 m, y: 0.0 m, z: 0.0 m)
      ├─ radius: 1.0 m
      └─ focallength: 1.0 m"""
    end
  end

  @testset "Cone" begin
    p = Plane(point(0, 0, 0), vec(0, 0, 1))
    d = Disk(p, T(2))
    a = point(0, 0, 1)
    c = Cone(d, a)
    @test embeddim(c) == 3
    @test paramdim(c) == 3
    @test Meshes.lentype(c) == ℳ
    @test boundary(c) == ConeSurface(d, a)

    c = rand(Cone)
    @test c isa Cone
    @test embeddim(c) == 3

    p = Plane(point(0, 0, 0), vec(0, 0, 1))
    d = Disk(p, T(2))
    a = point(0, 0, 1)
    c = Cone(d, a)
    @test sprint(show, c) ==
          "Cone(base: Disk(plane: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m)), radius: 2.0 m), apex: (x: 0.0 m, y: 0.0 m, z: 1.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Cone
      ├─ base: Disk(plane: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m)), radius: 2.0 m)
      └─ apex: Point(x: 0.0 m, y: 0.0 m, z: 1.0 m)"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Cone
      ├─ base: Disk(plane: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m)), radius: 2.0 m)
      └─ apex: Point(x: 0.0 m, y: 0.0 m, z: 1.0 m)"""
    end

    # cone: apex at (5,4,3); base center at (5,1,3)
    # halfangle: 30° -> radius: sqrt(3)
    # axis of the cone is parallel to y axis
    p = Plane(point(5, 1, 3), vec(0, 1, 0))
    d = Disk(p, sqrt(T(3)))
    a = point(5, 4, 3)
    c = Cone(d, a)

    @test rad2deg(Meshes.halfangle(c)) ≈ T(30)
    @test Meshes.height(c) ≈ T(3) * u"m"

    @test point(5, 1, 3) ∈ c
    @test point(5, 4, 3) ∈ c
    @test point(5, 1, 3 - sqrt(3)) ∈ c
    @test point(5, 1, 3 + sqrt(3)) ∈ c
    @test point(5 - sqrt(3), 1, 3) ∈ c
    @test point(5 + sqrt(3), 1, 3) ∈ c
    @test point(5, 2.5, 3) ∈ c
    @test point(5 + sqrt(3) / 2, 2.5, 3) ∈ c
    @test point(5 - sqrt(3) / 2, 2.5, 3) ∈ c

    @test point(5, 0.9, 3) ∉ c
    @test point(5, 4.1, 3) ∉ c
    @test point(5, 1, 1) ∉ c
    @test point(5 + sqrt(3) + 0.01, 1, 3) ∉ c
    @test point(5 + sqrt(3) / 2 + 0.01, 2.5, 3) ∉ c
    @test point(5 - sqrt(3) / 2 - 0.01, 2.5, 3) ∉ c
  end

  @testset "ConeSurface" begin
    p = Plane(point(0, 0, 0), vec(0, 0, 1))
    d = Disk(p, T(2))
    a = point(0, 0, 1)
    s = ConeSurface(d, a)
    @test embeddim(s) == 3
    @test paramdim(s) == 2
    @test Meshes.lentype(s) == ℳ
    @test isnothing(boundary(s))

    c = rand(ConeSurface)
    @test c isa ConeSurface
    @test embeddim(c) == 3

    p = Plane(point(0, 0, 0), vec(0, 0, 1))
    d = Disk(p, T(2))
    a = point(0, 0, 1)
    s = ConeSurface(d, a)
    @test sprint(show, s) ==
          "ConeSurface(base: Disk(plane: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m)), radius: 2.0 m), apex: (x: 0.0 m, y: 0.0 m, z: 1.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), s) == """
      ConeSurface
      ├─ base: Disk(plane: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m)), radius: 2.0 m)
      └─ apex: Point(x: 0.0 m, y: 0.0 m, z: 1.0 m)"""
    else
      @test sprint(show, MIME("text/plain"), s) == """
      ConeSurface
      ├─ base: Disk(plane: Plane(p: (x: 0.0 m, y: 0.0 m, z: 0.0 m), u: (1.0 m, -0.0 m, -0.0 m), v: (-0.0 m, 1.0 m, -0.0 m)), radius: 2.0 m)
      └─ apex: Point(x: 0.0 m, y: 0.0 m, z: 1.0 m)"""
    end
  end

  @testset "Frustum" begin
    pb = Plane(point(0, 0, 0), vec(0, 0, 1))
    db = Disk(pb, T(1))
    pt = Plane(point(0, 0, 10), vec(0, 0, 1))
    dt = Disk(pt, T(2))
    f = Frustum(db, dt)
    @test embeddim(f) == 3
    @test Meshes.lentype(f) == ℳ
    @test boundary(f) == FrustumSurface(db, dt)

    @test_throws AssertionError Frustum(db, db)

    f = rand(Frustum)
    @test f isa Frustum

    f = Frustum(db, dt)
    @test point(0, 0, 0) ∈ f
    @test point(0, 0, 10) ∈ f
    @test point(1, 0, 0) ∈ f
    @test point(2, 0, 10) ∈ f
    @test point(1, 0, 5) ∈ f

    @test point(1, 1, 0) ∉ f
    @test point(2, 2, 10) ∉ f
    @test point(0, 0, -0.01) ∉ f
    @test point(0, 0, 10.01) ∉ f

    # reverse order, when top is larger than bottom
    # the frustum is the same geometry
    f = Frustum(dt, db)
    @test point(0, 0, 0) ∈ f
    @test point(0, 0, 10) ∈ f
    @test point(1, 0, 0) ∈ f
    @test point(2, 0, 10) ∈ f
    @test point(1, 0, 5) ∈ f

    @test point(1, 1, 0) ∉ f
    @test point(2, 2, 10) ∉ f
    @test point(0, 0, -0.01) ∉ f
    @test point(0, 0, 10.01) ∉ f
  end

  @testset "FrustumSurface" begin
    pb = Plane(point(0, 0, 0), vec(0, 0, 1))
    db = Disk(pb, T(1))
    pt = Plane(point(0, 0, 10), vec(0, 0, 1))
    dt = Disk(pt, T(2))
    f = FrustumSurface(db, dt)
    @test embeddim(f) == 3
    @test paramdim(f) == 2
    @test Meshes.lentype(f) == ℳ
    @test isnothing(boundary(f))

    @test_throws AssertionError FrustumSurface(db, db)

    f = rand(FrustumSurface)
    @test f isa FrustumSurface
  end

  @testset "Torus" begin
    t = Torus(T.((1, 1, 1)), T.((1, 0, 0)), 2, 1)
    @test point(1, 1, -1) ∈ t
    @test point(1, 1, 1) ∉ t
    @test paramdim(t) == 2
    @test Meshes.center(t) == point(1, 1, 1)
    @test normal(t) == vec(1, 0, 0)
    @test radii(t) == (T(2) * u"m", T(1) * u"m")
    @test axis(t) == Line(point(1, 1, 1), point(2, 1, 1))
    @test measure(t) ≈ 8 * T(π)^2 * u"m^2"
    @test_throws ArgumentError length(t)
    @test_throws ArgumentError volume(t)

    # torus passing through three points
    p₁ = point(0, 0, 0)
    p₂ = point(1, 2, 3)
    p₃ = point(3, 2, 1)
    t = Torus(p₁, p₂, p₃, T(1))
    c = center(t)
    R, r = radii(t)
    @test r == T(1) * u"m"
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

    t = rand(Torus)
    @test t isa Torus
    @test embeddim(t) == 3
    @test Meshes.lentype(t) == Meshes.Met{Float64}
    @test isnothing(boundary(t))

    t = Torus(point(1, 1, 1), vec(1, 0, 0), T(2), T(1))
    @test sprint(show, t) == "Torus(center: (x: 1.0 m, y: 1.0 m, z: 1.0 m), normal: (1.0 m, 0.0 m, 0.0 m), major: 2.0 m, minor: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), t) == """
      Torus
      ├─ center: Point(x: 1.0 m, y: 1.0 m, z: 1.0 m)
      ├─ normal: Vec(1.0f0 m, 0.0f0 m, 0.0f0 m)
      ├─ major: 2.0f0 m
      └─ minor: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), t) == """
      Torus
      ├─ center: Point(x: 1.0 m, y: 1.0 m, z: 1.0 m)
      ├─ normal: Vec(1.0 m, 0.0 m, 0.0 m)
      ├─ major: 2.0 m
      └─ minor: 1.0 m"""
    end
  end
end
