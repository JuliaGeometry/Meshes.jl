@testset "Coordinates" begin
  @testset "Cartesian" begin
    @test Cartesian(T(1)) == Cartesian(T(1) * u"m")
    @test Cartesian(T(1), T(1)) == Cartesian(T(1) * u"m", T(1) * u"m")
    @test Cartesian(T(1), T(1), T(1)) == Cartesian(T(1) * u"m", T(1) * u"m", T(1) * u"m")
    @test Cartesian(T(1) * u"m", 1 * u"m") == Cartesian(T(1) * u"m", T(1) * u"m")

    c = Cartesian(T(1))
    @test sprint(show, c) == "Cartesian(x: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Cartesian coordinates
      └─ x: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Cartesian coordinates
      └─ x: 1.0 m"""
    end

    c = Cartesian(T(1), T(1))
    @test sprint(show, c) == "Cartesian(x: 1.0 m, y: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Cartesian coordinates
      ├─ x: 1.0f0 m
      └─ y: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Cartesian coordinates
      ├─ x: 1.0 m
      └─ y: 1.0 m"""
    end

    c = Cartesian(T(1), T(1), T(1))
    @test sprint(show, c) == "Cartesian(x: 1.0 m, y: 1.0 m, z: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Cartesian coordinates
      ├─ x: 1.0f0 m
      ├─ y: 1.0f0 m
      └─ z: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Cartesian coordinates
      ├─ x: 1.0 m
      ├─ y: 1.0 m
      └─ z: 1.0 m"""
    end

    c = Cartesian(T(1), T(1), T(1), T(1))
    @test sprint(show, c) == "Cartesian(x1: 1.0 m, x2: 1.0 m, x3: 1.0 m, x4: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Cartesian coordinates
      ├─ x1: 1.0f0 m
      ├─ x2: 1.0f0 m
      ├─ x3: 1.0f0 m
      └─ x4: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Cartesian coordinates
      ├─ x1: 1.0 m
      ├─ x2: 1.0 m
      ├─ x3: 1.0 m
      └─ x4: 1.0 m"""
    end
  end

  @testset "Polar" begin
    @test Polar(T(1), T(1)) == Polar(T(1) * u"m", T(1) * u"rad")
    @test Polar(T(1) * u"m", T(45) * u"°") ≈ Polar(T(1) * u"m", T(π / 4) * u"rad")

    c = Polar(T(1), T(1))
    @test sprint(show, c) == "Polar(ρ: 1.0 m, ϕ: 1.0 rad)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Polar coordinates
      ├─ ρ: 1.0f0 m
      └─ ϕ: 1.0f0 rad"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Polar coordinates
      ├─ ρ: 1.0 m
      └─ ϕ: 1.0 rad"""
    end
  end

  @testset "Cylindrical" begin
    @test Cylindrical(T(1), T(1), T(1)) == Cylindrical(T(1) * u"m", T(1) * u"rad", T(1) * u"m")
    @test Cylindrical(T(1) * u"m", T(1) * u"rad", 1 * u"m") == Cylindrical(T(1) * u"m", T(1) * u"rad", T(1) * u"m")
    @test Cylindrical(T(1) * u"m", T(45) * u"°", T(1) * u"m") ≈ Cylindrical(T(1) * u"m", T(π / 4) * u"rad", T(1) * u"m")

    c = Cylindrical(T(1), T(1), T(1))
    @test sprint(show, c) == "Cylindrical(ρ: 1.0 m, ϕ: 1.0 rad, z: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Cylindrical coordinates
      ├─ ρ: 1.0f0 m
      ├─ ϕ: 1.0f0 rad
      └─ z: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Cylindrical coordinates
      ├─ ρ: 1.0 m
      ├─ ϕ: 1.0 rad
      └─ z: 1.0 m"""
    end
  end

  @testset "Spherical" begin
    @test Spherical(T(1), T(1), T(1)) == Spherical(T(1) * u"m", T(1) * u"rad", T(1) * u"rad")
    @test Spherical(T(1) * u"m", T(1) * u"rad", 1 * u"rad") == Spherical(T(1) * u"m", T(1) * u"rad", T(1) * u"rad")
    @test Spherical(T(1) * u"m", T(45) * u"°", T(45) * u"°") ≈ Spherical(T(1) * u"m", T(π / 4) * u"rad", T(π / 4) * u"rad")

    c = Spherical(T(1), T(1), T(1))
    @test sprint(show, c) == "Spherical(r: 1.0 m, θ: 1.0 rad, ϕ: 1.0 rad)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Spherical coordinates
      ├─ r: 1.0f0 m
      ├─ θ: 1.0f0 rad
      └─ ϕ: 1.0f0 rad"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Spherical coordinates
      ├─ r: 1.0 m
      ├─ θ: 1.0 rad
      └─ ϕ: 1.0 rad"""
    end
  end

  @testset "LatLon" begin
    @test LatLon(T(1), T(1)) == LatLon(T(1) * u"°", T(1) * u"°")
    @test LatLon(T(1) * u"°", 1 * u"°") == LatLon(T(1) * u"°", T(1) * u"°")
    @test LatLon(T(π / 4) * u"rad", T(π / 4) * u"rad") ≈ LatLon(T(45) * u"°", T(45) * u"°")

    c = LatLon(T(1), T(1))
    @test sprint(show, c) == "LatLon(lat: 1.0°, lon: 1.0°)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      LatLon coordinates
      ├─ lat: 1.0f0°
      └─ lon: 1.0f0°"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      LatLon coordinates
      ├─ lat: 1.0°
      └─ lon: 1.0°"""
    end
  end

  @testset "LatLonAlt" begin
    @test LatLonAlt(T(1), T(1), T(1)) == LatLonAlt(T(1) * u"°", T(1) * u"°", T(1) * u"m")
    @test LatLonAlt(T(1) * u"°", 1 * u"°", T(1) * u"m") == LatLonAlt(T(1) * u"°", T(1) * u"°", T(1) * u"m")
    @test LatLonAlt(T(π / 4) * u"rad", T(π / 4) * u"rad", T(1) * u"m") ≈ LatLonAlt(T(45) * u"°", T(45) * u"°", T(1) * u"m")

    c = LatLonAlt(T(1), T(1), T(1))
    @test sprint(show, c) == "LatLonAlt(lat: 1.0°, lon: 1.0°, alt: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      LatLonAlt coordinates
      ├─ lat: 1.0f0°
      ├─ lon: 1.0f0°
      └─ alt: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      LatLonAlt coordinates
      ├─ lat: 1.0°
      ├─ lon: 1.0°
      └─ alt: 1.0 m"""
    end
  end

  @testset "EastNorth" begin
    @test EastNorth(T(1), T(1)) == EastNorth(T(1) * u"m", T(1) * u"m")
    @test EastNorth(T(1) * u"m", 1 * u"m") == EastNorth(T(1) * u"m", T(1) * u"m")

    c = EastNorth(T(1), T(1))
    @test sprint(show, c) == "EastNorth(east: 1.0 m, north: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      EastNorth coordinates
      ├─ east: 1.0f0 m
      └─ north: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      EastNorth coordinates
      ├─ east: 1.0 m
      └─ north: 1.0 m"""
    end
  end

  @testset "WebMercator" begin
    @test WebMercator(T(1), T(1)) == WebMercator(T(1) * u"m", T(1) * u"m")
    @test WebMercator(T(1) * u"m", 1 * u"m") == WebMercator(T(1) * u"m", T(1) * u"m")

    c = WebMercator(T(1), T(1))
    @test sprint(show, c) == "WebMercator(x: 1.0 m, y: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      WebMercator coordinates
      ├─ x: 1.0f0 m
      └─ y: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      WebMercator coordinates
      ├─ x: 1.0 m
      └─ y: 1.0 m"""
    end
  end

  @testset "conversions" begin
    Q = typeof(T(1) * u"m")
    @testset "Cartesian <-> Polar" begin
      c1 = Cartesian(T(1), T(1))
      c2 = convert(Polar, c1)
      @test c2 ≈ Polar(T(√2), T(π / 4))
      c3 = convert(Cartesian, c2)
      @test c3 ≈ c1

      c1 = Cartesian(-T(1), T(1))
      c2 = convert(Polar, c1)
      @test c2 ≈ Polar(T(√2), T(3π / 4))
      c3 = convert(Cartesian, c2)
      @test c3 ≈ c1

      c1 = Cartesian(-T(1), -T(1))
      c2 = convert(Polar, c1)
      @test c2 ≈ Polar(T(√2), T(5π / 4))
      c3 = convert(Cartesian, c2)
      @test c3 ≈ c1

      c1 = Cartesian(T(1), -T(1))
      c2 = convert(Polar, c1)
      @test c2 ≈ Polar(T(√2), T(7π / 4))
      c3 = convert(Cartesian, c2)
      @test c3 ≈ c1

      c1 = Cartesian(T(0), T(1))
      c2 = convert(Polar, c1)
      @test c2 ≈ Polar(T(1), T(π / 2))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(Q))

      c1 = Cartesian(T(0), -T(1))
      c2 = convert(Polar, c1)
      @test c2 ≈ Polar(T(1), T(3π / 2))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(Q))

      c1 = Cartesian(T(1), T(0))
      c2 = convert(Polar, c1)
      @test c2 ≈ Polar(T(1), T(0))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(Q))

      c1 = Cartesian(-T(1), T(0))
      c2 = convert(Polar, c1)
      @test c2 ≈ Polar(T(1), T(π))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(Q))
    end

    @testset "Cartesian <-> Cylindrical" begin
      c1 = Cartesian(T(1), T(1), T(1))
      c2 = convert(Cylindrical, c1)
      @test c2 ≈ Cylindrical(T(√2), T(π / 4), T(1))
      c3 = convert(Cartesian, c2)
      @test c3 ≈ c1

      c1 = Cartesian(-T(1), T(1), T(1))
      c2 = convert(Cylindrical, c1)
      @test c2 ≈ Cylindrical(T(√2), T(3π / 4), T(1))
      c3 = convert(Cartesian, c2)
      @test c3 ≈ c1

      c1 = Cartesian(-T(1), -T(1), T(1))
      c2 = convert(Cylindrical, c1)
      @test c2 ≈ Cylindrical(T(√2), T(5π / 4), T(1))
      c3 = convert(Cartesian, c2)
      @test c3 ≈ c1

      c1 = Cartesian(T(1), -T(1), T(1))
      c2 = convert(Cylindrical, c1)
      @test c2 ≈ Cylindrical(T(√2), T(7π / 4), T(1))
      c3 = convert(Cartesian, c2)
      @test c3 ≈ c1

      c1 = Cartesian(T(0), T(1), T(1))
      c2 = convert(Cylindrical, c1)
      @test c2 ≈ Cylindrical(T(1), T(π / 2), T(1))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(Q))

      c1 = Cartesian(T(0), -T(1), T(1))
      c2 = convert(Cylindrical, c1)
      @test c2 ≈ Cylindrical(T(1), T(3π / 2), T(1))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(Q))

      c1 = Cartesian(T(1), T(0), T(1))
      c2 = convert(Cylindrical, c1)
      @test c2 ≈ Cylindrical(T(1), T(0), T(1))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(Q))

      c1 = Cartesian(-T(1), T(0), T(1))
      c2 = convert(Cylindrical, c1)
      @test c2 ≈ Cylindrical(T(1), T(π), T(1))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(Q))
    end

    @testset "Cartesian <-> Spherical" begin
      c1 = Cartesian(T(1), T(1), T(1))
      c2 = convert(Spherical, c1)
      @test c2 ≈ Spherical(T(√3), atan(T(√2)), T(π / 4))
      c3 = convert(Cartesian, c2)
      @test c3 ≈ c1

      c1 = Cartesian(-T(1), T(1), T(1))
      c2 = convert(Spherical, c1)
      @test c2 ≈ Spherical(T(√3), atan(T(√2)), T(3π / 4))
      c3 = convert(Cartesian, c2)
      @test c3 ≈ c1

      c1 = Cartesian(-T(1), -T(1), T(1))
      c2 = convert(Spherical, c1)
      @test c2 ≈ Spherical(T(√3), atan(T(√2)), T(5π / 4))
      c3 = convert(Cartesian, c2)
      @test c3 ≈ c1

      c1 = Cartesian(T(1), -T(1), T(1))
      c2 = convert(Spherical, c1)
      @test c2 ≈ Spherical(T(√3), atan(T(√2)), T(7π / 4))
      c3 = convert(Cartesian, c2)
      @test c3 ≈ c1

      c1 = Cartesian(T(0), T(1), T(1))
      c2 = convert(Spherical, c1)
      @test c2 ≈ Spherical(T(√2), T(π / 4), T(π / 2))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(Q))

      c1 = Cartesian(T(0), -T(1), T(1))
      c2 = convert(Spherical, c1)
      @test c2 ≈ Spherical(T(√2), T(π / 4), T(3π / 2))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(Q))

      c1 = Cartesian(T(1), T(0), T(1))
      c2 = convert(Spherical, c1)
      @test c2 ≈ Spherical(T(√2), T(π / 4), T(0))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(Q))

      c1 = Cartesian(-T(1), T(0), T(1))
      c2 = convert(Spherical, c1)
      @test c2 ≈ Spherical(T(√2), T(π / 4), T(π))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(Q))
    end
  end
end
