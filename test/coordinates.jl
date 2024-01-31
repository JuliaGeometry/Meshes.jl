@testset "Coordinates" begin
  @testset "Cartesian" begin
    @test Cartesian(T(1)) == Cartesian((T(1),))
    @test Cartesian(T(1), T(1)) == Cartesian((T(1), T(1)))
    @test Cartesian(T(1), T(1), T(1)) == Cartesian((T(1), T(1), T(1)))
    @test Cartesian(T(1), 1) == Cartesian((T(1), T(1)))

    c = Cartesian(T(1), T(1))
    @test sprint(show, c) == "Cartesian(coords: (1.0, 1.0))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Cartesian coordinates
      └─ coords: (1.0f0, 1.0f0)"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Cartesian coordinates
      └─ coords: (1.0, 1.0)"""
    end
  end

  @testset "Polar" begin
    @test Polar(T(1), 1) == Polar(T(1), 1.0)
    @test Polar(T(1) * u"m", 1 * u"°") == Polar(T(1) * u"m", 1.0u"°")

    c = Polar(T(1), T(1))
    @test sprint(show, c) == "Polar(ρ: 1.0, ϕ: 1.0)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Polar coordinates
      ├─ ρ: 1.0f0
      └─ ϕ: 1.0f0"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Polar coordinates
      ├─ ρ: 1.0
      └─ ϕ: 1.0"""
    end
  end

  @testset "Cylindrical" begin
    @test Cylindrical(T(1), 1, 1) == Cylindrical(T(1), 1.0, T(1))
    @test Cylindrical(T(1), T(1), 1) == Cylindrical(T(1), T(1), T(1))
    @test Cylindrical(T(1) * u"m", T(1) * u"°", 1 * u"m") == Cylindrical(T(1) * u"m", T(1) * u"°", T(1) * u"m")

    c = Cylindrical(T(1), T(1), T(1))
    @test sprint(show, c) == "Cylindrical(ρ: 1.0, ϕ: 1.0, z: 1.0)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Cylindrical coordinates
      ├─ ρ: 1.0f0
      ├─ ϕ: 1.0f0
      └─ z: 1.0f0"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Cylindrical coordinates
      ├─ ρ: 1.0
      ├─ ϕ: 1.0
      └─ z: 1.0"""
    end
  end

  @testset "Spherical" begin
    @test Spherical(T(1), 1, 1) == Spherical(T(1), 1.0, 1.0)
    @test Spherical(T(1), T(1), 1) == Spherical(T(1), T(1), T(1))
    @test Spherical(T(1) * u"m", T(1) * u"°", 1 * u"°") == Spherical(T(1) * u"m", T(1) * u"°", T(1) * u"°")

    c = Spherical(T(1), T(1), T(1))
    @test sprint(show, c) == "Spherical(r: 1.0, θ: 1.0, ϕ: 1.0)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Spherical coordinates
      ├─ r: 1.0f0
      ├─ θ: 1.0f0
      └─ ϕ: 1.0f0"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Spherical coordinates
      ├─ r: 1.0
      ├─ θ: 1.0
      └─ ϕ: 1.0"""
    end
  end

  @testset "LatLon" begin
    @test LatLon(T(1) * u"°", 1 * u"°") == LatLon(T(1) * u"°", T(1) * u"°")
    @test LatLon(T(1), 1) == LatLon(T(1) * u"°", T(1) * u"°")

    # error: the units of "lat" and "lon" must be degrees
    @test_throws ArgumentError LatLon(T(1) * u"s", T(1) * u"s")

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
    @test LatLonAlt(T(1) * u"°", 1 * u"°", T(1) * u"km") == LatLonAlt(T(1) * u"°", T(1) * u"°", T(1) * u"km")
    @test LatLonAlt(T(1), 1, T(1)) == LatLonAlt(T(1) * u"°", T(1) * u"°", T(1) * u"m")

    # error: the units of "lat" and "lon" must be degrees
    @test_throws ArgumentError LatLonAlt(T(1) * u"s", T(1) * u"s", T(1) * u"km")
    # error: the unit of "alt" must be a length unit
    @test_throws ArgumentError LatLonAlt(T(1) * u"°", T(1) * u"°", T(1) * u"s")

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
    @test EastNorth(T(1) * u"km", 1 * u"km") == EastNorth(T(1) * u"km", T(1) * u"km")
    @test EastNorth(T(1), 1) == EastNorth(T(1) * u"m", T(1) * u"m")

    # error: the units of "east" and "north" must be length units
    @test_throws ArgumentError EastNorth(T(1) * u"s", T(1) * u"s")

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
    @test WebMercator(T(1) * u"km", 1 * u"km") == WebMercator(T(1) * u"km", T(1) * u"km")
    @test WebMercator(T(1), 1) == WebMercator(T(1) * u"m", T(1) * u"m")

    # error: the units of "x" and "y" must be length units
    @test_throws ArgumentError WebMercator(T(1) * u"s", T(1) * u"s")

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
    @testset "Cartesian <--> Polar" begin
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
      @test isapprox(c3, c1, atol=atol(T))

      c1 = Cartesian(T(0), -T(1))
      c2 = convert(Polar, c1)
      @test c2 ≈ Polar(T(1), T(3π / 2))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(T))

      c1 = Cartesian(T(1), T(0))
      c2 = convert(Polar, c1)
      @test c2 ≈ Polar(T(1), T(0))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(T))

      c1 = Cartesian(-T(1), T(0))
      c2 = convert(Polar, c1)
      @test c2 ≈ Polar(T(1), T(π))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(T))
    end

    @testset "Cartesian <--> Cylindrical" begin
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
      @test isapprox(c3, c1, atol=atol(T))

      c1 = Cartesian(T(0), -T(1), T(1))
      c2 = convert(Cylindrical, c1)
      @test c2 ≈ Cylindrical(T(1), T(3π / 2), T(1))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(T))

      c1 = Cartesian(T(1), T(0), T(1))
      c2 = convert(Cylindrical, c1)
      @test c2 ≈ Cylindrical(T(1), T(0), T(1))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(T))

      c1 = Cartesian(-T(1), T(0), T(1))
      c2 = convert(Cylindrical, c1)
      @test c2 ≈ Cylindrical(T(1), T(π), T(1))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(T))
    end

    @testset "Cartesian <--> Spherical" begin
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
      @test isapprox(c3, c1, atol=atol(T))

      c1 = Cartesian(T(0), -T(1), T(1))
      c2 = convert(Spherical, c1)
      @test c2 ≈ Spherical(T(√2), T(π / 4), T(3π / 2))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(T))

      c1 = Cartesian(T(1), T(0), T(1))
      c2 = convert(Spherical, c1)
      @test c2 ≈ Spherical(T(√2), T(π / 4), T(0))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(T))

      c1 = Cartesian(-T(1), T(0), T(1))
      c2 = convert(Spherical, c1)
      @test c2 ≈ Spherical(T(√2), T(π / 4), T(π))
      c3 = convert(Cartesian, c2)
      @test isapprox(c3, c1, atol=atol(T))
    end
  end
end
