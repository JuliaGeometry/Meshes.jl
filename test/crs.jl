@testset "CRS" begin
  @testset "Datum" begin
    c = Cartesian(T(1), T(1))
    @test datum(c) === NoDatum
    @test isnothing(latitudeₒ(c))
    @test isnothing(longitudeₒ(c))
    @test isnothing(altitudeₒ(c))

    c = LatLon(T(1), T(1))
    @test datum(c) === WGS84
    @test latitudeₒ(c) == latitudeₒ(WGS84)
    @test longitudeₒ(c) == longitudeₒ(WGS84)
    @test altitudeₒ(c) == altitudeₒ(WGS84)
  end

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

    # error: invalid units for coordinates
    @test_throws ArgumentError Cartesian(T(1), T(1) * u"m")
    @test_throws ArgumentError Cartesian(T(1) * u"s", T(1) * u"m")
    @test_throws ArgumentError Cartesian(T(1) * u"m", T(1) * u"s")
    @test_throws ArgumentError Cartesian(T(1) * u"s", T(1) * u"s")
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

    # error: invalid units for coordinates
    @test_throws ArgumentError Polar(T(1), T(1) * u"rad")
    @test_throws ArgumentError Polar(T(1) * u"s", T(1) * u"rad")
    @test_throws ArgumentError Polar(T(1) * u"m", T(1) * u"s")
    @test_throws ArgumentError Polar(T(1) * u"s", T(1) * u"s")
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

    # error: invalid units for coordinates
    @test_throws ArgumentError Cylindrical(T(1), T(1) * u"rad", T(1))
    @test_throws ArgumentError Cylindrical(T(1) * u"s", T(1) * u"rad", T(1) * u"m")
    @test_throws ArgumentError Cylindrical(T(1) * u"m", T(1) * u"s", T(1) * u"m")
    @test_throws ArgumentError Cylindrical(T(1) * u"m", T(1) * u"rad", T(1) * u"s")
    @test_throws ArgumentError Cylindrical(T(1) * u"s", T(1) * u"s", T(1) * u"s")
  end

  @testset "Spherical" begin
    @test Spherical(T(1), T(1), T(1)) == Spherical(T(1) * u"m", T(1) * u"rad", T(1) * u"rad")
    @test Spherical(T(1) * u"m", T(1) * u"rad", 1 * u"rad") == Spherical(T(1) * u"m", T(1) * u"rad", T(1) * u"rad")
    @test Spherical(T(1) * u"m", T(45) * u"°", T(45) * u"°") ≈
          Spherical(T(1) * u"m", T(π / 4) * u"rad", T(π / 4) * u"rad")

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

    # error: invalid units for coordinates
    @test_throws ArgumentError Spherical(T(1) * u"m", T(1), T(1))
    @test_throws ArgumentError Spherical(T(1) * u"s", T(1) * u"rad", T(1) * u"rad")
    @test_throws ArgumentError Spherical(T(1) * u"m", T(1) * u"s", T(1) * u"rad")
    @test_throws ArgumentError Spherical(T(1) * u"m", T(1) * u"rad", T(1) * u"s")
    @test_throws ArgumentError Spherical(T(1) * u"s", T(1) * u"s", T(1) * u"s")
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

    # error: invalid units for coordinates
    @test_throws ArgumentError LatLon(T(1), T(1) * u"°")
    @test_throws ArgumentError LatLon(T(1) * u"s", T(1) * u"°")
    @test_throws ArgumentError LatLon(T(1) * u"°", T(1) * u"s")
    @test_throws ArgumentError LatLon(T(1) * u"s", T(1) * u"s")
  end

  @testset "Mercator" begin
    @test Mercator(T(1), T(1)) == Mercator(T(1) * u"m", T(1) * u"m")
    @test Mercator(T(1) * u"m", 1 * u"m") == Mercator(T(1) * u"m", T(1) * u"m")
    @test Mercator(T(1) * u"km", T(1) * u"km") == Mercator(T(1000) * u"m", T(1000) * u"m")

    c = Mercator(T(1), T(1))
    @test sprint(show, c) == "Mercator(x: 1.0 m, y: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Mercator coordinates
      ├─ x: 1.0f0 m
      └─ y: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Mercator coordinates
      ├─ x: 1.0 m
      └─ y: 1.0 m"""
    end

    # error: invalid units for coordinates
    @test_throws ArgumentError Mercator(T(1), T(1) * u"m")
    @test_throws ArgumentError Mercator(T(1) * u"s", T(1) * u"m")
    @test_throws ArgumentError Mercator(T(1) * u"m", T(1) * u"s")
    @test_throws ArgumentError Mercator(T(1) * u"s", T(1) * u"s")
  end

  @testset "WebMercator" begin
    @test WebMercator(T(1), T(1)) == WebMercator(T(1) * u"m", T(1) * u"m")
    @test WebMercator(T(1) * u"m", 1 * u"m") == WebMercator(T(1) * u"m", T(1) * u"m")
    @test WebMercator(T(1) * u"km", T(1) * u"km") == WebMercator(T(1000) * u"m", T(1000) * u"m")

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

    # error: invalid units for coordinates
    @test_throws ArgumentError WebMercator(T(1), T(1) * u"m")
    @test_throws ArgumentError WebMercator(T(1) * u"s", T(1) * u"m")
    @test_throws ArgumentError WebMercator(T(1) * u"m", T(1) * u"s")
    @test_throws ArgumentError WebMercator(T(1) * u"s", T(1) * u"s")
  end

  @testset "PlateCarree" begin
    @test PlateCarree(T(1), T(1)) == PlateCarree(T(1) * u"m", T(1) * u"m")
    @test PlateCarree(T(1) * u"m", 1 * u"m") == PlateCarree(T(1) * u"m", T(1) * u"m")
    @test PlateCarree(T(1) * u"km", T(1) * u"km") == PlateCarree(T(1000) * u"m", T(1000) * u"m")

    c = PlateCarree(T(1), T(1))
    @test sprint(show, c) == "PlateCarree(x: 1.0 m, y: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      PlateCarree coordinates
      ├─ x: 1.0f0 m
      └─ y: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      PlateCarree coordinates
      ├─ x: 1.0 m
      └─ y: 1.0 m"""
    end

    # error: invalid units for coordinates
    @test_throws ArgumentError PlateCarree(T(1), T(1) * u"m")
    @test_throws ArgumentError PlateCarree(T(1) * u"s", T(1) * u"m")
    @test_throws ArgumentError PlateCarree(T(1) * u"m", T(1) * u"s")
    @test_throws ArgumentError PlateCarree(T(1) * u"s", T(1) * u"s")
  end

  @testset "conversions" begin
    Q = typeof(T(1) * u"m")
    @testset "Cartesian <> Polar" begin
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

      # type stability
      c1 = Cartesian(T(1), T(1))
      c2 = Polar(T(√2), T(π / 4))
      @inferred convert(Polar, c1)
      @inferred convert(Cartesian, c2)
    end

    @testset "Cartesian <> Cylindrical" begin
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

      # type stability
      c1 = Cartesian(T(1), T(1), T(1))
      c2 = Cylindrical(T(√2), T(π / 4), T(1))
      @inferred convert(Cylindrical, c1)
      @inferred convert(Cartesian, c2)
    end

    @testset "Cartesian <> Spherical" begin
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

      # type stability
      c1 = Cartesian(T(1), T(1), T(1))
      c2 = Spherical(T(√3), atan(T(√2)), T(π / 4))
      @inferred convert(Spherical, c1)
      @inferred convert(Cartesian, c2)
    end

    @testset "LatLon <> Mercator" begin
      c1 = LatLon(T(45), T(90))
      c2 = convert(Mercator, c1)
      @test c2 ≈ Mercator(T(10018754.171394622), T(5591295.9185533915))

      c1 = LatLon(-T(45), T(90))
      c2 = convert(Mercator, c1)
      @test c2 ≈ Mercator(T(10018754.171394622), -T(5591295.9185533915))

      c1 = LatLon(T(45), -T(90))
      c2 = convert(Mercator, c1)
      @test c2 ≈ Mercator(-T(10018754.171394622), T(5591295.9185533915))

      c1 = LatLon(-T(45), -T(90))
      c2 = convert(Mercator, c1)
      @test c2 ≈ Mercator(-T(10018754.171394622), -T(5591295.9185533915))

      # EPSG fallback
      c1 = LatLon(T(45), T(90))
      c2 = convert(EPSG{3395}, c1)
      @test c2 ≈ Mercator(T(10018754.171394622), T(5591295.9185533915))

      # type stability
      c1 = LatLon(T(45), T(90))
      @inferred convert(Mercator, c1)
      @inferred convert(EPSG{3395}, c1)
    end

    @testset "LatLon <> WebMercator" begin
      c1 = LatLon(T(45), T(90))
      c2 = convert(WebMercator, c1)
      @test c2 ≈ WebMercator(T(10018754.171394622), T(5621521.486192066))
      c3 = convert(LatLon, c2)
      @test c3 ≈ c1

      c1 = LatLon(-T(45), T(90))
      c2 = convert(WebMercator, c1)
      @test c2 ≈ WebMercator(T(10018754.171394622), -T(5621521.486192066))
      c3 = convert(LatLon, c2)
      @test c3 ≈ c1

      c1 = LatLon(T(45), -T(90))
      c2 = convert(WebMercator, c1)
      @test c2 ≈ WebMercator(-T(10018754.171394622), T(5621521.486192066))
      c3 = convert(LatLon, c2)
      @test c3 ≈ c1

      c1 = LatLon(-T(45), -T(90))
      c2 = convert(WebMercator, c1)
      @test c2 ≈ WebMercator(-T(10018754.171394622), -T(5621521.486192066))
      c3 = convert(LatLon, c2)
      @test c3 ≈ c1

      # EPSG fallback
      c1 = LatLon(T(45), T(90))
      c2 = convert(EPSG{3857}, c1)
      @test c2 ≈ WebMercator(T(10018754.171394622), T(5621521.486192066))
      c3 = convert(EPSG{4326}, c2)
      @test c3 ≈ c1

      # type stability
      c1 = LatLon(T(45), T(90))
      c2 = WebMercator(T(10018754.171394622), T(5621521.486192066))
      @inferred convert(WebMercator, c1)
      @inferred convert(LatLon, c2)
      @inferred convert(EPSG{3857}, c1)
      @inferred convert(EPSG{4326}, c2)
    end

    @testset "LatLon <> PlateCarree" begin
      c1 = LatLon(T(45), T(90))
      c2 = convert(PlateCarree, c1)
      @test c2 ≈ PlateCarree(T(10018754.171394622), T(5009377.085697311))
      c3 = convert(LatLon, c2)
      @test c3 ≈ c1

      c1 = LatLon(-T(45), T(90))
      c2 = convert(PlateCarree, c1)
      @test c2 ≈ PlateCarree(T(10018754.171394622), -T(5009377.085697311))
      c3 = convert(LatLon, c2)
      @test c3 ≈ c1

      c1 = LatLon(T(45), -T(90))
      c2 = convert(PlateCarree, c1)
      @test c2 ≈ PlateCarree(-T(10018754.171394622), T(5009377.085697311))
      c3 = convert(LatLon, c2)
      @test c3 ≈ c1

      c1 = LatLon(-T(45), -T(90))
      c2 = convert(PlateCarree, c1)
      @test c2 ≈ PlateCarree(-T(10018754.171394622), -T(5009377.085697311))
      c3 = convert(LatLon, c2)
      @test c3 ≈ c1

      # EPSG fallback
      c1 = LatLon(T(45), T(90))
      c2 = convert(EPSG{32662}, c1)
      @test c2 ≈ PlateCarree(T(10018754.171394622), T(5009377.085697311))
      c3 = convert(EPSG{4326}, c2)
      @test c3 ≈ c1

      # type stability
      c1 = LatLon(T(45), T(90))
      c2 = PlateCarree(T(10018754.171394622), T(5009377.085697311))
      @inferred convert(PlateCarree, c1)
      @inferred convert(LatLon, c2)
      @inferred convert(EPSG{32662}, c1)
      @inferred convert(EPSG{4326}, c2)
    end
  end
end
