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

  @testset "Lambert" begin
    @test Lambert(T(1), T(1)) == Lambert(T(1) * u"m", T(1) * u"m")
    @test Lambert(T(1) * u"m", 1 * u"m") == Lambert(T(1) * u"m", T(1) * u"m")
    @test Lambert(T(1) * u"km", T(1) * u"km") == Lambert(T(1000) * u"m", T(1000) * u"m")

    c = Lambert(T(1), T(1))
    @test sprint(show, c) == "Lambert(x: 1.0 m, y: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Lambert coordinates
      ├─ x: 1.0f0 m
      └─ y: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Lambert coordinates
      ├─ x: 1.0 m
      └─ y: 1.0 m"""
    end

    # error: invalid units for coordinates
    @test_throws ArgumentError Lambert(T(1), T(1) * u"m")
    @test_throws ArgumentError Lambert(T(1) * u"s", T(1) * u"m")
    @test_throws ArgumentError Lambert(T(1) * u"m", T(1) * u"s")
    @test_throws ArgumentError Lambert(T(1) * u"s", T(1) * u"s")
  end

  @testset "Behrmann" begin
    @test Behrmann(T(1), T(1)) == Behrmann(T(1) * u"m", T(1) * u"m")
    @test Behrmann(T(1) * u"m", 1 * u"m") == Behrmann(T(1) * u"m", T(1) * u"m")
    @test Behrmann(T(1) * u"km", T(1) * u"km") == Behrmann(T(1000) * u"m", T(1000) * u"m")

    c = Behrmann(T(1), T(1))
    @test sprint(show, c) == "Behrmann(x: 1.0 m, y: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Behrmann coordinates
      ├─ x: 1.0f0 m
      └─ y: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Behrmann coordinates
      ├─ x: 1.0 m
      └─ y: 1.0 m"""
    end

    # error: invalid units for coordinates
    @test_throws ArgumentError Behrmann(T(1), T(1) * u"m")
    @test_throws ArgumentError Behrmann(T(1) * u"s", T(1) * u"m")
    @test_throws ArgumentError Behrmann(T(1) * u"m", T(1) * u"s")
    @test_throws ArgumentError Behrmann(T(1) * u"s", T(1) * u"s")
  end

  @testset "GallPeters" begin
    @test GallPeters(T(1), T(1)) == GallPeters(T(1) * u"m", T(1) * u"m")
    @test GallPeters(T(1) * u"m", 1 * u"m") == GallPeters(T(1) * u"m", T(1) * u"m")
    @test GallPeters(T(1) * u"km", T(1) * u"km") == GallPeters(T(1000) * u"m", T(1000) * u"m")

    c = GallPeters(T(1), T(1))
    @test sprint(show, c) == "GallPeters(x: 1.0 m, y: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      GallPeters coordinates
      ├─ x: 1.0f0 m
      └─ y: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      GallPeters coordinates
      ├─ x: 1.0 m
      └─ y: 1.0 m"""
    end

    # error: invalid units for coordinates
    @test_throws ArgumentError GallPeters(T(1), T(1) * u"m")
    @test_throws ArgumentError GallPeters(T(1) * u"s", T(1) * u"m")
    @test_throws ArgumentError GallPeters(T(1) * u"m", T(1) * u"s")
    @test_throws ArgumentError GallPeters(T(1) * u"s", T(1) * u"s")
  end

  @testset "WinkelTripel" begin
    @test WinkelTripel(T(1), T(1)) == WinkelTripel(T(1) * u"m", T(1) * u"m")
    @test WinkelTripel(T(1) * u"m", 1 * u"m") == WinkelTripel(T(1) * u"m", T(1) * u"m")
    @test WinkelTripel(T(1) * u"km", T(1) * u"km") == WinkelTripel(T(1000) * u"m", T(1000) * u"m")

    c = WinkelTripel(T(1), T(1))
    @test sprint(show, c) == "WinkelTripel(x: 1.0 m, y: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      WinkelTripel coordinates
      ├─ x: 1.0f0 m
      └─ y: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      WinkelTripel coordinates
      ├─ x: 1.0 m
      └─ y: 1.0 m"""
    end

    # error: invalid units for coordinates
    @test_throws ArgumentError WinkelTripel(T(1), T(1) * u"m")
    @test_throws ArgumentError WinkelTripel(T(1) * u"s", T(1) * u"m")
    @test_throws ArgumentError WinkelTripel(T(1) * u"m", T(1) * u"s")
    @test_throws ArgumentError WinkelTripel(T(1) * u"s", T(1) * u"s")
  end

  @testset "Robinson" begin
    @test Robinson(T(1), T(1)) == Robinson(T(1) * u"m", T(1) * u"m")
    @test Robinson(T(1) * u"m", 1 * u"m") == Robinson(T(1) * u"m", T(1) * u"m")
    @test Robinson(T(1) * u"km", T(1) * u"km") == Robinson(T(1000) * u"m", T(1000) * u"m")

    c = Robinson(T(1), T(1))
    @test sprint(show, c) == "Robinson(x: 1.0 m, y: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      Robinson coordinates
      ├─ x: 1.0f0 m
      └─ y: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      Robinson coordinates
      ├─ x: 1.0 m
      └─ y: 1.0 m"""
    end

    # error: invalid units for coordinates
    @test_throws ArgumentError Robinson(T(1), T(1) * u"m")
    @test_throws ArgumentError Robinson(T(1) * u"s", T(1) * u"m")
    @test_throws ArgumentError Robinson(T(1) * u"m", T(1) * u"s")
    @test_throws ArgumentError Robinson(T(1) * u"s", T(1) * u"s")
  end

  @testset "OrthoNorth" begin
    @test OrthoNorth(T(1), T(1)) == OrthoNorth(T(1) * u"m", T(1) * u"m")
    @test OrthoNorth(T(1) * u"m", 1 * u"m") == OrthoNorth(T(1) * u"m", T(1) * u"m")
    @test OrthoNorth(T(1) * u"km", T(1) * u"km") == OrthoNorth(T(1000) * u"m", T(1000) * u"m")

    c = OrthoNorth(T(1), T(1))
    @test sprint(show, c) == "OrthoNorth(x: 1.0 m, y: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      OrthoNorth coordinates
      ├─ x: 1.0f0 m
      └─ y: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      OrthoNorth coordinates
      ├─ x: 1.0 m
      └─ y: 1.0 m"""
    end

    # error: invalid units for coordinates
    @test_throws ArgumentError OrthoNorth(T(1), T(1) * u"m")
    @test_throws ArgumentError OrthoNorth(T(1) * u"s", T(1) * u"m")
    @test_throws ArgumentError OrthoNorth(T(1) * u"m", T(1) * u"s")
    @test_throws ArgumentError OrthoNorth(T(1) * u"s", T(1) * u"s")
  end

  @testset "OrthoSouth" begin
    @test OrthoSouth(T(1), T(1)) == OrthoSouth(T(1) * u"m", T(1) * u"m")
    @test OrthoSouth(T(1) * u"m", 1 * u"m") == OrthoSouth(T(1) * u"m", T(1) * u"m")
    @test OrthoSouth(T(1) * u"km", T(1) * u"km") == OrthoSouth(T(1000) * u"m", T(1000) * u"m")

    c = OrthoSouth(T(1), T(1))
    @test sprint(show, c) == "OrthoSouth(x: 1.0 m, y: 1.0 m)"
    if T === Float32
      @test sprint(show, MIME("text/plain"), c) == """
      OrthoSouth coordinates
      ├─ x: 1.0f0 m
      └─ y: 1.0f0 m"""
    else
      @test sprint(show, MIME("text/plain"), c) == """
      OrthoSouth coordinates
      ├─ x: 1.0 m
      └─ y: 1.0 m"""
    end

    # error: invalid units for coordinates
    @test_throws ArgumentError OrthoSouth(T(1), T(1) * u"m")
    @test_throws ArgumentError OrthoSouth(T(1) * u"s", T(1) * u"m")
    @test_throws ArgumentError OrthoSouth(T(1) * u"m", T(1) * u"s")
    @test_throws ArgumentError OrthoSouth(T(1) * u"s", T(1) * u"s")
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
      c2 = convert(Mercator{WGS84}, c1)
      @test c2 ≈ Mercator(T(10018754.171394622), T(5591295.9185533915))

      c1 = LatLon(-T(45), T(90))
      c2 = convert(Mercator{WGS84}, c1)
      @test c2 ≈ Mercator(T(10018754.171394622), -T(5591295.9185533915))

      c1 = LatLon(T(45), -T(90))
      c2 = convert(Mercator{WGS84}, c1)
      @test c2 ≈ Mercator(-T(10018754.171394622), T(5591295.9185533915))

      c1 = LatLon(-T(45), -T(90))
      c2 = convert(Mercator{WGS84}, c1)
      @test c2 ≈ Mercator(-T(10018754.171394622), -T(5591295.9185533915))

      # type stability
      c1 = LatLon(T(45), T(90))
      @inferred convert(Mercator{WGS84}, c1)
    end

    @testset "LatLon <> WebMercator" begin
      c1 = LatLon(T(45), T(90))
      c2 = convert(WebMercator{WGS84}, c1)
      @test c2 ≈ WebMercator(T(10018754.171394622), T(5621521.486192066))
      c3 = convert(LatLon{WGS84}, c2)
      @test c3 ≈ c1

      c1 = LatLon(-T(45), T(90))
      c2 = convert(WebMercator{WGS84}, c1)
      @test c2 ≈ WebMercator(T(10018754.171394622), -T(5621521.486192066))
      c3 = convert(LatLon{WGS84}, c2)
      @test c3 ≈ c1

      c1 = LatLon(T(45), -T(90))
      c2 = convert(WebMercator{WGS84}, c1)
      @test c2 ≈ WebMercator(-T(10018754.171394622), T(5621521.486192066))
      c3 = convert(LatLon{WGS84}, c2)
      @test c3 ≈ c1

      c1 = LatLon(-T(45), -T(90))
      c2 = convert(WebMercator{WGS84}, c1)
      @test c2 ≈ WebMercator(-T(10018754.171394622), -T(5621521.486192066))
      c3 = convert(LatLon{WGS84}, c2)
      @test c3 ≈ c1

      # type stability
      c1 = LatLon(T(45), T(90))
      c2 = WebMercator(T(10018754.171394622), T(5621521.486192066))
      @inferred convert(WebMercator{WGS84}, c1)
      @inferred convert(LatLon{WGS84}, c2)
    end

    @testset "LatLon <> PlateCarree" begin
      c1 = LatLon(T(45), T(90))
      c2 = convert(PlateCarree{WGS84}, c1)
      @test c2 ≈ PlateCarree(T(10018754.171394622), T(5009377.085697311))
      c3 = convert(LatLon{WGS84}, c2)
      @test c3 ≈ c1

      c1 = LatLon(-T(45), T(90))
      c2 = convert(PlateCarree{WGS84}, c1)
      @test c2 ≈ PlateCarree(T(10018754.171394622), -T(5009377.085697311))
      c3 = convert(LatLon{WGS84}, c2)
      @test c3 ≈ c1

      c1 = LatLon(T(45), -T(90))
      c2 = convert(PlateCarree{WGS84}, c1)
      @test c2 ≈ PlateCarree(-T(10018754.171394622), T(5009377.085697311))
      c3 = convert(LatLon{WGS84}, c2)
      @test c3 ≈ c1

      c1 = LatLon(-T(45), -T(90))
      c2 = convert(PlateCarree{WGS84}, c1)
      @test c2 ≈ PlateCarree(-T(10018754.171394622), -T(5009377.085697311))
      c3 = convert(LatLon{WGS84}, c2)
      @test c3 ≈ c1

      # type stability
      c1 = LatLon(T(45), T(90))
      c2 = PlateCarree(T(10018754.171394622), T(5009377.085697311))
      @inferred convert(PlateCarree{WGS84}, c1)
      @inferred convert(LatLon{WGS84}, c2)
    end

    @testset "LatLon <> Lambert" begin
      c1 = LatLon(T(45), T(90))
      c2 = convert(Lambert{WGS84}, c1)
      @test c2 ≈ Lambert(T(10018754.171394622), T(4489858.8869480025))

      c1 = LatLon(-T(45), T(90))
      c2 = convert(Lambert{WGS84}, c1)
      @test c2 ≈ Lambert(T(10018754.171394622), -T(4489858.8869480025))

      c1 = LatLon(T(45), -T(90))
      c2 = convert(Lambert{WGS84}, c1)
      @test c2 ≈ Lambert(-T(10018754.171394622), T(4489858.8869480025))

      c1 = LatLon(-T(45), -T(90))
      c2 = convert(Lambert{WGS84}, c1)
      @test c2 ≈ Lambert(-T(10018754.171394622), -T(4489858.8869480025))

      # type stability
      c1 = LatLon(T(45), T(90))
      @inferred convert(Lambert{WGS84}, c1)
    end

    @testset "LatLon <> Behrmann" begin
      c1 = LatLon(T(45), T(90))
      c2 = convert(Behrmann{WGS84}, c1)
      @test c2 ≈ Behrmann(T(8683765.222580686), T(5180102.328839251))

      c1 = LatLon(-T(45), T(90))
      c2 = convert(Behrmann{WGS84}, c1)
      @test c2 ≈ Behrmann(T(8683765.222580686), -T(5180102.328839251))

      c1 = LatLon(T(45), -T(90))
      c2 = convert(Behrmann{WGS84}, c1)
      @test c2 ≈ Behrmann(-T(8683765.222580686), T(5180102.328839251))

      c1 = LatLon(-T(45), -T(90))
      c2 = convert(Behrmann{WGS84}, c1)
      @test c2 ≈ Behrmann(-T(8683765.222580686), -T(5180102.328839251))

      # type stability
      c1 = LatLon(T(45), T(90))
      @inferred convert(Behrmann{WGS84}, c1)
    end

    @testset "LatLon <> GallPeters" begin
      c1 = LatLon(T(45), T(90))
      c2 = convert(GallPeters{WGS84}, c1)
      @test c2 ≈ GallPeters(T(7096215.158458031), T(6338983.732612475))

      c1 = LatLon(-T(45), T(90))
      c2 = convert(GallPeters{WGS84}, c1)
      @test c2 ≈ GallPeters(T(7096215.158458031), -T(6338983.732612475))

      c1 = LatLon(T(45), -T(90))
      c2 = convert(GallPeters{WGS84}, c1)
      @test c2 ≈ GallPeters(-T(7096215.158458031), T(6338983.732612475))

      c1 = LatLon(-T(45), -T(90))
      c2 = convert(GallPeters{WGS84}, c1)
      @test c2 ≈ GallPeters(-T(7096215.158458031), -T(6338983.732612475))

      # type stability
      c1 = LatLon(T(45), T(90))
      @inferred convert(GallPeters{WGS84}, c1)
    end

    @testset "LatLon <> WinkelTripel" begin
      c1 = LatLon(T(45), T(90))
      c2 = convert(WinkelTripel{WGS84}, c1)
      @test c2 ≈ WinkelTripel(T(7044801.6979576545), T(5231448.051548355))

      c1 = LatLon(-T(45), T(90))
      c2 = convert(WinkelTripel{WGS84}, c1)
      @test c2 ≈ WinkelTripel(T(7044801.6979576545), -T(5231448.051548355))

      c1 = LatLon(T(45), -T(90))
      c2 = convert(WinkelTripel{WGS84}, c1)
      @test c2 ≈ WinkelTripel(-T(7044801.6979576545), T(5231448.051548355))

      c1 = LatLon(-T(45), -T(90))
      c2 = convert(WinkelTripel{WGS84}, c1)
      @test c2 ≈ WinkelTripel(-T(7044801.6979576545), -T(5231448.051548355))

      c1 = LatLon(T(0), T(0))
      c2 = convert(WinkelTripel{WGS84}, c1)
      @test c2 ≈ WinkelTripel(T(0), T(0))

      # type stability
      c1 = LatLon(T(45), T(90))
      @inferred convert(WinkelTripel{WGS84}, c1)
    end

    @testset "LatLon <> Robinson" begin
      c1 = LatLon(T(45), T(90))
      c2 = convert(Robinson{WGS84}, c1)
      @test c2 ≈ Robinson(T(7620313.925950073), T(4805073.646653474))

      c1 = LatLon(-T(45), T(90))
      c2 = convert(Robinson{WGS84}, c1)
      @test c2 ≈ Robinson(T(7620313.925950073), -T(4805073.646653474))

      c1 = LatLon(T(45), -T(90))
      c2 = convert(Robinson{WGS84}, c1)
      @test c2 ≈ Robinson(-T(7620313.925950073), T(4805073.646653474))

      c1 = LatLon(-T(45), -T(90))
      c2 = convert(Robinson{WGS84}, c1)
      @test c2 ≈ Robinson(-T(7620313.925950073), -T(4805073.646653474))

      # type stability
      c1 = LatLon(T(45), T(90))
      @inferred convert(Robinson{WGS84}, c1)
    end

    @testset "LatLon <> OrthoNorth" begin
      c1 = LatLon(T(30), T(60))
      c2 = convert(OrthoNorth{WGS84}, c1)
      @test c2 ≈ OrthoNorth(T(4787610.688267582), T(-2764128.319646418))

      c1 = LatLon(T(30), -T(60))
      c2 = convert(OrthoNorth{WGS84}, c1)
      @test c2 ≈ OrthoNorth(-T(4787610.688267582), T(-2764128.319646418))

      # type stability
      c1 = LatLon(T(30), T(60))
      @inferred convert(OrthoNorth{WGS84}, c1)
    end

    @testset "LatLon <> OrthoSouth" begin
      c1 = LatLon(-T(30), T(60))
      c2 = convert(OrthoSouth{WGS84}, c1)
      @test c2 ≈ OrthoSouth(T(4787610.688267582), T(2764128.319646418))

      c1 = LatLon(-T(30), -T(60))
      c2 = convert(OrthoSouth{WGS84}, c1)
      @test c2 ≈ OrthoSouth(-T(4787610.688267582), T(2764128.319646418))

      # type stability
      c1 = LatLon(T(30), T(60))
      @inferred convert(OrthoSouth{WGS84}, c1)
    end

    @testset "LatLon <> OrthoSpherical" begin
      OrthoNorthSpherical = Meshes.Orthographic{90.0u"°",0.0u"°",true,WGS84}
      OrthoSouthSpherical = Meshes.Orthographic{-90.0u"°",0.0u"°",true,WGS84} 

      c1 = LatLon(T(30), T(60))
      c2 = convert(OrthoNorthSpherical, c1)
      @test c2 ≈ OrthoNorthSpherical(T(4783602.75), T(-2761814.335408735))

      c1 = LatLon(T(30), -T(60))
      c2 = convert(OrthoNorthSpherical, c1)
      @test c2 ≈ OrthoNorthSpherical(-T(4783602.75), T(-2761814.335408735))

      c1 = LatLon(-T(30), T(60))
      c2 = convert(OrthoSouthSpherical, c1)
      @test c2 ≈ OrthoSouthSpherical(T(4783602.75), T(2761814.335408735))

      c1 = LatLon(-T(30), -T(60))
      c2 = convert(OrthoSouthSpherical, c1)
      @test c2 ≈ OrthoSouthSpherical(-T(4783602.75), T(2761814.335408735))

      # type stability
      c1 = LatLon(T(30), T(60))
      @inferred convert(OrthoNorthSpherical, c1)
      @inferred convert(OrthoSouthSpherical, c1)
    end
  end
end
