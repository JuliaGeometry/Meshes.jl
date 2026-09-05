@testitem "GeometricDistance" setup = [Setup] begin
  # ----------------
  # EUCLIDEAN SPACE
  # ----------------

  d = EuclideanDistance()
  @test d(cart(0, 0), cart(3, 4)) == T(5) * u"m"
  @test d(cart(0, 0, 0), cart(1, 2, 2)) == T(3) * u"m"
  @test d(cart(1, 1), cart(1, 1)) == T(0) * u"m"

  # the geodesic of Euclidean space is the straight line
  d = GeodesicDistance()
  @test d(cart(0, 0), cart(3, 4)) == EuclideanDistance()(cart(0, 0), cart(3, 4))
  @test d(cart(0, 0, 0), cart(1, 2, 2)) == EuclideanDistance()(cart(0, 0, 0), cart(1, 2, 2))

  # ----------
  # ELLIPSOID
  # ----------

  d = GeodesicDistance()
  e = EuclideanDistance()

  # the Euclidean distance is the chord, which is shorter than the geodesic
  p₁ = latlon(-33.8688, 151.2093)
  p₂ = latlon(51.5074, -0.1278)
  @test e(p₁, p₂) < d(p₁, p₂)

  # identity and symmetry
  @test d(p₁, p₁) == T(0) * u"m"
  @test d(p₁, p₂) == d(p₂, p₁)

  # units are inherited from the ellipsoid of the datum
  @test unit(d(p₁, p₂)) == u"m"
  @test Unitful.numtype(d(p₁, p₂)) === T

  if T === Float64
    # reference values from the test set of Karney (2013), which is
    # computed with high precision arithmetic on the WGS84 ellipsoid
    # https://geographiclib.sourceforge.io/C++/doc/geodesic.html#testgeod
    τ = 1e-6u"m"
    @test isapprox(
      d(latlon(36.530042355041, 0), latlon(-48.164270779097769, 5.7623446946765105)),
      9398502.0434687u"m",
      atol=τ
    )
    @test isapprox(
      d(latlon(31.371144087006, 0), latlon(28.909631785856762, 17.071526465331372)),
      1665550.2846815u"m",
      atol=τ
    )
    @test isapprox(
      d(latlon(19.707097385334, 0), latlon(19.707739582810257, 0.00026025642810105)),
      76.1478894u"m",
      atol=τ
    )

    # nearly antipodal points, where the series of Vincenty fail to converge
    @test isapprox(
      d(latlon(41.211099672963, 0), latlon(-41.567531864088980, 179.92621956183250)),
      19964100.0439999u"m",
      atol=τ
    )
    @test isapprox(
      d(latlon(0.925973600837, 0), latlon(-0.607398400683232, 179.95083146293813)),
      19968558.8681118u"m",
      atol=τ
    )
    @test isapprox(
      d(latlon(74.098356057014, 0), latlon(-74.073175364040975, 179.99585371619742)),
      20001120.0416458u"m",
      atol=τ
    )

    # a quarter of the equator is a quarter of the equatorial circumference
    a = CoordRefSystems.majoraxis(ellipsoid(WGS84Latest))
    @test isapprox(d(latlon(0, 0), latlon(0, 90)), π * a / 2, atol=τ)

    # pole to pole is half of the meridian, obtained by numerical
    # integration of the meridional radius of curvature
    @test isapprox(d(latlon(90, 0), latlon(-90, 0)), 20003931.4586254u"m", atol=τ)

    # on a datum with a spherical ellipsoid the geodesic is the great circle
    R = CoordRefSystems.majoraxis(ellipsoid(GRS80S))
    for (lat₁, lon₁, lat₂, lon₂) in [(-33.8688, 151.2093, 51.5074, -0.1278), (0, 0, 0.5, 179.7), (60, -45, -60, 135)]
      s₁, c₁ = sincosd(lat₁)
      s₂, c₂ = sincosd(lat₂)
      sΔ, cΔ = sincosd(lon₂ - lon₁)
      greatcircle = R * atan(hypot(c₂ * sΔ, c₁ * s₂ - s₁ * c₂ * cΔ), s₁ * s₂ + c₁ * c₂ * cΔ)
      geodesic = d(Point(LatLon{GRS80S}(lat₁, lon₁)), Point(LatLon{GRS80S}(lat₂, lon₂)))
      @test isapprox(geodesic, greatcircle, atol=τ)
    end

    # the ellipsoid comes from the datum, so a different datum gives a different answer
    q₁ = Point(LatLon{ITRF{2008}}(-33.8688, 151.2093))
    q₂ = Point(LatLon{ITRF{2008}}(51.5074, -0.1278))
    @test d(q₁, q₂) ≠ d(p₁, p₂)
    @test isapprox(d(q₁, q₂), d(p₁, p₂), atol=1e-3u"m")
  end
end

@testitem "Distances.jl" setup = [Setup] begin
  p = cart(0, 1)
  l = Line(cart(0, 0), cart(1, 0))
  @test evaluate(Euclidean(), p, l) == T(1) * u"m"
  @test evaluate(Euclidean(), l, p) == T(1) * u"m"
  p = cart(-3, 4)
  s = Segment(cart(0, 0), cart(1, 0))
  @test evaluate(Euclidean(), p, s) == T(5) * u"m"
  @test evaluate(Euclidean(), s, p) == T(5) * u"m"
  @test evaluate(Euclidean(), p, l) != T(5) * u"m"

  p = cart(68, 259)
  l = Line(cart(68, 260), cart(69, 261))
  @test evaluate(Euclidean(), p, l) ≤ T(0.8) * u"m"
  line1 = Line(cart(-1, 0, 0), cart(1, 0, 0))
  line2 = Line(cart(0, -1, 1), cart(0, 1, 1))  # line2 ⟂ line1, z++
  line3 = Line(cart(-1, 1, 0), cart(1, 1, 0))  # line3 ∥ line1
  line4 = Line(cart(-2, 0, 0), cart(2, 0, 0))  # line4 colinear with line1
  line5 = Line(cart(0, -1, 0), cart(0, 1, 0))  # line5 intersects line1
  line6 = Line(cart(0, -1, 0), cart(0, -2, 0))  # line6 intersects line1, if infinite
  @test evaluate(Euclidean(), line1, line2) ≈ T(1) * u"m"
  @test evaluate(Euclidean(), line1, line3) ≈ T(1) * u"m"
  @test evaluate(Euclidean(), line1, line4) ≈ T(0) * u"m"
  @test evaluate(Euclidean(), line1, line5) ≈ T(0) * u"m"
  @test evaluate(Euclidean(), line1, line6) ≈ T(0) * u"m"
  seg1 = Segment(cart(-1, 0, 0), cart(1, 0, 0))
  seg2 = Segment(cart(0, -1, 1), cart(0, 1, 1))  # seg2 ⟂ seg1, z++
  seg3 = Segment(cart(-1, 1, 0), cart(1, 1, 0))  # seg3 ∥ seg1
  seg4 = Segment(cart(-0, 0, 0), cart(2, 0, 0))  # seg4 colinear with seg1
  seg5 = Segment(cart(0, -1, 0), cart(0, 1, 0))  # seg5 intersects seg1
  seg6 = Segment(cart(0, -1, 0), cart(0, -2, 0))  # seg6 intersects seg1, if infinite
  seg7 = Segment(cart(2, 0, 0), cart(4, 0, 0))  # seg7 colinear with seg1 but shifted (gap of 1)
  seg8 = Segment(cart(3, -2, 0), cart(5, -2, 0))  # seg8 ∥ seg1 but offset (gap of 2 by 2=√8)
  @test evaluate(Euclidean(), seg1, seg2) ≈ T(1) * u"m"
  @test evaluate(Euclidean(), seg1, seg3) ≈ T(1) * u"m"
  @test evaluate(Euclidean(), seg1, seg4) ≈ T(0) * u"m"
  @test evaluate(Euclidean(), seg1, seg5) ≈ T(0) * u"m"
  @test evaluate(Euclidean(), seg1, seg6) ≈ T(1) * u"m"
  @test evaluate(Euclidean(), seg1, seg7) ≈ T(1) * u"m"
  @test evaluate(Euclidean(), seg1, seg8) ≈ T(√8) * u"m"

  p1, p2 = cart(1, 0), cart(0, 1)
  @test evaluate(Chebyshev(), p1, p2) == T(1) * u"m"
  @test evaluate(Euclidean(), p1, p2) == T(√2) * u"m"

  latlon1 = LatLon(T(0), T(0))
  latlon2 = LatLon(T(1), T(0))
  cart1 = convert(Cartesian, latlon1)
  cart2 = convert(Cartesian, latlon2)
  p1 = Point(latlon1)
  p2 = Point(latlon2)
  p3 = Point(cart1)
  p4 = Point(cart2)
  @test evaluate(Haversine(), p1, p2) ≈ T(111194.92664455874) * u"m"
  @test evaluate(Haversine(), p3, p4) ≈ T(111194.92664455874) * u"m"
  @test evaluate(Haversine(6371000u"m"), p1, p2) ≈ T(111194.92664455874) * u"m"
  @test evaluate(Haversine(6371000u"m"), p3, p4) ≈ T(111194.92664455874) * u"m"
  @test evaluate(Haversine(6371u"km"), p1, p2) ≈ T(111.19492664455874) * u"km"
  @test evaluate(Haversine(6371u"km"), p3, p4) ≈ T(111.19492664455874) * u"km"
  @test evaluate(SphericalAngle(), p1, p2) ≈ deg2rad(T(1) * u"°")
  @test evaluate(SphericalAngle(), p3, p4) ≈ deg2rad(T(1) * u"°")
end
