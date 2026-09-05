@testitem "Geodesics" setup = [Setup] begin
  # the geodesics are accurate to nanometres, which is far coarser than the
  # picometre default tolerance that Meshes uses to compare points
  τ = T === Float64 ? 1e-7u"m" : 10u"m"

  # the azimuth is only defined for points on the ellipsoid
  @test_throws MethodError geodesicfwd(cart(0, 0), 90, 1000)
  @test_throws MethodError geodesicbwd(cart(0, 0), cart(1, 1))

  # walking nowhere leaves the point where it is
  p = latlon(45, 10)
  @test geodesicfwd(p, 30, 0) ≈ p

  # a quarter of the equator to the east
  a = CoordRefSystems.majoraxis(ellipsoid(WGS84Latest))
  @test isapprox(geodesicfwd(latlon(0, 0), 90u"°", π * a / 2), latlon(0, 90), atol=τ)

  # the azimuth of a meridian is due north or due south
  @test geodesicbwd(latlon(10, 20), latlon(30, 20)) ≈ T(0) * u"°" atol = 1e-6u"°"
  @test geodesicbwd(latlon(30, 20), latlon(10, 20)) ≈ T(180) * u"°" atol = 1e-6u"°"

  # the azimuth along the equator is due east or due west
  @test geodesicbwd(latlon(0, 0), latlon(0, 10)) ≈ T(90) * u"°" atol = 1e-6u"°"
  @test geodesicbwd(latlon(0, 10), latlon(0, 0)) ≈ T(-90) * u"°" atol = 1e-6u"°"

  # units are optional, and the result carries them
  @test geodesicfwd(p, 30u"°", 100u"km") ≈ geodesicfwd(p, 30, 100000)
  ϕ = geodesicbwd(p, latlon(46, 11))
  @test unit(ϕ) == u"°"
  @test Unitful.numtype(ϕ) === T

  # walking along the azimuth that connects two points lands on the second one
  p₁ = latlon(-33.8688, 151.2093)
  p₂ = latlon(51.5074, -0.1278)
  @test isapprox(geodesicfwd(p₁, geodesicbwd(p₁, p₂), GeodesicDistance()(p₁, p₂)), p₂, atol=τ)

  # the ellipsoid comes from the datum
  q₁ = Point(LatLon{ITRF{2008}}(-33.8688, 151.2093))
  q₂ = Point(LatLon{ITRF{2008}}(51.5074, -0.1278))
  @test geodesicbwd(q₁, q₂) ≠ geodesicbwd(p₁, p₂)

  if T === Float64
    # reference values from the test set of Karney (2013)
    # https://geographiclib.sourceforge.io/C++/doc/geodesic.html#testgeod
    τϕ = 1e-9u"°"
    cases = [
      (
        2.881248229541,
        0.0,
        27.763592972746,
        53.997072295385487,
        44.520619105667620,
        52.159486739947740,
        6958264.1576889
      ),
      (
        65.656162297631,
        0.0,
        176.971135321064,
        -6.529066987956306,
        2.895923948124536,
        178.740350145953805,
        8009999.3798375
      ),
      (
        75.511482283510,
        0.0,
        83.078727908415,
        55.600487151982554,
        75.128743229495482,
        153.896688535571762,
        3723062.6140266
      )
    ]
    for (lat₁, lon₁, azi₁, lat₂, lon₂, azi₂, s₁₂) in cases
      p₁ = latlon(lat₁, lon₁)
      p₂ = latlon(lat₂, lon₂)
      # inverse problem: the azimuth at each end
      @test isapprox(geodesicbwd(p₁, p₂), azi₁ * u"°", atol=τϕ)
      @test isapprox(geodesicbwd(p₂, p₁), (azi₂ - 180) * u"°", atol=τϕ)
      # direct problem: the point reached and the distance to it
      @test isapprox(geodesicfwd(p₁, azi₁, s₁₂), p₂, atol=τ)
      @test isapprox(GeodesicDistance()(p₁, geodesicfwd(p₁, azi₁, s₁₂)), s₁₂ * u"m", atol=1e-7u"m")
    end

    # on a datum with a spherical ellipsoid the azimuth is the great circle one
    r₁ = Point(LatLon{GRS80S}(10, 20))
    r₂ = Point(LatLon{GRS80S}(30, 50))
    Δ = deg2rad(50 - 20)
    greatcircle = atand(
      cos(deg2rad(30)) * sin(Δ),
      cos(deg2rad(10)) * sin(deg2rad(30)) - sin(deg2rad(10)) * cos(deg2rad(30)) * cos(Δ)
    )
    @test isapprox(geodesicbwd(r₁, r₂), greatcircle * u"°", atol=1e-9u"°")
  end
end
