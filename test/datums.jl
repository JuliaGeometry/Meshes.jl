@testset "Datums" begin
  @testset "NoDatum" begin
    @test isnothing(ellipsoid(NoDatum))
    @test isnothing(latitudeₒ(NoDatum))
    @test isnothing(longitudeₒ(NoDatum))
    @test isnothing(altitudeₒ(NoDatum))
  end

  @testset "WGS84" begin
    🌎 = ellipsoid(WGS84)
    @test majoraxis(🌎) == 6378137.0u"m"
    @test minoraxis(🌎) == 6356752.314245179u"m"
    @test eccentricity(🌎) == 0.08181919084262149
    @test eccentricity²(🌎) == 0.0066943799901413165
    @test flattening(🌎) == 0.0033528106647474805
    @test flattening⁻¹(🌎) == 298.257223563

    @test latitudeₒ(WGS84) == 0.0u"°"
    @test longitudeₒ(WGS84) == 0.0u"°"
    @test altitudeₒ(WGS84) == 0.0u"m"
  end
end
