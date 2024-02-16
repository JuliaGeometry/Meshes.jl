@testset "EPSG/ESRI Codes" begin
  @test Meshes.crs(EPSG{3395}) === Mercator{WGS84}
  @test Meshes.crs(EPSG{3857}) === WebMercator{WGS84}
  @test Meshes.crs(EPSG{4326}) === LatLon{WGS84}
  @test Meshes.crs(EPSG{32662}) === PlateCarree{WGS84}
  @test Meshes.crs(ESRI{54017}) === Behrmann{WGS84}
  @test Meshes.crs(ESRI{54030}) === Robinson{WGS84}
  @test Meshes.crs(ESRI{54034}) === Lambert{WGS84}
  @test Meshes.crs(ESRI{54042}) === WinkelTripel{WGS84}
  @test Meshes.crs(ESRI{102035}) === Meshes.Orthographic{90.0u"°",0.0u"°",true,WGS84}
  @test Meshes.crs(ESRI{102037}) === Meshes.Orthographic{-90.0u"°",0.0u"°",true,WGS84}
end
