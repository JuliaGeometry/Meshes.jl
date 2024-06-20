@testset "Tesselation" begin
  @testset "Delaunay" begin
    pts = randpoint2(10)
    pset = PointSet(pts)
    mesh1 = tesselate(pts, DelaunayTesselation(StableRNG(2024)))
    mesh2 = tesselate(pset, DelaunayTesselation(StableRNG(2024)))
    @test mesh1 == mesh2

    # datum and unit propagation
    tuples = [(rand(T) * u"km", rand(T) * u"km") for _ in 1:10]
    pset = PointSet(Point.(Cartesian{WGS84Latest}.(tuples)))
    mesh = tesselate(pset, DelaunayTesselation(StableRNG(2024)))
    @test datum(Meshes.crs(mesh)) === WGS84Latest
    @test unit(Meshes.lentype(mesh)) == u"km"
  end

  @testset "Voronoi" begin
    pts = randpoint2(10)
    pset = PointSet(pts)
    mesh1 = tesselate(pts, VoronoiTesselation(StableRNG(2024)))
    mesh2 = tesselate(pset, VoronoiTesselation(StableRNG(2024)))
    @test mesh1 == mesh2

    # datum and unit propagation
    tuples = [(rand(T) * u"km", rand(T) * u"km") for _ in 1:10]
    pset = PointSet(Point.(Cartesian{WGS84Latest}.(tuples)))
    mesh = tesselate(pset, VoronoiTesselation(StableRNG(2024)))
    @test datum(Meshes.crs(mesh)) === WGS84Latest
    @test unit(Meshes.lentype(mesh)) == u"km"
  end
end
