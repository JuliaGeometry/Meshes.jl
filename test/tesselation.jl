@testset "Tesselation" begin
  @testset "Delaunay" begin
    pts = randpoint2(10)
    pset = PointSet(pts)
    mesh1 = tesselate(pts, DelaunayTesselation(StableRNG(2024)))
    mesh2 = tesselate(pset, DelaunayTesselation(StableRNG(2024)))
    @test mesh1 == mesh2

    # CRS propagation
    tuples = [(rand(T) * u"km", rand(T) * u"km") for _ in 1:10]
    pset = PointSet(Point.(Cartesian{WGS84Latest}.(tuples)))
    mesh = tesselate(pset, DelaunayTesselation(StableRNG(2024)))
    @test crs(mesh) === crs(pset)

    coords = [LatLon(rand(-90:T(0.1):90), rand(-180:T(0.1):180)) for _ in 1:10]
    pset = PointSet(Point.(coords))
    mesh = tesselate(pset, DelaunayTesselation(StableRNG(2024)))
    @test crs(mesh) === crs(pset)

    # error: the number of coordinates of the points must be 2
    pts = randpoint3(10)
    pset = PointSet(pts)
    @test_throws AssertionError tesselate(pset, DelaunayTesselation(StableRNG(2024)))
  end

  @testset "Voronoi" begin
    pts = randpoint2(10)
    pset = PointSet(pts)
    mesh1 = tesselate(pts, VoronoiTesselation(StableRNG(2024)))
    mesh2 = tesselate(pset, VoronoiTesselation(StableRNG(2024)))
    @test mesh1 == mesh2

    # CRS propagation
    tuples = [(rand(T) * u"km", rand(T) * u"km") for _ in 1:10]
    pset = PointSet(Point.(Cartesian{WGS84Latest}.(tuples)))
    mesh = tesselate(pset, VoronoiTesselation(StableRNG(2024)))
    @test crs(mesh) === crs(pset)

    coords = [LatLon(rand(-90:T(0.1):90), rand(-180:T(0.1):180)) for _ in 1:10]
    pset = PointSet(Point.(coords))
    mesh = tesselate(pset, VoronoiTesselation(StableRNG(2024)))
    @test crs(mesh) === crs(pset)

    # error: the number of coordinates of the points must be 2
    pts = randpoint3(10)
    pset = PointSet(pts)
    @test_throws AssertionError tesselate(pset, VoronoiTesselation(StableRNG(2024)))

    # Test polygon order is the same as input points order
    pts = randpoint2(10)
    mesh = tesselate(pts, VoronoiTesselation(StableRNG(2024)))
    @test all(zip(pts, mesh)) do (pt, poly)
      pt in poly && return true
      # Point is not in poly, might be to a rounding error. We check if the target polygon's centroid is the closest of all
      centroid_dists = map(mesh) do element
        norm(centroid(element) - pt)
      end
      this_dist = norm(centroid(poly) - pt)
      all(<=(this_dist), centroid_dists)
    end
  end
end
