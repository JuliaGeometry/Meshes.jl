@testset "Tesselation" begin
  @testset "Delaunay" begin
    pts = randpoint2(100)
    pset = PointSet(pts)
    mesh1 = tesselate(pset, DelaunayTesselation(StableRNG(2024)))
    mesh2 = tesselate(pts, DelaunayTesselation(StableRNG(2024)))
    @test mesh1 == mesh2
  end
end
