@testset "Cleaning" begin
  @testset "Removal of unused vertices" begin
    points = P3[(0, 0, 0), (0, 0, 1), (5, 5, 5), (0, 1, 0), (1, 0, 0)]
    triangles = [(1, 2, 4), (1, 2, 5), (1, 4, 5), (2, 4, 5)] # point 3 is not used
    mesh = SimpleMesh(points, connect.(triangles, Triangle))
    cmesh = clean(mesh)
    newpoints = vertices(cmesh)

    # only point 3 has been removed
    @test nvertices(cmesh) == 4
    @test points[3] ∉ newpoints
    @test all(points[[1,2,4,5]] .∈ (newpoints,))

    # number of faces has not changed
    @test nelements(cmesh) == 4
  end
end