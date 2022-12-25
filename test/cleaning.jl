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

  @testset "Merging of duplicated vertices" begin
    # make a tetrahedron with duplicated vertices
    p1 = P3(0, 1, 1)
    p2 = P3(-1, 2, 3)
    p3 = P3(0, 3, 2)
    p4 = P3(2, 2, 2)
    points = [p1, p2, p3, p3, p2, p4, p4, p2, p1, p1, p3, p4]
    triangles = [(1, 2, 3), (4, 5, 6), (7, 8, 9), (10, 11, 12)]
    mesh = SimpleMesh(points, connect.(triangles, Triangle))

    # merge the duplicated vertices
    cmesh = gather(mesh)

    # the new mesh has four vertices
    @test nvertices(cmesh) == 4

    # the new mesh has four faces
    @test nelements(cmesh) == 4
  end
end