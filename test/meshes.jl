@testset "Meshes" begin
  @testset "Unstructured" begin
    points = Point2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = UnstructuredMesh(points, connec)
    triangles = Triangle.([
      [(0.0,0.0), (1.0,0.0), (0.5,0.5)],
      [(1.0,0.0), (1.0,1.0), (0.5,0.5)],
      [(1.0,1.0), (0.0,1.0), (0.5,0.5)],
      [(0.0,1.0), (0.0,0.0), (0.5,0.5)]
    ])
    bytes = @allocated faces(mesh, 2)
    @test bytes < 100
    cells = faces(mesh, 2)
    bytes = @allocated collect(cells)
    @test bytes < 800
    @test collect(cells) == triangles
  end
end
