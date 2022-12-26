@testset "Connected components" begin
  # make a mesh of two disconnected tetrahedra
  connec = connect.([(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)], Triangle)
  vs1 = [P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0), P3(0, 0, 1)]
  mesh1 = SimpleMesh(vs1, connec)
  vs2 = [P3(0, 0, 5), P3(1, 0, 5), P3(0, 1, 5), P3(0, 0, 1)]
  mesh2 = SimpleMesh(vs2, connec)
  mesh = merge(mesh1, mesh2)
  # compute the connected components
  meshes = connected_components(mesh)
  # check there are two connected components
  @test length(meshes) == 2
  # check each component is a tetrahedron
  @test nvertices(meshes[1]) == 4
  @test nelements(meshes[1]) == 4
  @test nvertices(meshes[2]) == 4
  @test nelements(meshes[2]) == 4
end
