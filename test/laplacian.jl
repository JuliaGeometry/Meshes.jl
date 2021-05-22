@testset "Laplacian" begin
  points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
  connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
  mesh = SimpleMesh(points, connec)
  L = laplacematrix(mesh, weights=:uniform, normalize=true)
  @test L == [
     -1 1/3 1/3   0 1/3
    1/3  -1   0 1/3 1/3
    1/3   0  -1 1/3 1/3
      0 1/3 1/3  -1 1/3
    1/4 1/4 1/4 1/4  -1
  ]
end
