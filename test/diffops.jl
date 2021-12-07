@testset "DiffOps" begin
  # uniform weights for simple mesh
  points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
  connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
  mesh = SimpleMesh(points, connec)
  L = laplacematrix(mesh, weights=:uniform)
  @test L == [
     -1 1/3 1/3   0 1/3
    1/3  -1   0 1/3 1/3
    1/3   0  -1 1/3 1/3
      0 1/3 1/3  -1 1/3
    1/4 1/4 1/4 1/4  -1
  ]

  # cotangent weights only defined for triangle meshes
  points = P2[(0,0), (1,0), (1,1), (0,1)]
  connec = connect.([(1,2,3,4)], Quadrangle)
  mesh = SimpleMesh(points, connec)
  @test_throws AssertionError laplacematrix(mesh, weights=:cotangent)

  # full Laplace-Beltrami operator
  sphere = Sphere(P3(0,0,0), T(1))
  mesh = triangulate(sphere)
  L = laplacematrix(mesh)
  M = measurematrix(mesh)
  @test issymmetric(L)
  @test issparse(L)
  @test isdiag(M)
end
