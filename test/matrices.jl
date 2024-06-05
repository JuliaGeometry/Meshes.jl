@testset "Matrices" begin
  # uniform weights for simple mesh
  points = point.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
  connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
  mesh = SimpleMesh(points, connec)
  L = laplacematrix(mesh, kind=:uniform)
  @test L == [
    -1 1/3 1/3 0 1/3
    1/3 -1 0 1/3 1/3
    1/3 0 -1 1/3 1/3
    0 1/3 1/3 -1 1/3
    1/4 1/4 1/4 1/4 -1
  ]

  # cotangent weights only defined for triangle meshes
  points = point.([(0, 0), (1, 0), (1, 1), (0, 1)])
  connec = connect.([(1, 2, 3, 4)], Quadrangle)
  mesh = SimpleMesh(points, connec)
  @test_throws AssertionError laplacematrix(mesh, kind=:cotangent)

  # full Laplace-Beltrami operator
  sphere = Sphere(point(0, 0, 0), T(1))
  mesh = simplexify(sphere)
  L = laplacematrix(mesh)
  M = measurematrix(mesh)
  @test issymmetric(L)
  @test issparse(L)
  @test isdiag(M)

  # adjacency of CartesianGrid
  grid = cartgrid(100, 100)
  A = adjacencymatrix(grid)
  d = sum(A, dims=2)
  @test size(A) == (10000, 10000)
  @test issymmetric(A)
  @test issparse(A)
  @test minimum(d) == 2
  @test maximum(d) == 4
  @test length(findall(==(2), d)) == 4
  A = adjacencymatrix(grid, rank=0)
  @test size(A) == (101*101, 101*101)

  # adjacency of SimpleMesh
  points = point.([(0, 0), (1, -1), (1, 1), (2, -1), (2, 1)])
  connec = connect.([(1, 2, 3), (3, 2, 4, 5)])
  mesh = SimpleMesh(points, connec, relations=true)
  A = adjacencymatrix(mesh)
  @test A == [0 1; 1 0]
  A = adjacencymatrix(mesh, rank=0)
  @test A == [
    0 1 1 0 0
    1 0 1 1 0
    1 1 0 0 1
    0 1 0 0 1
    0 0 1 1 0
  ]
end
