@testset "Views" begin
  pset = PointSet(rand(P3, 100))
  inds = rand(1:100, 3)
  v = view(pset, inds)
  @test nelements(v) == 3
  for i in 1:3
    p = pset[inds[i]]
    @test v[i] == p
    @test coordinates(v, i) == coordinates(p)
  end

  grid = CartesianGrid{T}(10, 10)
  inds = rand(1:100, 3)
  v = view(grid, inds)
  @test nelements(v) == 3
  for i in 1:3
    e = grid[inds[i]]
    @test v[i] == e
    @test coordinates(v, i) == coordinates(centroid(e))
  end

  points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
  connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
  mesh = UnstructuredMesh(points, connec)
  inds = rand(1:4, 3)
  v = view(mesh, inds)
  @test nelements(v) == 3
  for i in 1:3
    e = mesh[inds[i]]
    @test v[i] == e
    @test coordinates(v, i) == coordinates(centroid(e))
  end

  if visualtests
    d = CartesianGrid{T}(10, 10)
    v = view(d, 1:50)
    @test_ref_plot "data/domain-view-$T.png" plot(v)
  end
end
