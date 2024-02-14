@testset "SubDomains" begin
  pset = PointSet(rand(P3, 100))
  inds = rand(1:100, 3)
  v = view(pset, inds)
  @test nelements(v) == 3
  for i in 1:3
    p = pset[inds[i]]
    @test v[i] == p
    @test centroid(v, i) == p
  end

  grid = CartesianGrid{T}(10, 10)
  inds = rand(1:100, 3)
  v = view(grid, inds)
  @test nelements(v) == 3
  for i in 1:3
    e = grid[inds[i]]
    @test v[i] == e
    @test centroid(v, i) == centroid(e)
  end

  points = P2[(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
  connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
  mesh = SimpleMesh(points, connec)
  inds = rand(1:4, 3)
  v = view(mesh, inds)
  @test nelements(v) == 3
  for i in 1:3
    e = mesh[inds[i]]
    @test v[i] == e
    @test centroid(v, i) == centroid(e)
  end

  # view of view stores the correct domain
  g = CartesianGrid{T}(10, 10)
  v = view(view(g, 11:20), 1:3)
  @test v isa Meshes.SubDomain{2,T,<:CartesianGrid}
  @test v[1] == g[11]
  @test v[2] == g[12]
  @test v[3] == g[13]

  # centroid of view of PointSet
  points = P2[(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
  pview = view(PointSet(points), 1:4)
  @test centroid(pview) == P2(0.5, 0.5)

  # measure of view
  g = CartesianGrid{T}(10, 10)
  v = view(g, 1:3)
  @test measure(v) ≈ T(3)

  # concatenation with same parent
  g = CartesianGrid{T}(10, 10)
  vg = vcat(view(g, 50:70), view(g, 10:30))
  @test vg isa Meshes.SubDomain
  @test vg == view(g, [50:70; 10:30])
  # concatenation with different parents
  g1 = CartesianGrid{T}(10, 10)
  g2 = CartesianGrid{T}(20, 20)
  vg = vcat(view(g1, 50:70), view(g2, 10:30))
  @test vg isa GeometrySet
  @test vg == GeometrySet([g1[50:70]; g2[10:30]])

  # eltype
  d1 = CartesianGrid{T}(10, 10)
  d2 = GeometrySet([P2(0, 0), Box(P2(0, 0), P2(1, 1)), P2(2, 2)])
  v1 = view(d1, 1:10)
  v2 = view(d2, [1, 3])
  @test eltype(v1) === Quadrangle{2,T}
  @test eltype(v2) === Primitive{2,T}

  # show
  pset = PointSet(P2.(1:100, 1:100))
  v1 = view(pset, 1:10)
  v2 = view(pset, [4, 8, 10, 7, 9, 1, 2, 3, 6, 5])
  @test sprint(show, v1) == "10 view(::PointSet{2,$T}, 1:10)"
  @test sprint(show, v2) == "10 view(::PointSet{2,$T}, [4, 8, 10, 7, ..., 2, 3, 6, 5])"
  if T === Float32
    @test sprint(show, MIME"text/plain"(), v1) == """
    10 view(::PointSet{2,Float32}, 1:10)
    ├─ Point(1.0f0, 1.0f0)
    ├─ Point(2.0f0, 2.0f0)
    ├─ Point(3.0f0, 3.0f0)
    ├─ Point(4.0f0, 4.0f0)
    ├─ Point(5.0f0, 5.0f0)
    ├─ Point(6.0f0, 6.0f0)
    ├─ Point(7.0f0, 7.0f0)
    ├─ Point(8.0f0, 8.0f0)
    ├─ Point(9.0f0, 9.0f0)
    └─ Point(10.0f0, 10.0f0)"""
    @test sprint(show, MIME"text/plain"(), v2) == """
    10 view(::PointSet{2,Float32}, [4, 8, 10, 7, ..., 2, 3, 6, 5])
    ├─ Point(4.0f0, 4.0f0)
    ├─ Point(8.0f0, 8.0f0)
    ├─ Point(10.0f0, 10.0f0)
    ├─ Point(7.0f0, 7.0f0)
    ├─ Point(9.0f0, 9.0f0)
    ├─ Point(1.0f0, 1.0f0)
    ├─ Point(2.0f0, 2.0f0)
    ├─ Point(3.0f0, 3.0f0)
    ├─ Point(6.0f0, 6.0f0)
    └─ Point(5.0f0, 5.0f0)"""
  else
    @test sprint(show, MIME"text/plain"(), v1) == """
    10 view(::PointSet{2,Float64}, 1:10)
    ├─ Point(1.0, 1.0)
    ├─ Point(2.0, 2.0)
    ├─ Point(3.0, 3.0)
    ├─ Point(4.0, 4.0)
    ├─ Point(5.0, 5.0)
    ├─ Point(6.0, 6.0)
    ├─ Point(7.0, 7.0)
    ├─ Point(8.0, 8.0)
    ├─ Point(9.0, 9.0)
    └─ Point(10.0, 10.0)"""
    @test sprint(show, MIME"text/plain"(), v2) == """
    10 view(::PointSet{2,Float64}, [4, 8, 10, 7, ..., 2, 3, 6, 5])
    ├─ Point(4.0, 4.0)
    ├─ Point(8.0, 8.0)
    ├─ Point(10.0, 10.0)
    ├─ Point(7.0, 7.0)
    ├─ Point(9.0, 9.0)
    ├─ Point(1.0, 1.0)
    ├─ Point(2.0, 2.0)
    ├─ Point(3.0, 3.0)
    ├─ Point(6.0, 6.0)
    └─ Point(5.0, 5.0)"""
  end
end
