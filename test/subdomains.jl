@testset "SubDomains" begin
  pset = PointSet(randpoint3(100))
  inds = rand(1:100, 3)
  v = view(pset, inds)
  @test nelements(v) == 3
  @test Meshes.crs(v) <: Cartesian{NoDatum}
  @test Meshes.lentype(v) == ℳ
  for i in 1:3
    p = pset[inds[i]]
    @test v[i] == p
    @test centroid(v, i) == p
  end

  grid = cartgrid(10, 10)
  inds = rand(1:100, 3)
  v = view(grid, inds)
  @test nelements(v) == 3
  @test Meshes.crs(v) <: Cartesian{NoDatum}
  @test Meshes.lentype(v) == ℳ
  for i in 1:3
    e = grid[inds[i]]
    @test v[i] == e
    @test centroid(v, i) == centroid(e)
  end

  points = point.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
  connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
  mesh = SimpleMesh(points, connec)
  inds = rand(1:4, 3)
  v = view(mesh, inds)
  @test nelements(v) == 3
  @test Meshes.crs(v) <: Cartesian{NoDatum}
  @test Meshes.lentype(v) == ℳ
  for i in 1:3
    e = mesh[inds[i]]
    @test v[i] == e
    @test centroid(v, i) == centroid(e)
  end

  # view of view stores the correct domain
  g = cartgrid(10, 10)
  v = view(view(g, 11:20), 1:3)
  @test v isa SubGrid{2}
  @test v[1] == g[11]
  @test v[2] == g[12]
  @test v[3] == g[13]

  # centroid of view of PointSet
  points = point.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
  pview = view(PointSet(points), 1:4)
  @test centroid(pview) == point(0.5, 0.5)

  # measure of view
  g = cartgrid(10, 10)
  v = view(g, 1:3)
  @test measure(v) ≈ T(3) * u"m^2"

  # concatenation with same parent
  g = cartgrid(10, 10)
  vg = vcat(view(g, 50:70), view(g, 10:30))
  @test vg isa SubGrid{2}
  @test vg == view(g, [50:70; 10:30])
  # concatenation with different parents
  g1 = cartgrid(10, 10)
  g2 = cartgrid(20, 20)
  vg = vcat(view(g1, 50:70), view(g2, 10:30))
  @test vg isa GeometrySet
  @test vg == GeometrySet([g1[50:70]; g2[10:30]])

  # eltype
  d1 = cartgrid(1000, 1000)
  d2 = cartgrid(1000, 1000, 1000)
  d3 = GeometrySet([point(0, 0), Box(point(0, 0), point(1, 1)), point(2, 2)])
  v1 = view(d1, 1:500000)
  v2 = view(d2, 1:500000000)
  v3 = view(d3, [1, 3])
  @test eltype(v1) <: Quadrangle{2}
  @test eltype(v2) <: Hexahedron{3}
  @test eltype(v3) <: Primitive{2}

  # show
  pset = PointSet(point.(1:100, 1:100))
  v1 = view(pset, 1:10)
  v2 = view(pset, [4, 8, 10, 7, 9, 1, 2, 3, 6, 5])
  @test sprint(show, v1) == "10 view(::PointSet, 1:10)"
  @test sprint(show, v2) == "10 view(::PointSet, [4, 8, 10, 7, ..., 2, 3, 6, 5])"
  if T === Float32
    @test sprint(show, MIME"text/plain"(), v1) == """
    10 view(::PointSet, 1:10)
    ├─ Point(x: 1.0f0 m, y: 1.0f0 m)
    ├─ Point(x: 2.0f0 m, y: 2.0f0 m)
    ├─ Point(x: 3.0f0 m, y: 3.0f0 m)
    ├─ Point(x: 4.0f0 m, y: 4.0f0 m)
    ├─ Point(x: 5.0f0 m, y: 5.0f0 m)
    ├─ Point(x: 6.0f0 m, y: 6.0f0 m)
    ├─ Point(x: 7.0f0 m, y: 7.0f0 m)
    ├─ Point(x: 8.0f0 m, y: 8.0f0 m)
    ├─ Point(x: 9.0f0 m, y: 9.0f0 m)
    └─ Point(x: 10.0f0 m, y: 10.0f0 m)"""
    @test sprint(show, MIME"text/plain"(), v2) == """
    10 view(::PointSet, [4, 8, 10, 7, ..., 2, 3, 6, 5])
    ├─ Point(x: 4.0f0 m, y: 4.0f0 m)
    ├─ Point(x: 8.0f0 m, y: 8.0f0 m)
    ├─ Point(x: 10.0f0 m, y: 10.0f0 m)
    ├─ Point(x: 7.0f0 m, y: 7.0f0 m)
    ├─ Point(x: 9.0f0 m, y: 9.0f0 m)
    ├─ Point(x: 1.0f0 m, y: 1.0f0 m)
    ├─ Point(x: 2.0f0 m, y: 2.0f0 m)
    ├─ Point(x: 3.0f0 m, y: 3.0f0 m)
    ├─ Point(x: 6.0f0 m, y: 6.0f0 m)
    └─ Point(x: 5.0f0 m, y: 5.0f0 m)"""
  else
    @test sprint(show, MIME"text/plain"(), v1) == """
    10 view(::PointSet, 1:10)
    ├─ Point(x: 1.0 m, y: 1.0 m)
    ├─ Point(x: 2.0 m, y: 2.0 m)
    ├─ Point(x: 3.0 m, y: 3.0 m)
    ├─ Point(x: 4.0 m, y: 4.0 m)
    ├─ Point(x: 5.0 m, y: 5.0 m)
    ├─ Point(x: 6.0 m, y: 6.0 m)
    ├─ Point(x: 7.0 m, y: 7.0 m)
    ├─ Point(x: 8.0 m, y: 8.0 m)
    ├─ Point(x: 9.0 m, y: 9.0 m)
    └─ Point(x: 10.0 m, y: 10.0 m)"""
    @test sprint(show, MIME"text/plain"(), v2) == """
    10 view(::PointSet, [4, 8, 10, 7, ..., 2, 3, 6, 5])
    ├─ Point(x: 4.0 m, y: 4.0 m)
    ├─ Point(x: 8.0 m, y: 8.0 m)
    ├─ Point(x: 10.0 m, y: 10.0 m)
    ├─ Point(x: 7.0 m, y: 7.0 m)
    ├─ Point(x: 9.0 m, y: 9.0 m)
    ├─ Point(x: 1.0 m, y: 1.0 m)
    ├─ Point(x: 2.0 m, y: 2.0 m)
    ├─ Point(x: 3.0 m, y: 3.0 m)
    ├─ Point(x: 6.0 m, y: 6.0 m)
    └─ Point(x: 5.0 m, y: 5.0 m)"""
  end
end
