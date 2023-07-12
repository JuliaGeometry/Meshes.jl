@testset "Meshes" begin
  @testset "CartesianGrid" begin
    grid = CartesianGrid{T}(100)
    @test embeddim(grid) == 1
    @test coordtype(grid) == T
    @test size(grid) == (100,)
    @test minimum(grid) == P1(0)
    @test maximum(grid) == P1(100)
    @test extrema(grid) == (P1(0), P1(100))
    @test spacing(grid) == T.((1,))
    @test nelements(grid) == 100
    @test eltype(grid) <: Segment{1,T}
    @test measure(grid) ≈ T(100)

    grid = CartesianGrid{T}(200, 100)
    @test embeddim(grid) == 2
    @test coordtype(grid) == T
    @test size(grid) == (200, 100)
    @test minimum(grid) == P2(0, 0)
    @test maximum(grid) == P2(200, 100)
    @test extrema(grid) == (P2(0, 0), P2(200, 100))
    @test spacing(grid) == T.((1, 1))
    @test nelements(grid) == 200 * 100
    @test eltype(grid) <: Quadrangle{2,T}
    @test measure(grid) ≈ T(200 * 100)

    grid = CartesianGrid((200, 100, 50), T.((0, 0, 0)), T.((1, 1, 1)))
    @test embeddim(grid) == 3
    @test coordtype(grid) == T
    @test size(grid) == (200, 100, 50)
    @test minimum(grid) == P3(0, 0, 0)
    @test maximum(grid) == P3(200, 100, 50)
    @test extrema(grid) == (P3(0, 0, 0), P3(200, 100, 50))
    @test spacing(grid) == T.((1, 1, 1))
    @test nelements(grid) == 200 * 100 * 50
    @test eltype(grid) <: Hexahedron{3,T}
    @test measure(grid) ≈ T(200 * 100 * 50)

    grid = CartesianGrid(T.((0, 0, 0)), T.((1, 1, 1)), T.((0.1, 0.1, 0.1)))
    @test embeddim(grid) == 3
    @test coordtype(grid) == T
    @test size(grid) == (10, 10, 10)
    @test minimum(grid) == P3(0, 0, 0)
    @test maximum(grid) == P3(1, 1, 1)
    @test spacing(grid) == T.((0.1, 0.1, 0.1))

    grid = CartesianGrid(T.((-1.0, -1.0)), T.((1.0, 1.0)), dims=(200, 100))
    @test embeddim(grid) == 2
    @test coordtype(grid) == T
    @test size(grid) == (200, 100)
    @test minimum(grid) == P2(-1.0, -1.0)
    @test maximum(grid) == P2(1.0, 1.0)
    @test spacing(grid) == T.((2 / 200, 2 / 100))
    @test nelements(grid) == 200 * 100
    @test eltype(grid) <: Quadrangle{2,T}

    grid = CartesianGrid((20, 10, 5), T.((0, 0, 0)), T.((5, 5, 5)))
    @test embeddim(grid) == 3
    @test coordtype(grid) == T
    @test size(grid) == (20, 10, 5)
    @test minimum(grid) == P3(0, 0, 0)
    @test maximum(grid) == P3(100, 50, 25)
    @test extrema(grid) == (P3(0, 0, 0), P3(100, 50, 25))
    @test spacing(grid) == T.((5, 5, 5))
    @test nelements(grid) == 20 * 10 * 5
    @test eltype(grid) <: Hexahedron{3,T}
    @test vertices(grid[1]) == 
          (P3(0, 0, 0), P3(5, 0, 0), P3(5, 5, 0), P3(0, 5, 0), P3(0, 0, 5), P3(5, 0, 5), P3(5, 5, 5), P3(0, 5, 5))
    @test all(centroid(grid, i) == centroid(grid[i]) for i in 1:nelements(grid))

    # constructor with offset
    grid = CartesianGrid((10, 10), T.((1.0, 1.0)), T.((1.0, 1.0)), (2, 2))
    @test embeddim(grid) == 2
    @test coordtype(grid) == T
    @test size(grid) == (10, 10)
    @test minimum(grid) == P2(0.0, 0.0)
    @test maximum(grid) == P2(10.0, 10.0)
    @test spacing(grid) == T.((1, 1))
    @test nelements(grid) == 10 * 10
    @test eltype(grid) <: Quadrangle{2,T}

    # indexing into a subgrid
    grid = CartesianGrid{T}(10, 10)
    sub = grid[1:2, 1:2]
    @test size(sub) == (2, 2)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == minimum(grid)
    @test maximum(sub) == P2(2, 2)
    sub = grid[1:1, 2:3]
    @test size(sub) == (1, 2)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == P2(0, 1)
    @test maximum(sub) == P2(1, 3)
    sub = grid[2:4, 3:7]
    @test size(sub) == (3, 5)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == P2(1, 2)
    @test maximum(sub) == P2(4, 7)
    grid = CartesianGrid(P2(1, 1), P2(11, 11), dims=(10, 10))
    sub = grid[2:4, 3:7]
    @test size(sub) == (3, 5)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == P2(2, 3)
    @test maximum(sub) == P2(5, 8)

    # subgrid with comparable vertices of grid
    grid = CartesianGrid((10, 10), P2(0.0, 0.0), T.((1.2, 1.2)))
    sub = grid[2:4, 5:7]
    @test sub == CartesianGrid((3, 3), P2(0.0, 0.0), T.((1.2, 1.2)), (0, -3))
    ind = reshape(reshape(1:121, 11, 11)[2:5, 5:8], :)
    @test vertices(grid)[ind] == vertices(sub)

    # subgrid from Cartesian ranges
    grid = CartesianGrid{T}(10, 10)
    sub1 = grid[1:2, 4:6]
    sub2 = grid[CartesianIndex(1, 4):CartesianIndex(2, 6)]
    @test sub1 == sub2

    grid = CartesianGrid{T}(200, 100)
    @test centroid(grid, 1) == P2(0.5, 0.5)
    @test centroid(grid, 2) == P2(1.5, 0.5)
    @test centroid(grid, 200 * 100) == P2(199.5, 99.5)
    @test nelements(grid) == 200 * 100
    @test eltype(grid) <: Quadrangle{2,T}
    @test grid[1] == Quadrangle(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    @test grid[2] == Quadrangle(P2[(1, 0), (2, 0), (2, 1), (1, 1)])

    # expand CartesianGrid with comparable vertices
    grid = CartesianGrid((10, 10), P2(0.0, 0.0), T.((1.0, 1.0)))
    left, right = (1, 1), (1, 1)
    newdim = size(grid) .+ left .+ right
    newoffset = offset(grid) .+ left
    grid2 = CartesianGrid(newdim, minimum(grid), spacing(grid), newoffset)
    @test issubset(vertices(grid), vertices(grid2))

    # GridTopology from CartesianGrid
    grid = CartesianGrid{T}(5, 5)
    topo = topology(grid)
    vs = vertices(grid)
    for i in 1:nelements(grid)
      inds = indices(element(topo, i))
      v = vertices(element(grid, i))
      @test vs[[inds...]] == collect(v)
    end

    # convert topology
    grid = CartesianGrid{T}(10, 10)
    mesh = topoconvert(HalfEdgeTopology, grid)
    @test mesh isa SimpleMesh
    @test nvertices(mesh) == 121
    @test nelements(mesh) == 100
    @test eltype(mesh) <: Quadrangle

    # single vertex access
    grid = CartesianGrid{T}(10, 10)
    @test vertex(grid, 1) == P2(0, 0)
    @test vertex(grid, 121) == P2(10, 10)

    grid = CartesianGrid{T}(200, 100)
    if T == Float32
      @test sprint(show, MIME"text/plain"(), grid) ==
            "200×100 CartesianGrid{2,Float32}\n  minimum: Point(0.0f0, 0.0f0)\n  maximum: Point(200.0f0, 100.0f0)\n  spacing: (1.0f0, 1.0f0)"
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), grid) ==
            "200×100 CartesianGrid{2,Float64}\n  minimum: Point(0.0, 0.0)\n  maximum: Point(200.0, 100.0)\n  spacing: (1.0, 1.0)"
    end
  end

  @testset "RectilinearGrid" begin
    x = range(zero(T), stop=one(T), length=6)
    y = T[0.0, 0.1, 0.3, 0.7, 0.9, 1.0]
    grid = RectilinearGrid(x, y)
    @test embeddim(grid) == 2
    @test coordtype(grid) == T
    @test size(grid) == (5, 5)
    @test minimum(grid) == P2(0, 0)
    @test maximum(grid) == P2(1, 1)
    @test extrema(grid) == (P2(0, 0), P2(1, 1))
    @test nelements(grid) == 25
    @test eltype(grid) <: Quadrangle{2,T}
    @test measure(grid) ≈ T(1)
    @test centroid(grid, 1) ≈ P2(0.1, 0.05)
    @test centroid(grid[1]) ≈ P2(0.1, 0.05)
    @test centroid(grid, 2) ≈ P2(0.3, 0.05)
    @test centroid(grid[2]) ≈ P2(0.3, 0.05)

    # single vertex access
    grid = RectilinearGrid(T.(0:10), T.(0:10))
    @test vertex(grid, 1) == P2(0, 0)
    @test vertex(grid, 121) == P2(10, 10)
  end

  @testset "SimpleMesh" begin
    points = P2[(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    mesh = SimpleMesh(points, connec)
    triangles =
      Triangle.([
        P2[(0.0, 0.0), (1.0, 0.0), (0.5, 0.5)],
        P2[(1.0, 0.0), (1.0, 1.0), (0.5, 0.5)],
        P2[(1.0, 1.0), (0.0, 1.0), (0.5, 0.5)],
        P2[(0.0, 1.0), (0.0, 0.0), (0.5, 0.5)]
      ])
    @test vertices(mesh) == points
    @test collect(faces(mesh, 2)) == triangles
    @test collect(elements(mesh)) == triangles
    @test nelements(mesh) == 4
    for i in 1:length(triangles)
      @test mesh[i] == triangles[i]
    end
    @test eltype(mesh) <: Triangle{2,T}
    @test measure(mesh) ≈ T(1)
    @test extrema(mesh) == (P2(0, 0), P2(1, 1))

    # test constructors
    coords = [T.((0, 0)), T.((1, 0)), T.((0, 1)), T.((1, 1)), T.((0.5, 0.5))]
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    mesh = SimpleMesh(coords, SimpleTopology(connec))
    @test eltype(mesh) <: Triangle{2,T}
    @test topology(mesh) isa SimpleTopology
    @test nvertices(mesh) == 5
    @test nelements(mesh) == 4
    mesh = SimpleMesh(coords, connec)
    @test eltype(mesh) <: Triangle{2,T}
    @test topology(mesh) isa SimpleTopology
    @test nvertices(mesh) == 5
    @test nelements(mesh) == 4
    mesh = SimpleMesh(coords, connec, relations=true)
    @test eltype(mesh) <: Triangle{2,T}
    @test topology(mesh) isa HalfEdgeTopology
    @test nvertices(mesh) == 5
    @test nelements(mesh) == 4

    points = P2[(0, 0), (1, 0), (0, 1), (1, 1), (0.25, 0.5), (0.75, 0.5)]
    Δs = connect.([(3, 1, 5), (4, 6, 2)], Triangle)
    □s = connect.([(1, 2, 6, 5), (5, 6, 4, 3)], Quadrangle)
    mesh = SimpleMesh(points, [Δs; □s])
    elms = [
      Triangle(P2[(0.0, 1.0), (0.0, 0.0), (0.25, 0.5)]),
      Triangle(P2[(1.0, 1.0), (0.75, 0.5), (1.0, 0.0)]),
      Quadrangle(P2[(0.0, 0.0), (1.0, 0.0), (0.75, 0.5), (0.25, 0.5)]),
      Quadrangle(P2[(0.25, 0.5), (0.75, 0.5), (1.0, 1.0), (0.0, 1.0)])
    ]
    @test collect(elements(mesh)) == elms
    @test nelements(mesh) == 4
    for i in 1:length(elms)
      @test mesh[i] == elms[i]
    end
    @test eltype(mesh) <: Polygon{2,T}

    # test for https://github.com/JuliaGeometry/Meshes.jl/issues/177
    points = P3[(0, 0, 0), (1, 0, 0), (1, 1, 1), (0, 1, 0)]
    connec = connect.([(1, 2, 3, 4), (3, 4, 1)], [Tetrahedron, Triangle])
    mesh = SimpleMesh(points, connec)
    topo = topology(mesh)
    @test collect(faces(topo, 2)) == [connect((3, 4, 1), Triangle)]
    @test collect(faces(topo, 3)) == [connect((1, 2, 3, 4), Tetrahedron)]

    # test for https://github.com/JuliaGeometry/Meshes.jl/issues/187
    points = P3[(0, 0, 0), (1, 0, 0), (1, 1, 1), (0, 1, 0)]
    connec = connect.([(1, 2, 3, 4), (3, 4, 1)], [Tetrahedron, Triangle])
    mesh = SimpleMesh(points[4:-1:1], connec)
    meshvp = SimpleMesh(view(points, 4:-1:1), connec)
    @test mesh == meshvp

    points = P2[(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    mesh = SimpleMesh(points, connec)
    bytes = @allocated faces(mesh, 2)
    @test bytes < 100
    cells = faces(mesh, 2)
    bytes = @allocated collect(cells)
    @test bytes < 800

    points = P2[(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    mesh = SimpleMesh(points, connec)
    @test centroid(mesh, 1) == centroid(Triangle(P2[(0, 0), (1, 0), (0.5, 0.5)]))
    @test centroid(mesh, 2) == centroid(Triangle(P2[(1, 0), (1, 1), (0.5, 0.5)]))
    @test centroid(mesh, 3) == centroid(Triangle(P2[(1, 1), (0, 1), (0.5, 0.5)]))
    @test centroid(mesh, 4) == centroid(Triangle(P2[(0, 1), (0, 0), (0.5, 0.5)]))

    # merge operation with 2D geometries
    mesh₁ = SimpleMesh(P2[(0, 0), (1, 0), (0, 1)], connect.([(1, 2, 3)]))
    mesh₂ = SimpleMesh(P2[(1, 0), (1, 1), (0, 1)], connect.([(1, 2, 3)]))
    mesh = merge(mesh₁, mesh₂)
    @test vertices(mesh) == [vertices(mesh₁); vertices(mesh₂)]
    @test collect(elements(topology(mesh))) == connect.([(1, 2, 3), (4, 5, 6)])

    # merge operation with 3D geometries
    mesh₁ = SimpleMesh(P3[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)], connect.([(1, 2, 3, 4)], Tetrahedron))
    mesh₂ = SimpleMesh(P3[(1, 0, 0), (1, 1, 0), (0, 1, 0), (1, 1, 1)], connect.([(1, 2, 3, 4)], Tetrahedron))
    mesh = merge(mesh₁, mesh₂)
    @test vertices(mesh) == [vertices(mesh₁); vertices(mesh₂)]
    @test collect(elements(topology(mesh))) == connect.([(1, 2, 3, 4), (5, 6, 7, 8)], Tetrahedron)

    # convert any mesh to SimpleMesh
    grid = CartesianGrid{T}(10, 10)
    mesh = convert(SimpleMesh, grid)
    @test mesh isa SimpleMesh
    @test topology(mesh) == GridTopology(10, 10)
    @test nvertices(mesh) == 121
    @test nelements(mesh) == 100
    @test eltype(mesh) <: Quadrangle

    # test for https://github.com/JuliaGeometry/Meshes.jl/issues/261
    points = rand(P2, 5)
    connec = [connect((1, 2, 3))]
    mesh = SimpleMesh(points, connec)
    @test nvertices(mesh) == length(vertices(mesh)) == 5

    # single vertex access
    points = rand(P2, 5)
    connec = [connect((1, 2, 3))]
    mesh = SimpleMesh(points, connec)
    @test vertex(mesh, 1) == points[1]
    @test vertex(mesh, 2) == points[2]
    @test vertex(mesh, 3) == points[3]
    @test vertex(mesh, 4) == points[4]
    @test vertex(mesh, 5) == points[5]

    points = P2[(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    mesh = SimpleMesh(points, connec)
    if T == Float32
      @test sprint(show, MIME"text/plain"(), mesh) ==
            "4 SimpleMesh{2,Float32}\n  5 vertices\n    └─Point(0.0f0, 0.0f0)\n    └─Point(1.0f0, 0.0f0)\n    └─Point(0.0f0, 1.0f0)\n    └─Point(1.0f0, 1.0f0)\n    └─Point(0.5f0, 0.5f0)\n  4 elements\n    └─Triangle(1, 2, 5)\n    └─Triangle(2, 4, 5)\n    └─Triangle(4, 3, 5)\n    └─Triangle(3, 1, 5)"
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), mesh) ==
            "4 SimpleMesh{2,Float64}\n  5 vertices\n    └─Point(0.0, 0.0)\n    └─Point(1.0, 0.0)\n    └─Point(0.0, 1.0)\n    └─Point(1.0, 1.0)\n    └─Point(0.5, 0.5)\n  4 elements\n    └─Triangle(1, 2, 5)\n    └─Triangle(2, 4, 5)\n    └─Triangle(4, 3, 5)\n    └─Triangle(3, 1, 5)"
    end
  end
end
