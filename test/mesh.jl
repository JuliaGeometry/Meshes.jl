@testset "Meshes" begin
  @testset "CartesianGrid" begin
    @test isgrid(CartesianGrid{1,T})
    @test isgrid(CartesianGrid{2,T})
    @test isgrid(CartesianGrid{3,T})

    grid = CartesianGrid{T}(100)
    @test embeddim(grid) == 1
    @test coordtype(grid) == T
    @test size(grid) == (100,)
    @test minimum(grid) == P1(0)
    @test maximum(grid) == P1(100)
    @test extrema(grid) == (P1(0), P1(100))
    @test spacing(grid) == T[1]
    @test nelements(grid) == 100
    @test eltype(grid) <: Segment{1,T}

    grid = CartesianGrid{T}(200,100)
    @test embeddim(grid) == 2
    @test coordtype(grid) == T
    @test size(grid) == (200, 100)
    @test minimum(grid) == P2(0, 0)
    @test maximum(grid) == P2(200, 100)
    @test extrema(grid) == (P2(0, 0), P2(200, 100))
    @test spacing(grid) == T[1, 1]
    @test nelements(grid) == 200*100
    @test eltype(grid) <: Quadrangle{2,T}

    grid = CartesianGrid((200,100,50), T.((0,0,0)), T.((1,1,1)))
    @test embeddim(grid) == 3
    @test coordtype(grid) == T
    @test size(grid) == (200, 100, 50)
    @test minimum(grid) == P3(0, 0, 0)
    @test maximum(grid) == P3(200, 100, 50)
    @test extrema(grid) == (P3(0, 0, 0), P3(200, 100, 50))
    @test spacing(grid) == T[1, 1, 1]
    @test nelements(grid) == 200*100*50
    @test eltype(grid) <: Hexahedron{3,T}

    grid = CartesianGrid(T.((-1.,-1.)), T.((1.,1.)), dims=(200,100))
    @test embeddim(grid) == 2
    @test coordtype(grid) == T
    @test size(grid) == (200, 100)
    @test minimum(grid) == P2(-1., -1.)
    @test maximum(grid) == P2(1., 1.)
    @test spacing(grid) == T[2/200, 2/100]
    @test nelements(grid) == 200*100
    @test eltype(grid) <: Quadrangle{2,T}

    grid = CartesianGrid((20,10,5), T.((0,0,0)), T.((5,5,5)))
    @test embeddim(grid) == 3
    @test coordtype(grid) == T
    @test size(grid) == (20, 10, 5)
    @test minimum(grid) == P3(0, 0, 0)
    @test maximum(grid) == P3(100, 50, 25)
    @test extrema(grid) == (P3(0, 0, 0), P3(100, 50, 25))
    @test spacing(grid) == T[5, 5, 5]
    @test nelements(grid) == 20*10*5
    @test eltype(grid) <: Hexahedron{3,T}
    @test vertices(grid[1]) == P3[(0, 0, 0), (5, 0, 0), (5, 5, 0), (0, 5, 0), (0, 0, 5), (5, 0, 5), (5, 5, 5), (0, 5, 5)]
    @test all(centroid(grid, i) == centroid(grid[i]) for i in 1:nelements(grid))

    # indexing into a subgrid
    grid = CartesianGrid{T}(10,10)
    sub  = grid[1:2,1:2]
    @test size(sub) == (2,2)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == minimum(grid)
    @test maximum(sub) == P2(2,2)
    sub  = grid[1:1,2:3]
    @test size(sub) == (1,2)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == P2(0,1)
    @test maximum(sub) == P2(1,3)
    sub  = grid[2:4,3:7]
    @test size(sub) == (3,5)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == P2(1,2)
    @test maximum(sub) == P2(4,7)
    grid = CartesianGrid(P2(1,1), P2(11,11), dims=(10,10))
    sub = grid[2:4,3:7]
    @test size(sub) == (3,5)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == P2(2,3)
    @test maximum(sub) == P2(5,8)

    # subgrid from Cartesian ranges
    grid = CartesianGrid{T}(10,10)
    sub1 = grid[1:2,4:6]
    sub2 = grid[CartesianIndex(1,4):CartesianIndex(2,6)]
    @test sub1 == sub2

    grid = CartesianGrid{T}(200,100)
    @test centroid(grid, 1) == P2(0.5, 0.5)
    @test centroid(grid, 2) == P2(1.5, 0.5)
    @test centroid(grid, 200*100) == P2(199.5, 99.5)
    @test nelements(grid) == 200*100
    @test eltype(grid) <: Quadrangle{2,T}
    @test grid[1] == Quadrangle(P2[(0,0), (1,0), (1,1), (0,1)])
    @test grid[2] == Quadrangle(P2[(1,0), (2,0), (2,1), (1,1)])

    grid = CartesianGrid{T}(200,100)
    if T == Float32
      @test sprint(show, MIME"text/plain"(), grid) == "200×100 CartesianGrid{2,Float32}\n  minimum: Point(0.0f0, 0.0f0)\n  maximum: Point(200.0f0, 100.0f0)\n  spacing: (1.0f0, 1.0f0)"
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), grid) == "200×100 CartesianGrid{2,Float64}\n  minimum: Point(0.0, 0.0)\n  maximum: Point(200.0, 100.0)\n  spacing: (1.0, 1.0)"
    end

    if visualtests
      @test_reference "data/grid-1D-$T.png" plot(CartesianGrid{T}(10))
      @test_reference "data/grid-2D-$T.png" plot(CartesianGrid{T}(10,20))
      @test_reference "data/grid-3D-$T.png" plot(CartesianGrid{T}(10,20,30))
      @test_reference "data/grid-1D-$T-data.png" plot(CartesianGrid{T}(10),[1,2,3,4,5,5,4,3,2,1])
      @test_reference "data/grid-2D-$T-data.png" plot(CartesianGrid{T}(10,10),1:100)
      # @test_reference "data/grid3D-data.png" plot(RegularGrid(10,10,10),1:1000)
    end
  end

  @testset "SimpleMesh" begin
    @test !isgrid(SimpleMesh)

    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = SimpleMesh(points, connec)
    triangles = Triangle.([
      P2[(0.0,0.0), (1.0,0.0), (0.5,0.5)],
      P2[(1.0,0.0), (1.0,1.0), (0.5,0.5)],
      P2[(1.0,1.0), (0.0,1.0), (0.5,0.5)],
      P2[(0.0,1.0), (0.0,0.0), (0.5,0.5)]
    ])
    @test vertices(mesh) == points
    @test collect(faces(mesh, 2)) == triangles
    @test collect(elements(mesh)) == triangles
    @test nelements(mesh) == 4
    for i in 1:length(triangles)
      @test mesh[i] == triangles[i]
    end
    @test eltype(mesh) <: Triangle{2,T}

    points = P2[(0,0), (1,0), (0,1), (1,1), (0.25,0.5), (0.75,0.5)]
    Δs = connect.([(3,1,5),(4,6,2)], Triangle)
    □s = connect.([(1,2,6,5),(5,6,4,3)], Quadrangle)
    mesh = SimpleMesh(points, [Δs; □s])
    elms = [
      Triangle(P2[(0.0,1.0), (0.0,0.0), (0.25,0.5)]),
      Triangle(P2[(1.0,1.0), (0.75,0.5), (1.0,0.0)]),
      Quadrangle(P2[(0.0,0.0), (1.0,0.0), (0.75,0.5), (0.25,0.5)]),
      Quadrangle(P2[(0.25,0.5), (0.75,0.5), (1.0,1.0), (0.0,1.0)])
    ]
    @test collect(elements(mesh)) == elms
    @test nelements(mesh) == 4
    for i in 1:length(elms)
      @test mesh[i] == elms[i]
    end
    @test eltype(mesh) <: Polygon{2,T}

    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = SimpleMesh(points, connec)
    bytes = @allocated faces(mesh, 2)
    @test bytes < 100
    cells = faces(mesh, 2)
    bytes = @allocated collect(cells)
    @test bytes < 800

    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = SimpleMesh(points, connec)
    @test centroid(mesh, 1) == centroid(Triangle(P2[(0,0), (1,0), (0.5,0.5)]))
    @test centroid(mesh, 2) == centroid(Triangle(P2[(1,0), (1,1), (0.5,0.5)]))
    @test centroid(mesh, 3) == centroid(Triangle(P2[(1,1), (0,1), (0.5,0.5)]))
    @test centroid(mesh, 4) == centroid(Triangle(P2[(0,1), (0,0), (0.5,0.5)]))

    # merge operation
    mesh₁ = SimpleMesh(P2[(0,0), (1,0), (0,1)], connect.([(1,2,3)]))
    mesh₂ = SimpleMesh(P2[(1,0), (1,1), (0,1)], connect.([(1,2,3)]))
    mesh  = merge(mesh₁, mesh₂)
    @test vertices(mesh) == [vertices(mesh₁); vertices(mesh₂)]
    @test collect(elements(topology(mesh))) == connect.([(1,2,3),(4,5,6)])

    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = SimpleMesh(points, connec)
    if T == Float32
      @test sprint(show, MIME"text/plain"(), mesh) == "4 SimpleMesh{2,Float32}\n  5 vertices\n    └─Point(0.0f0, 0.0f0)\n    └─Point(1.0f0, 0.0f0)\n    └─Point(0.0f0, 1.0f0)\n    └─Point(1.0f0, 1.0f0)\n    └─Point(0.5f0, 0.5f0)\n  4 elements\n    └─Triangle(1, 2, 5)\n    └─Triangle(2, 4, 5)\n    └─Triangle(4, 3, 5)\n    └─Triangle(3, 1, 5)"
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), mesh) == "4 SimpleMesh{2,Float64}\n  5 vertices\n    └─Point(0.0, 0.0)\n    └─Point(1.0, 0.0)\n    └─Point(0.0, 1.0)\n    └─Point(1.0, 1.0)\n    └─Point(0.5, 0.5)\n  4 elements\n    └─Triangle(1, 2, 5)\n    └─Triangle(2, 4, 5)\n    └─Triangle(4, 3, 5)\n    └─Triangle(3, 1, 5)"
    end
  end
end
