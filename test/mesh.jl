@testset "Meshes" begin
  @testset "CartesianGrid" begin
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

    grid = CartesianGrid{T}(200,100)
    @test coordinates(grid, 1) == T[0.5, 0.5]
    @test coordinates(grid, 2) == T[1.5, 0.5]
    @test coordinates(grid, 200*100) == T[199.5, 99.5]
    @test coordinates(grid, 1:2) == T[0.5 1.5; 0.5 0.5]
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
      @test_ref_plot "data/grid-1D-$T.png" plot(CartesianGrid{T}(10))
      @test_ref_plot "data/grid-2D-$T.png" plot(CartesianGrid{T}(10,20))
      @test_ref_plot "data/grid-3D-$T.png" plot(CartesianGrid{T}(10,20,30))
      @test_ref_plot "data/grid-1D-$T-data.png" plot(CartesianGrid{T}(10),[1,2,3,4,5,5,4,3,2,1])
      @test_ref_plot "data/grid-2D-$T-data.png" plot(CartesianGrid{T}(10,10),1:100)
      # @test_ref_plot "data/grid3D-data.png" plot(RegularGrid(10,10,10),1:1000)
    end
  end

  @testset "UnstructuredMesh" begin
    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = UnstructuredMesh(points, connec)
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
    □s = connect.([(1,2,5,6),(5,6,3,4)], Quadrangle)
    mesh = UnstructuredMesh(points, [Δs; □s])
    elms = [
      Triangle(P2[(0.0,1.0), (0.0,0.0), (0.25,0.5)]),
      Triangle(P2[(1.0,1.0), (0.75,0.5), (1.0,0.0)]),
      Quadrangle(P2[(0.0,0.0), (1.0,0.0), (0.25,0.5), (0.75,0.5)]),
      Quadrangle(P2[(0.25,0.5), (0.75,0.5), (0.0,1.0), (1.0,1.0)])
    ]
    @test collect(elements(mesh)) == elms
    @test nelements(mesh) == 4
    for i in 1:length(elms)
      @test mesh[i] == elms[i]
    end
    @test eltype(mesh) <: Polygon{2,T}

    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = UnstructuredMesh(points, connec)
    bytes = @allocated faces(mesh, 2)
    @test bytes < 100
    cells = faces(mesh, 2)
    bytes = @allocated collect(cells)
    @test bytes < 800

    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = UnstructuredMesh(points, connec)
    @test coordinates(mesh, 1) == coordinates(centroid(Triangle(P2[(0,0), (1,0), (0.5,0.5)])))
    @test coordinates(mesh, 2) == coordinates(centroid(Triangle(P2[(1,0), (1,1), (0.5,0.5)])))
    @test coordinates(mesh, 3) == coordinates(centroid(Triangle(P2[(1,1), (0,1), (0.5,0.5)])))
    @test coordinates(mesh, 4) == coordinates(centroid(Triangle(P2[(0,1), (0,0), (0.5,0.5)])))

    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = UnstructuredMesh(points, connec)
    if T == Float32
      @test sprint(show, MIME"text/plain"(), mesh) == "4 UnstructuredMesh{2,Float32}\n  5 vertices\n    └─Point(0.0f0, 0.0f0)\n    └─Point(1.0f0, 0.0f0)\n    └─Point(0.0f0, 1.0f0)\n    └─Point(1.0f0, 1.0f0)\n    └─Point(0.5f0, 0.5f0)\n  4 faces\n    └─Triangle(1, 2, 5)\n    └─Triangle(2, 4, 5)\n    └─Triangle(4, 3, 5)\n    └─Triangle(3, 1, 5)"
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), mesh) == "4 UnstructuredMesh{2,Float64}\n  5 vertices\n    └─Point(0.0, 0.0)\n    └─Point(1.0, 0.0)\n    └─Point(0.0, 1.0)\n    └─Point(1.0, 1.0)\n    └─Point(0.5, 0.5)\n  4 faces\n    └─Triangle(1, 2, 5)\n    └─Triangle(2, 4, 5)\n    └─Triangle(4, 3, 5)\n    └─Triangle(3, 1, 5)"
    end
  end
end
