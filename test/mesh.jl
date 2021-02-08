@testset "Meshes" begin
  @testset "CartesianGrid" begin
    grid = CartesianGrid{T}(200,100)
    @test embeddim(grid) == 2
    @test coordtype(grid) == T
    @test size(grid) == (200, 100)
    @test minimum(grid) == P2(0, 0)
    @test maximum(grid) == P2(200, 100)
    @test spacing(grid) == T[1, 1]

    grid = CartesianGrid((200,100,50), T.((0,0,0)), T.((1,1,1)))
    @test embeddim(grid) == 3
    @test coordtype(grid) == T
    @test size(grid) == (200, 100, 50)
    @test minimum(grid) == P3(0, 0, 0)
    @test maximum(grid) == P3(200, 100, 50)
    @test spacing(grid) == T[1, 1, 1]

    grid = CartesianGrid(T.((-1.,-1.)), T.((1.,1.)), dims=(200,100))
    @test embeddim(grid) == 2
    @test coordtype(grid) == T
    @test size(grid) == (200, 100)
    @test minimum(grid) == P2(-1., -1.)
    @test maximum(grid) == P2(1., 1.)
    @test spacing(grid) == T[2/200, 2/100]

    grid = CartesianGrid{T}(200,100)
    @test coordinates(grid, 1) == T[0.5, 0.5]
    @test coordinates(grid, 2) == T[1.5, 0.5]
    @test coordinates(grid, 200*100) == T[199.5, 99.5]
    @test coordinates(grid, 1:2) == T[0.5 1.5; 0.5 0.5]

    grid = CartesianGrid{T}(200,100)
    if T == Float32
      @test sprint(show, MIME"text/plain"(), grid) == "200×100 CartesianGrid{2,Float32}\n  minimum: Point(0.0f0, 0.0f0)\n  maximum: Point(200.0f0, 100.0f0)\n  spacing: (1.0f0, 1.0f0)"
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), grid) == "200×100 CartesianGrid{2,Float64}\n  minimum: Point(0.0, 0.0)\n  maximum: Point(200.0, 100.0)\n  spacing: (1.0, 1.0)"
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
    @test coordinates(mesh, 1) == coordinates(center(Triangle(P2[(0,0), (1,0), (0.5,0.5)])))
    @test coordinates(mesh, 2) == coordinates(center(Triangle(P2[(1,0), (1,1), (0.5,0.5)])))
    @test coordinates(mesh, 3) == coordinates(center(Triangle(P2[(1,1), (0,1), (0.5,0.5)])))
    @test coordinates(mesh, 4) == coordinates(center(Triangle(P2[(0,1), (0,0), (0.5,0.5)])))

    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = UnstructuredMesh(points, connec)
    if T == Float32
      @test sprint(show, MIME"text/plain"(), mesh) == "UnstructuredMesh{2,Float32}\n  5 vertices\n    └─Point(0.0f0, 0.0f0)\n    └─Point(1.0f0, 0.0f0)\n    └─Point(0.0f0, 1.0f0)\n    └─Point(1.0f0, 1.0f0)\n    └─Point(0.5f0, 0.5f0)\n  4 faces\n    └─Triangle(1, 2, 5)\n    └─Triangle(2, 4, 5)\n    └─Triangle(4, 3, 5)\n    └─Triangle(3, 1, 5)"
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), mesh) == "UnstructuredMesh{2,Float64}\n  5 vertices\n    └─Point(0.0, 0.0)\n    └─Point(1.0, 0.0)\n    └─Point(0.0, 1.0)\n    └─Point(1.0, 1.0)\n    └─Point(0.5, 0.5)\n  4 faces\n    └─Triangle(1, 2, 5)\n    └─Triangle(2, 4, 5)\n    └─Triangle(4, 3, 5)\n    └─Triangle(3, 1, 5)"
    end
  end
end
