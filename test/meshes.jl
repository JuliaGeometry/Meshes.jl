@testset "Meshes ($T)" begin
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
  end

  @testset "UnstructuredMesh" begin
    points = Point2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = UnstructuredMesh(points, connec)
    triangles = Triangle.([
      [(0.0,0.0), (1.0,0.0), (0.5,0.5)],
      [(1.0,0.0), (1.0,1.0), (0.5,0.5)],
      [(1.0,1.0), (0.0,1.0), (0.5,0.5)],
      [(0.0,1.0), (0.0,0.0), (0.5,0.5)]
    ])
    @test vertices(mesh) == points
    bytes = @allocated faces(mesh, 2)
    @test bytes < 100
    cells = faces(mesh, 2)
    bytes = @allocated collect(cells)
    @test bytes < 800
    @test collect(cells) == triangles
    cells = elements(mesh)
    @test collect(cells) == triangles

    points = Point2[(0,0), (1,0), (0,1), (1,1), (0.25,0.5), (0.75,0.5)]
    Δs = connect.([(3,1,5),(4,6,2)], Triangle)
    □s = connect.([(1,2,5,6),(5,6,3,4)], Quadrangle)
    mesh = UnstructuredMesh(points, [Δs; □s])
    elms = [
      Triangle((0.0,1.0), (0.0,0.0), (0.25,0.5)),
      Triangle((1.0,1.0), (0.75,0.5), (1.0,0.0)),
      Quadrangle((0.0,0.0), (1.0,0.0), (0.25,0.5), (0.75,0.5)),
      Quadrangle((0.25,0.5), (0.75,0.5), (0.0,1.0), (1.0,1.0))
    ]
    @test collect(elements(mesh)) == elms
  end
end
