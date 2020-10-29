@testset "Meshes" begin
  @testset "CartesianGrid" begin
    grid = CartesianGrid{Float32}(200,100)
    @test embeddim(grid) == 2
    @test coordtype(grid) == Float32
    @test size(grid) == (200,100)
    @test minimum(grid) == Point(0f0, 0f0)
    @test maximum(grid) == Point(199f0, 99f0)
    @test spacing(grid) == [1f0, 1f0]

    grid = CartesianGrid((200,100,50), (0.,0.,0.), (1.,1.,1.))
    @test embeddim(grid) == 3
    @test coordtype(grid) == Float64
    @test size(grid) == (200,100,50)
    @test minimum(grid) == Point(0., 0., 0.)
    @test maximum(grid) == Point(199., 99., 49.)
    @test spacing(grid) == [1.,1.,1.]

    grid = CartesianGrid((-1.,-1.), (1.,1.), dims=(200,100))
    @test embeddim(grid) == 2
    @test coordtype(grid) == Float64
    @test size(grid) == (200,100)
    @test minimum(grid) == Point(-1., -1.)
    @test maximum(grid) == Point(1., 1.)
    @test spacing(grid) == [2/199, 2/99]
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
