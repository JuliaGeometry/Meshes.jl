@testset "Meshes" begin
  @testset "CartesianGrid" begin
    grid = cartgrid(100)
    @test embeddim(grid) == 1
    @test Meshes.crs(grid) <: Cartesian{NoDatum}
    @test Meshes.lentype(grid) == ℳ
    @test size(grid) == (100,)
    @test minimum(grid) == cart(0)
    @test maximum(grid) == cart(100)
    @test extrema(grid) == (cart(0), cart(100))
    @test spacing(grid) == (T(1) * u"m",)
    @test nelements(grid) == 100
    @test eltype(grid) <: Segment{1}
    @test measure(grid) ≈ T(100) * u"m"
    @test vertex(grid, 1) == vertex(grid, ntuple(i -> 1, embeddim(grid)))
    @test vertex(grid, nvertices(grid)) == vertex(grid, size(grid) .+ 1)
    @test grid[1] == Segment(cart(0), cart(1))
    @test grid[100] == Segment(cart(99), cart(100))

    grid = cartgrid(200, 100)
    @test embeddim(grid) == 2
    @test Meshes.crs(grid) <: Cartesian{NoDatum}
    @test Meshes.lentype(grid) == ℳ
    @test size(grid) == (200, 100)
    @test minimum(grid) == cart(0, 0)
    @test maximum(grid) == cart(200, 100)
    @test extrema(grid) == (cart(0, 0), cart(200, 100))
    @test spacing(grid) == (T(1) * u"m", T(1) * u"m")
    @test nelements(grid) == 200 * 100
    @test eltype(grid) <: Quadrangle{2}
    @test measure(grid) ≈ T(200 * 100) * u"m^2"
    @test vertex(grid, 1) == vertex(grid, ntuple(i -> 1, embeddim(grid)))
    @test vertex(grid, nvertices(grid)) == vertex(grid, size(grid) .+ 1)
    @test grid[1, 1] == grid[1]
    @test grid[200, 100] == grid[20000]

    grid = CartesianGrid((200, 100, 50), T.((0, 0, 0)), T.((1, 1, 1)))
    @test embeddim(grid) == 3
    @test Meshes.crs(grid) <: Cartesian{NoDatum}
    @test Meshes.lentype(grid) == ℳ
    @test size(grid) == (200, 100, 50)
    @test minimum(grid) == cart(0, 0, 0)
    @test maximum(grid) == cart(200, 100, 50)
    @test extrema(grid) == (cart(0, 0, 0), cart(200, 100, 50))
    @test spacing(grid) == (T(1) * u"m", T(1) * u"m", T(1) * u"m")
    @test nelements(grid) == 200 * 100 * 50
    @test eltype(grid) <: Hexahedron{3}
    @test measure(grid) ≈ T(200 * 100 * 50) * u"m^3"
    @test vertex(grid, 1) == vertex(grid, ntuple(i -> 1, embeddim(grid)))
    @test vertex(grid, nvertices(grid)) == vertex(grid, size(grid) .+ 1)
    @test grid[1, 1, 1] == grid[1]
    @test grid[200, 100, 50] == grid[1000000]

    grid = CartesianGrid(T.((0, 0, 0)), T.((1, 1, 1)), T.((0.1, 0.1, 0.1)))
    @test embeddim(grid) == 3
    @test Meshes.crs(grid) <: Cartesian{NoDatum}
    @test Meshes.lentype(grid) == ℳ
    @test size(grid) == (10, 10, 10)
    @test minimum(grid) == cart(0, 0, 0)
    @test maximum(grid) == cart(1, 1, 1)
    @test spacing(grid) == (T(0.1) * u"m", T(0.1) * u"m", T(0.1) * u"m")

    grid = CartesianGrid(T.((-1.0, -1.0)), T.((1.0, 1.0)), dims=(200, 100))
    @test embeddim(grid) == 2
    @test Meshes.crs(grid) <: Cartesian{NoDatum}
    @test Meshes.lentype(grid) == ℳ
    @test size(grid) == (200, 100)
    @test minimum(grid) == cart(-1.0, -1.0)
    @test maximum(grid) == cart(1.0, 1.0)
    @test spacing(grid) == (T(2 / 200) * u"m", T(2 / 100) * u"m")
    @test nelements(grid) == 200 * 100
    @test eltype(grid) <: Quadrangle{2}

    grid = CartesianGrid((20, 10, 5), T.((0, 0, 0)), T.((5, 5, 5)))
    @test embeddim(grid) == 3
    @test Meshes.crs(grid) <: Cartesian{NoDatum}
    @test Meshes.lentype(grid) == ℳ
    @test size(grid) == (20, 10, 5)
    @test minimum(grid) == cart(0, 0, 0)
    @test maximum(grid) == cart(100, 50, 25)
    @test extrema(grid) == (cart(0, 0, 0), cart(100, 50, 25))
    @test spacing(grid) == (T(5) * u"m", T(5) * u"m", T(5) * u"m")
    @test nelements(grid) == 20 * 10 * 5
    @test eltype(grid) <: Hexahedron{3}
    @test vertices(grid[1]) == (
      cart(0, 0, 0),
      cart(5, 0, 0),
      cart(5, 5, 0),
      cart(0, 5, 0),
      cart(0, 0, 5),
      cart(5, 0, 5),
      cart(5, 5, 5),
      cart(0, 5, 5)
    )
    @test all(centroid(grid, i) == centroid(grid[i]) for i in 1:nelements(grid))

    # constructor with offset
    grid = CartesianGrid((10, 10), T.((1.0, 1.0)), T.((1.0, 1.0)), (2, 2))
    @test embeddim(grid) == 2
    @test Meshes.crs(grid) <: Cartesian{NoDatum}
    @test Meshes.lentype(grid) == ℳ
    @test size(grid) == (10, 10)
    @test minimum(grid) == cart(0.0, 0.0)
    @test maximum(grid) == cart(10.0, 10.0)
    @test spacing(grid) == (T(1) * u"m", T(1) * u"m")
    @test nelements(grid) == 10 * 10
    @test eltype(grid) <: Quadrangle{2}

    # mixed units
    grid = CartesianGrid((10, 10), (T(0) * u"m", T(0) * u"cm"), (T(100) * u"cm", T(1) * u"m"))
    @test unit(Meshes.lentype(grid)) == u"m"
    grid = CartesianGrid((T(0) * u"cm", T(0) * u"m"), (T(10) * u"m", T(1000) * u"cm"), (T(100) * u"cm", T(1) * u"m"))
    @test unit(Meshes.lentype(grid)) == u"m"
    grid = CartesianGrid((T(0) * u"cm", T(0) * u"m"), (T(10) * u"m", T(1000) * u"cm"), dims=(10, 10))
    @test unit(Meshes.lentype(grid)) == u"m"

    # indexing into a subgrid
    grid = cartgrid(10, 10)
    sub = grid[1:2, 1:2]
    @test size(sub) == (2, 2)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == minimum(grid)
    @test maximum(sub) == cart(2, 2)
    sub = grid[1:1, 2:3]
    @test size(sub) == (1, 2)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == cart(0, 1)
    @test maximum(sub) == cart(1, 3)
    sub = grid[2:4, 3:7]
    @test size(sub) == (3, 5)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == cart(1, 2)
    @test maximum(sub) == cart(4, 7)
    grid = CartesianGrid(cart(1, 1), cart(11, 11), dims=(10, 10))
    sub = grid[2:4, 3:7]
    @test size(sub) == (3, 5)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == cart(2, 3)
    @test maximum(sub) == cart(5, 8)
    sub = grid[2, 3:7]
    @test size(sub) == (1, 5)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == cart(2, 3)
    @test maximum(sub) == cart(3, 8)
    sub = grid[:, 3:7]
    @test size(sub) == (10, 5)
    @test spacing(sub) == spacing(grid)
    @test minimum(sub) == cart(1, 3)
    @test maximum(sub) == cart(11, 8)
    @test_throws BoundsError grid[3:11, :]

    # subgrid with comparable vertices of grid
    grid = CartesianGrid((10, 10), cart(0.0, 0.0), T.((1.2, 1.2)))
    sub = grid[2:4, 5:7]
    @test sub == CartesianGrid((3, 3), cart(0.0, 0.0), T.((1.2, 1.2)), (0, -3))
    ind = reshape(reshape(1:121, 11, 11)[2:5, 5:8], :)
    @test vertices(grid)[ind] == vertices(sub)

    # subgrid from Cartesian ranges
    grid = cartgrid(10, 10)
    sub1 = grid[1:2, 4:6]
    sub2 = grid[CartesianIndex(1, 4):CartesianIndex(2, 6)]
    @test sub1 == sub2

    grid = cartgrid(200, 100)
    @test centroid(grid, 1) == cart(0.5, 0.5)
    @test centroid(grid, 2) == cart(1.5, 0.5)
    @test centroid(grid, 200 * 100) == cart(199.5, 99.5)
    @test nelements(grid) == 200 * 100
    @test eltype(grid) <: Quadrangle{2}
    @test grid[1] == Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    @test grid[2] == Quadrangle(cart(1, 0), cart(2, 0), cart(2, 1), cart(1, 1))

    # expand CartesianGrid with comparable vertices
    grid = CartesianGrid((10, 10), cart(0.0, 0.0), T.((1.0, 1.0)))
    left, right = (1, 1), (1, 1)
    newdim = size(grid) .+ left .+ right
    newoffset = offset(grid) .+ left
    grid2 = CartesianGrid(newdim, minimum(grid), spacing(grid), newoffset)
    @test issubset(vertices(grid), vertices(grid2))

    # GridTopology from CartesianGrid
    grid = cartgrid(5, 5)
    topo = topology(grid)
    vs = vertices(grid)
    for i in 1:nelements(grid)
      inds = indices(element(topo, i))
      @test vs[[inds...]] == pointify(element(grid, i))
    end

    # convert topology
    grid = cartgrid(10, 10)
    mesh = topoconvert(HalfEdgeTopology, grid)
    @test mesh isa SimpleMesh
    @test nvertices(mesh) == 121
    @test nelements(mesh) == 100
    @test eltype(mesh) <: Quadrangle

    # single vertex access
    grid = cartgrid(10, 10)
    @test vertex(grid, 1) == cart(0, 0)
    @test vertex(grid, 121) == cart(10, 10)

    # xyz
    g1D = cartgrid(10)
    g2D = cartgrid(10, 10)
    g3D = cartgrid(10, 10, 10)
    @test Meshes.xyz(g1D) == (T.(0:10) * u"m",)
    @test Meshes.xyz(g2D) == (T.(0:10) * u"m", T.(0:10) * u"m")
    @test Meshes.xyz(g3D) == (T.(0:10) * u"m", T.(0:10) * u"m", T.(0:10) * u"m")

    # XYZ
    g1D = cartgrid(10)
    g2D = cartgrid(10, 10)
    g3D = cartgrid(10, 10, 10)
    x = T.(0:10) * u"m"
    y = T.(0:10)' * u"m"
    z = reshape(T.(0:10), 1, 1, 11) * u"m"
    @test Meshes.XYZ(g1D) == (x,)
    @test Meshes.XYZ(g2D) == (repeat(x, 1, 11), repeat(y, 11, 1))
    @test Meshes.XYZ(g3D) == (repeat(x, 1, 11, 11), repeat(y, 11, 1, 11), repeat(z, 11, 11, 1))

    # units
    grid = CartesianGrid((10, 10), cart(0, 0), (T(1) * u"m", T(1) * u"m"))
    o = minimum(grid)
    s = spacing(grid)
    @test unit(Meshes.lentype(o)) == u"m"
    @test Unitful.numtype(Meshes.lentype(o)) === T
    @test unit(eltype(s)) == u"m"
    @test Unitful.numtype(eltype(s)) === T

    grid = cartgrid(200, 100)
    if T == Float32
      @test sprint(show, MIME"text/plain"(), grid) == """
      200×100 CartesianGrid
      ├─ minimum: Point(x: 0.0f0 m, y: 0.0f0 m)
      ├─ maximum: Point(x: 200.0f0 m, y: 100.0f0 m)
      └─ spacing: (1.0f0 m, 1.0f0 m)"""
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), grid) == """
      200×100 CartesianGrid
      ├─ minimum: Point(x: 0.0 m, y: 0.0 m)
      ├─ maximum: Point(x: 200.0 m, y: 100.0 m)
      └─ spacing: (1.0 m, 1.0 m)"""
    end
  end

  @testset "RectilinearGrid" begin
    x = range(zero(T), stop=one(T), length=6)
    y = T[0.0, 0.1, 0.3, 0.7, 0.9, 1.0]
    grid = RectilinearGrid(x, y)
    @test embeddim(grid) == 2
    @test Meshes.crs(grid) <: Cartesian{NoDatum}
    @test Meshes.lentype(grid) == ℳ
    @test size(grid) == (5, 5)
    @test minimum(grid) == cart(0, 0)
    @test maximum(grid) == cart(1, 1)
    @test extrema(grid) == (cart(0, 0), cart(1, 1))
    @test nelements(grid) == 25
    @test eltype(grid) <: Quadrangle{2}
    @test measure(grid) ≈ T(1) * u"m^2"
    @test centroid(grid, 1) ≈ cart(0.1, 0.05)
    @test centroid(grid[1]) ≈ cart(0.1, 0.05)
    @test centroid(grid, 2) ≈ cart(0.3, 0.05)
    @test centroid(grid[2]) ≈ cart(0.3, 0.05)
    @test vertex(grid, 1) == vertex(grid, ntuple(i -> 1, embeddim(grid)))
    @test vertex(grid, nvertices(grid)) == vertex(grid, size(grid) .+ 1)
    @test grid[1, 1] == grid[1]
    @test grid[5, 5] == grid[25]
    sub = grid[2:4, 3:5]
    @test size(sub) == (3, 3)
    @test minimum(sub) == cart(0.2, 0.3)
    @test maximum(sub) == cart(0.8, 1.0)
    sub = grid[2, 3:5]
    @test size(sub) == (1, 3)
    @test minimum(sub) == cart(0.2, 0.3)
    @test maximum(sub) == cart(0.4, 1.0)
    sub = grid[:, 3:5]
    @test size(sub) == (5, 3)
    @test minimum(sub) == cart(0.0, 0.3)
    @test maximum(sub) == cart(1.0, 1.0)
    @test_throws BoundsError grid[2:6, :]
    @test Meshes.xyz(grid) == (x * u"m", y * u"m")
    @test Meshes.XYZ(grid) == (repeat(x, 1, 6) * u"m", repeat(y', 6, 1) * u"m")

    # single vertex access
    grid = RectilinearGrid(T.(0:10), T.(0:10))
    @test vertex(grid, 1) == cart(0, 0)
    @test vertex(grid, 121) == cart(10, 10)

    # constructor with datum & datum propagation
    grid = RectilinearGrid{WGS84Latest}(x, y)
    @test datum(Meshes.crs(grid)) === WGS84Latest
    @test datum(Meshes.crs(centroid(grid))) === WGS84Latest

    # conversion
    cg = cartgrid(10, 10)
    rg = convert(RectilinearGrid, cg)
    @test size(rg) == size(cg)
    @test nvertices(rg) == nvertices(cg)
    @test nelements(rg) == nelements(cg)
    @test topology(rg) == topology(cg)
    @test vertices(rg) == vertices(cg)

    cg = cartgrid(10, 20, 30)
    rg = convert(RectilinearGrid, cg)
    @test size(rg) == size(cg)
    @test nvertices(rg) == nvertices(cg)
    @test nelements(rg) == nelements(cg)
    @test topology(rg) == topology(cg)
    @test vertices(rg) == vertices(cg)

    x = range(zero(T), stop=one(T), length=6)
    y = T[0.0, 0.1, 0.3, 0.7, 0.9, 1.0]
    grid = RectilinearGrid(x, y)
    @test sprint(show, grid) == "5×5 RectilinearGrid"
    if T == Float32
      @test sprint(show, MIME"text/plain"(), grid) == """
      5×5 RectilinearGrid
        36 vertices
        ├─ Point(x: 0.0f0 m, y: 0.0f0 m)
        ├─ Point(x: 0.2f0 m, y: 0.0f0 m)
        ├─ Point(x: 0.4f0 m, y: 0.0f0 m)
        ├─ Point(x: 0.6f0 m, y: 0.0f0 m)
        ├─ Point(x: 0.8f0 m, y: 0.0f0 m)
        ⋮
        ├─ Point(x: 0.2f0 m, y: 1.0f0 m)
        ├─ Point(x: 0.4f0 m, y: 1.0f0 m)
        ├─ Point(x: 0.6f0 m, y: 1.0f0 m)
        ├─ Point(x: 0.8f0 m, y: 1.0f0 m)
        └─ Point(x: 1.0f0 m, y: 1.0f0 m)
        25 elements
        ├─ Quadrangle(1, 2, 8, 7)
        ├─ Quadrangle(2, 3, 9, 8)
        ├─ Quadrangle(3, 4, 10, 9)
        ├─ Quadrangle(4, 5, 11, 10)
        ├─ Quadrangle(5, 6, 12, 11)
        ⋮
        ├─ Quadrangle(25, 26, 32, 31)
        ├─ Quadrangle(26, 27, 33, 32)
        ├─ Quadrangle(27, 28, 34, 33)
        ├─ Quadrangle(28, 29, 35, 34)
        └─ Quadrangle(29, 30, 36, 35)"""
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), grid) == """
      5×5 RectilinearGrid
        36 vertices
        ├─ Point(x: 0.0 m, y: 0.0 m)
        ├─ Point(x: 0.2 m, y: 0.0 m)
        ├─ Point(x: 0.4 m, y: 0.0 m)
        ├─ Point(x: 0.6 m, y: 0.0 m)
        ├─ Point(x: 0.8 m, y: 0.0 m)
        ⋮
        ├─ Point(x: 0.2 m, y: 1.0 m)
        ├─ Point(x: 0.4 m, y: 1.0 m)
        ├─ Point(x: 0.6 m, y: 1.0 m)
        ├─ Point(x: 0.8 m, y: 1.0 m)
        └─ Point(x: 1.0 m, y: 1.0 m)
        25 elements
        ├─ Quadrangle(1, 2, 8, 7)
        ├─ Quadrangle(2, 3, 9, 8)
        ├─ Quadrangle(3, 4, 10, 9)
        ├─ Quadrangle(4, 5, 11, 10)
        ├─ Quadrangle(5, 6, 12, 11)
        ⋮
        ├─ Quadrangle(25, 26, 32, 31)
        ├─ Quadrangle(26, 27, 33, 32)
        ├─ Quadrangle(27, 28, 34, 33)
        ├─ Quadrangle(28, 29, 35, 34)
        └─ Quadrangle(29, 30, 36, 35)"""
    end
  end

  @testset "StructuredGrid" begin
    X = repeat(range(zero(T), stop=one(T), length=6), 1, 6)
    Y = repeat(T[0.0, 0.1, 0.3, 0.7, 0.9, 1.0]', 6, 1)
    grid = StructuredGrid(X, Y)
    @test embeddim(grid) == 2
    @test Meshes.crs(grid) <: Cartesian{NoDatum}
    @test Meshes.lentype(grid) == ℳ
    @test size(grid) == (5, 5)
    @test minimum(grid) == cart(0, 0)
    @test maximum(grid) == cart(1, 1)
    @test extrema(grid) == (cart(0, 0), cart(1, 1))
    @test nelements(grid) == 25
    @test eltype(grid) <: Quadrangle{2}
    @test measure(grid) ≈ T(1) * u"m^2"
    @test centroid(grid, 1) ≈ cart(0.1, 0.05)
    @test centroid(grid[1]) ≈ cart(0.1, 0.05)
    @test centroid(grid, 2) ≈ cart(0.3, 0.05)
    @test centroid(grid[2]) ≈ cart(0.3, 0.05)
    @test vertex(grid, 1) == vertex(grid, ntuple(i -> 1, embeddim(grid)))
    @test vertex(grid, nvertices(grid)) == vertex(grid, size(grid) .+ 1)
    @test grid[1, 1] == grid[1]
    @test grid[5, 5] == grid[25]
    sub = grid[2:4, 3:5]
    @test size(sub) == (3, 3)
    @test minimum(sub) == cart(0.2, 0.3)
    @test maximum(sub) == cart(0.8, 1.0)
    sub = grid[2, 3:5]
    @test size(sub) == (1, 3)
    @test minimum(sub) == cart(0.2, 0.3)
    @test maximum(sub) == cart(0.4, 1.0)
    sub = grid[:, 3:5]
    @test size(sub) == (5, 3)
    @test minimum(sub) == cart(0.0, 0.3)
    @test maximum(sub) == cart(1.0, 1.0)
    @test_throws BoundsError grid[2:6, :]
    @test Meshes.XYZ(grid) == (X * u"m", Y * u"m")

    # constructor with datum
    grid = StructuredGrid{WGS84Latest}(X, Y)
    @test datum(Meshes.crs(grid)) === WGS84Latest

    # conversion
    cg = cartgrid(10, 10)
    sg = convert(StructuredGrid, cg)
    @test size(sg) == size(cg)
    @test nvertices(sg) == nvertices(cg)
    @test nelements(sg) == nelements(cg)
    @test topology(sg) == topology(cg)
    @test vertices(sg) == vertices(cg)

    cg = cartgrid(10, 20, 30)
    sg = convert(StructuredGrid, cg)
    @test size(sg) == size(cg)
    @test nvertices(sg) == nvertices(cg)
    @test nelements(sg) == nelements(cg)
    @test topology(sg) == topology(cg)
    @test vertices(sg) == vertices(cg)

    rg = RectilinearGrid(T.(0:10), T.(0:10))
    sg = convert(StructuredGrid, rg)
    @test size(sg) == size(rg)
    @test nvertices(sg) == nvertices(rg)
    @test nelements(sg) == nelements(rg)
    @test topology(sg) == topology(rg)
    @test vertices(sg) == vertices(rg)

    rg = RectilinearGrid(T.(0:10), T.(0:20), T.(0:30))
    sg = convert(StructuredGrid, rg)
    @test size(sg) == size(rg)
    @test nvertices(sg) == nvertices(rg)
    @test nelements(sg) == nelements(rg)
    @test topology(sg) == topology(rg)
    @test vertices(sg) == vertices(rg)

    X = repeat(range(zero(T), stop=one(T), length=6), 1, 6)
    Y = repeat(T[0.0, 0.1, 0.3, 0.7, 0.9, 1.0]', 6, 1)
    grid = StructuredGrid(X, Y)
    @test sprint(show, grid) == "5×5 StructuredGrid"
    if T == Float32
      @test sprint(show, MIME"text/plain"(), grid) == """
      5×5 StructuredGrid
        36 vertices
        ├─ Point(x: 0.0f0 m, y: 0.0f0 m)
        ├─ Point(x: 0.2f0 m, y: 0.0f0 m)
        ├─ Point(x: 0.4f0 m, y: 0.0f0 m)
        ├─ Point(x: 0.6f0 m, y: 0.0f0 m)
        ├─ Point(x: 0.8f0 m, y: 0.0f0 m)
        ⋮
        ├─ Point(x: 0.2f0 m, y: 1.0f0 m)
        ├─ Point(x: 0.4f0 m, y: 1.0f0 m)
        ├─ Point(x: 0.6f0 m, y: 1.0f0 m)
        ├─ Point(x: 0.8f0 m, y: 1.0f0 m)
        └─ Point(x: 1.0f0 m, y: 1.0f0 m)
        25 elements
        ├─ Quadrangle(1, 2, 8, 7)
        ├─ Quadrangle(2, 3, 9, 8)
        ├─ Quadrangle(3, 4, 10, 9)
        ├─ Quadrangle(4, 5, 11, 10)
        ├─ Quadrangle(5, 6, 12, 11)
        ⋮
        ├─ Quadrangle(25, 26, 32, 31)
        ├─ Quadrangle(26, 27, 33, 32)
        ├─ Quadrangle(27, 28, 34, 33)
        ├─ Quadrangle(28, 29, 35, 34)
        └─ Quadrangle(29, 30, 36, 35)"""
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), grid) == """
      5×5 StructuredGrid
        36 vertices
        ├─ Point(x: 0.0 m, y: 0.0 m)
        ├─ Point(x: 0.2 m, y: 0.0 m)
        ├─ Point(x: 0.4 m, y: 0.0 m)
        ├─ Point(x: 0.6 m, y: 0.0 m)
        ├─ Point(x: 0.8 m, y: 0.0 m)
        ⋮
        ├─ Point(x: 0.2 m, y: 1.0 m)
        ├─ Point(x: 0.4 m, y: 1.0 m)
        ├─ Point(x: 0.6 m, y: 1.0 m)
        ├─ Point(x: 0.8 m, y: 1.0 m)
        └─ Point(x: 1.0 m, y: 1.0 m)
        25 elements
        ├─ Quadrangle(1, 2, 8, 7)
        ├─ Quadrangle(2, 3, 9, 8)
        ├─ Quadrangle(3, 4, 10, 9)
        ├─ Quadrangle(4, 5, 11, 10)
        ├─ Quadrangle(5, 6, 12, 11)
        ⋮
        ├─ Quadrangle(25, 26, 32, 31)
        ├─ Quadrangle(26, 27, 33, 32)
        ├─ Quadrangle(27, 28, 34, 33)
        ├─ Quadrangle(28, 29, 35, 34)
        └─ Quadrangle(29, 30, 36, 35)"""
    end
  end

  @testset "SimpleMesh" begin
    points = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    mesh = SimpleMesh(points, connec)
    triangles =
      Triangle.([
        (cart(0.0, 0.0), cart(1.0, 0.0), cart(0.5, 0.5)),
        (cart(1.0, 0.0), cart(1.0, 1.0), cart(0.5, 0.5)),
        (cart(1.0, 1.0), cart(0.0, 1.0), cart(0.5, 0.5)),
        (cart(0.0, 1.0), cart(0.0, 0.0), cart(0.5, 0.5))
      ])
    @test Meshes.crs(mesh) <: Cartesian{NoDatum}
    @test Meshes.lentype(mesh) == ℳ
    @test vertices(mesh) == points
    @test collect(faces(mesh, 2)) == triangles
    @test collect(elements(mesh)) == triangles
    @test nelements(mesh) == 4
    for i in 1:length(triangles)
      @test mesh[i] == triangles[i]
    end
    @test eltype(mesh) <: Triangle{2}
    @test measure(mesh) ≈ T(1) * u"m^2"
    @test area(mesh) ≈ T(1) * u"m^2"
    @test extrema(mesh) == (cart(0, 0), cart(1, 1))

    # test constructors
    coords = [T.((0, 0)), T.((1, 0)), T.((0, 1)), T.((1, 1)), T.((0.5, 0.5))]
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    mesh = SimpleMesh(coords, SimpleTopology(connec))
    @test eltype(mesh) <: Triangle{2}
    @test topology(mesh) isa SimpleTopology
    @test nvertices(mesh) == 5
    @test nelements(mesh) == 4
    mesh = SimpleMesh(coords, connec)
    @test eltype(mesh) <: Triangle{2}
    @test topology(mesh) isa SimpleTopology
    @test nvertices(mesh) == 5
    @test nelements(mesh) == 4
    mesh = SimpleMesh(coords, connec, relations=true)
    @test eltype(mesh) <: Triangle{2}
    @test topology(mesh) isa HalfEdgeTopology
    @test nvertices(mesh) == 5
    @test nelements(mesh) == 4

    points = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.25, 0.5), (0.75, 0.5)])
    Δs = connect.([(3, 1, 5), (4, 6, 2)], Triangle)
    □s = connect.([(1, 2, 6, 5), (5, 6, 4, 3)], Quadrangle)
    mesh = SimpleMesh(points, [Δs; □s])
    elms = [
      Triangle(cart(0.0, 1.0), cart(0.0, 0.0), cart(0.25, 0.5)),
      Triangle(cart(1.0, 1.0), cart(0.75, 0.5), cart(1.0, 0.0)),
      Quadrangle(cart(0.0, 0.0), cart(1.0, 0.0), cart(0.75, 0.5), cart(0.25, 0.5)),
      Quadrangle(cart(0.25, 0.5), cart(0.75, 0.5), cart(1.0, 1.0), cart(0.0, 1.0))
    ]
    @test collect(elements(mesh)) == elms
    @test nelements(mesh) == 4
    for i in 1:length(elms)
      @test mesh[i] == elms[i]
    end
    @test eltype(mesh) <: Polygon{2}

    # test for https://github.com/JuliaGeometry/Meshes.jl/issues/177
    points = cart.([(0, 0, 0), (1, 0, 0), (1, 1, 1), (0, 1, 0)])
    connec = connect.([(1, 2, 3, 4), (3, 4, 1)], [Tetrahedron, Triangle])
    mesh = SimpleMesh(points, connec)
    topo = topology(mesh)
    @test collect(faces(topo, 2)) == [connect((3, 4, 1), Triangle)]
    @test collect(faces(topo, 3)) == [connect((1, 2, 3, 4), Tetrahedron)]

    # test for https://github.com/JuliaGeometry/Meshes.jl/issues/187
    points = cart.([(0, 0, 0), (1, 0, 0), (1, 1, 1), (0, 1, 0)])
    connec = connect.([(1, 2, 3, 4), (3, 4, 1)], [Tetrahedron, Triangle])
    mesh = SimpleMesh(points[4:-1:1], connec)
    meshvp = SimpleMesh(view(points, 4:-1:1), connec)
    @test mesh == meshvp

    points = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    mesh = SimpleMesh(points, connec)
    bytes = @allocated faces(mesh, 2)
    @test bytes < 100
    cells = faces(mesh, 2)
    bytes = @allocated collect(cells)
    @test bytes < 800

    points = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    mesh = SimpleMesh(points, connec)
    @test centroid(mesh, 1) == centroid(Triangle(cart(0, 0), cart(1, 0), cart(0.5, 0.5)))
    @test centroid(mesh, 2) == centroid(Triangle(cart(1, 0), cart(1, 1), cart(0.5, 0.5)))
    @test centroid(mesh, 3) == centroid(Triangle(cart(1, 1), cart(0, 1), cart(0.5, 0.5)))
    @test centroid(mesh, 4) == centroid(Triangle(cart(0, 1), cart(0, 0), cart(0.5, 0.5)))

    # merge operation with 2D geometries
    mesh₁ = SimpleMesh(cart.([(0, 0), (1, 0), (0, 1)]), connect.([(1, 2, 3)]))
    mesh₂ = SimpleMesh(cart.([(1, 0), (1, 1), (0, 1)]), connect.([(1, 2, 3)]))
    mesh = merge(mesh₁, mesh₂)
    @test vertices(mesh) == [vertices(mesh₁); vertices(mesh₂)]
    @test collect(elements(topology(mesh))) == connect.([(1, 2, 3), (4, 5, 6)])

    # merge operation with 3D geometries
    mesh₁ = SimpleMesh(cart.([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]), connect.([(1, 2, 3, 4)], Tetrahedron))
    mesh₂ = SimpleMesh(cart.([(1, 0, 0), (1, 1, 0), (0, 1, 0), (1, 1, 1)]), connect.([(1, 2, 3, 4)], Tetrahedron))
    mesh = merge(mesh₁, mesh₂)
    @test vertices(mesh) == [vertices(mesh₁); vertices(mesh₂)]
    @test collect(elements(topology(mesh))) == connect.([(1, 2, 3, 4), (5, 6, 7, 8)], Tetrahedron)

    # convert any mesh to SimpleMesh
    grid = cartgrid(10, 10)
    mesh = convert(SimpleMesh, grid)
    @test mesh isa SimpleMesh
    @test topology(mesh) == GridTopology(10, 10)
    @test nvertices(mesh) == 121
    @test nelements(mesh) == 100
    @test eltype(mesh) <: Quadrangle
    # grid interface
    @test size(mesh) == (10, 10)
    @test minimum(mesh) == cart(0, 0)
    @test maximum(mesh) == cart(10, 10)
    @test extrema(mesh) == (cart(0, 0), cart(10, 10))
    @test vertex(mesh, 1) == vertex(mesh, ntuple(i -> 1, embeddim(mesh)))
    @test vertex(mesh, nvertices(mesh)) == vertex(mesh, size(mesh) .+ 1)
    @test mesh[1, 1] == mesh[1]
    @test mesh[10, 10] == mesh[100]
    sub = mesh[2:4, 3:7]
    @test size(sub) == (3, 5)
    @test minimum(sub) == cart(1, 2)
    @test maximum(sub) == cart(4, 7)
    sub = mesh[2, 3:7]
    @test size(sub) == (1, 5)
    @test minimum(sub) == cart(1, 2)
    @test maximum(sub) == cart(2, 7)
    sub = mesh[:, 3:7]
    @test size(sub) == (10, 5)
    @test minimum(sub) == cart(0, 2)
    @test maximum(sub) == cart(10, 7)
    @test_throws BoundsError grid[3:11, :]

    # test for https://github.com/JuliaGeometry/Meshes.jl/issues/261
    points = randpoint2(5)
    connec = [connect((1, 2, 3))]
    mesh = SimpleMesh(points, connec)
    @test nvertices(mesh) == length(vertices(mesh)) == 5

    # single vertex access
    points = randpoint2(5)
    connec = [connect((1, 2, 3))]
    mesh = SimpleMesh(points, connec)
    @test vertex(mesh, 1) == points[1]
    @test vertex(mesh, 2) == points[2]
    @test vertex(mesh, 3) == points[3]
    @test vertex(mesh, 4) == points[4]
    @test vertex(mesh, 5) == points[5]

    points = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    mesh = SimpleMesh(points, connec)
    @test sprint(show, mesh) == "4 SimpleMesh"
    if T == Float32
      @test sprint(show, MIME"text/plain"(), mesh) == """
      4 SimpleMesh
        5 vertices
        ├─ Point(x: 0.0f0 m, y: 0.0f0 m)
        ├─ Point(x: 1.0f0 m, y: 0.0f0 m)
        ├─ Point(x: 0.0f0 m, y: 1.0f0 m)
        ├─ Point(x: 1.0f0 m, y: 1.0f0 m)
        └─ Point(x: 0.5f0 m, y: 0.5f0 m)
        4 elements
        ├─ Triangle(1, 2, 5)
        ├─ Triangle(2, 4, 5)
        ├─ Triangle(4, 3, 5)
        └─ Triangle(3, 1, 5)"""
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), mesh) == """
      4 SimpleMesh
        5 vertices
        ├─ Point(x: 0.0 m, y: 0.0 m)
        ├─ Point(x: 1.0 m, y: 0.0 m)
        ├─ Point(x: 0.0 m, y: 1.0 m)
        ├─ Point(x: 1.0 m, y: 1.0 m)
        └─ Point(x: 0.5 m, y: 0.5 m)
        4 elements
        ├─ Triangle(1, 2, 5)
        ├─ Triangle(2, 4, 5)
        ├─ Triangle(4, 3, 5)
        └─ Triangle(3, 1, 5)"""
    end
  end

  @testset "TransformedMesh" begin
    grid = cartgrid(10, 10)
    rgrid = convert(RectilinearGrid, grid)
    sgrid = convert(StructuredGrid, grid)
    mesh = convert(SimpleMesh, grid)
    trans = Identity()
    tmesh = TransformedMesh(mesh, trans)
    @test Meshes.crs(tmesh) <: Cartesian{NoDatum}
    @test Meshes.lentype(tmesh) == ℳ
    @test parent(tmesh) === mesh
    @test Meshes.transform(tmesh) === trans
    @test TransformedMesh(grid, trans) == grid
    @test TransformedMesh(rgrid, trans) == rgrid
    @test TransformedMesh(sgrid, trans) == sgrid
    @test TransformedMesh(mesh, trans) == mesh
    trans = Translate(T(10), T(10)) → Translate(T(-10), T(-10))
    @test TransformedMesh(grid, trans) == grid
    @test TransformedMesh(rgrid, trans) == rgrid
    @test TransformedMesh(sgrid, trans) == sgrid
    @test TransformedMesh(mesh, trans) == mesh
    trans1 = Translate(T(10), T(10))
    trans2 = Translate(T(-10), T(-10))
    @test TransformedMesh(TransformedMesh(grid, trans1), trans2) == TransformedMesh(grid, trans1 → trans2)
    # grid interface
    trans = Identity()
    tgrid = TransformedMesh(grid, trans)
    @test tgrid isa TransformedGrid
    @test size(tgrid) == (10, 10)
    @test minimum(tgrid) == cart(0, 0)
    @test maximum(tgrid) == cart(10, 10)
    @test extrema(tgrid) == (cart(0, 0), cart(10, 10))
    @test vertex(tgrid, 1) == vertex(tgrid, ntuple(i -> 1, embeddim(tgrid)))
    @test vertex(tgrid, nvertices(tgrid)) == vertex(tgrid, size(tgrid) .+ 1)
    @test tgrid[1, 1] == tgrid[1]
    @test tgrid[10, 10] == tgrid[100]
    sub = tgrid[2:4, 3:7]
    @test size(sub) == (3, 5)
    @test minimum(sub) == cart(1, 2)
    @test maximum(sub) == cart(4, 7)
    sub = tgrid[2, 3:7]
    @test size(sub) == (1, 5)
    @test minimum(sub) == cart(1, 2)
    @test maximum(sub) == cart(2, 7)
    sub = tgrid[:, 3:7]
    @test size(sub) == (10, 5)
    @test minimum(sub) == cart(0, 2)
    @test maximum(sub) == cart(10, 7)
    @test sprint(show, tgrid) == "10×10 TransformedGrid"
    if T == Float32
      @test sprint(show, MIME"text/plain"(), tgrid) == """
      10×10 TransformedGrid
        121 vertices
        ├─ Point(x: 0.0f0 m, y: 0.0f0 m)
        ├─ Point(x: 1.0f0 m, y: 0.0f0 m)
        ├─ Point(x: 2.0f0 m, y: 0.0f0 m)
        ├─ Point(x: 3.0f0 m, y: 0.0f0 m)
        ├─ Point(x: 4.0f0 m, y: 0.0f0 m)
        ⋮
        ├─ Point(x: 6.0f0 m, y: 10.0f0 m)
        ├─ Point(x: 7.0f0 m, y: 10.0f0 m)
        ├─ Point(x: 8.0f0 m, y: 10.0f0 m)
        ├─ Point(x: 9.0f0 m, y: 10.0f0 m)
        └─ Point(x: 10.0f0 m, y: 10.0f0 m)
        100 elements
        ├─ Quadrangle(1, 2, 13, 12)
        ├─ Quadrangle(2, 3, 14, 13)
        ├─ Quadrangle(3, 4, 15, 14)
        ├─ Quadrangle(4, 5, 16, 15)
        ├─ Quadrangle(5, 6, 17, 16)
        ⋮
        ├─ Quadrangle(105, 106, 117, 116)
        ├─ Quadrangle(106, 107, 118, 117)
        ├─ Quadrangle(107, 108, 119, 118)
        ├─ Quadrangle(108, 109, 120, 119)
        └─ Quadrangle(109, 110, 121, 120)"""
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), tgrid) == """
      10×10 TransformedGrid
        121 vertices
        ├─ Point(x: 0.0 m, y: 0.0 m)
        ├─ Point(x: 1.0 m, y: 0.0 m)
        ├─ Point(x: 2.0 m, y: 0.0 m)
        ├─ Point(x: 3.0 m, y: 0.0 m)
        ├─ Point(x: 4.0 m, y: 0.0 m)
        ⋮
        ├─ Point(x: 6.0 m, y: 10.0 m)
        ├─ Point(x: 7.0 m, y: 10.0 m)
        ├─ Point(x: 8.0 m, y: 10.0 m)
        ├─ Point(x: 9.0 m, y: 10.0 m)
        └─ Point(x: 10.0 m, y: 10.0 m)
        100 elements
        ├─ Quadrangle(1, 2, 13, 12)
        ├─ Quadrangle(2, 3, 14, 13)
        ├─ Quadrangle(3, 4, 15, 14)
        ├─ Quadrangle(4, 5, 16, 15)
        ├─ Quadrangle(5, 6, 17, 16)
        ⋮
        ├─ Quadrangle(105, 106, 117, 116)
        ├─ Quadrangle(106, 107, 118, 117)
        ├─ Quadrangle(107, 108, 119, 118)
        ├─ Quadrangle(108, 109, 120, 119)
        └─ Quadrangle(109, 110, 121, 120)"""
    end
    @test_throws BoundsError grid[3:11, :]
  end
end
