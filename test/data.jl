@testset "Data" begin
  dummydata(domain, table) = DummyData(domain, Dict(paramdim(domain) => table))
  dummymeta(domain, table) = meshdata(domain, Dict(paramdim(domain) => table))

  for (dummy, DummyType) in [(dummydata, DummyData), (dummymeta, MeshData)]
    # fallback constructor with spatial table
    dom = CartesianGrid{T}(2, 2)
    tab = dummydata(dom, (a=[1, 2, 3, 4], b=[5, 6, 7, 8]))
    dat = DummyType(tab)
    @test domain(dat) == domain(tab)
    @test values(dat) == values(tab)

    # equality of data sets
    data₁ = dummy(CartesianGrid{T}(2, 2), (a=[1, 2, 3, 4], b=[5, 6, 7, 8]))
    data₂ = dummy(CartesianGrid{T}(2, 2), (a=[1, 2, 3, 4], b=[5, 6, 7, 8]))
    data₃ = dummy(PointSet(rand(P2, 4)), (a=[1, 2, 3, 4], b=[5, 6, 7, 8]))
    @test data₁ == data₂
    @test data₁ != data₃
    @test data₂ != data₃

    # equality with missing data
    data₁ = dummy(PointSet(T[1 2 3; 4 5 6]), (a=[1, missing, 3], b=[3, 2, 1]))
    data₂ = dummy(PointSet(T[1 2 3; 4 5 6]), (a=[1, missing, 3], b=[3, 2, 1]))
    @test data₁ == data₂

    # Tables interface
    dom = CartesianGrid{T}(2, 2)
    dat = dummy(dom, (a=[1, 2, 3, 4], b=[5, 6, 7, 8]))
    @test Tables.istable(dat)
    sch = Tables.schema(dat)
    @test sch.names == (:a, :b, :geometry)
    @test sch.types == (Int, Int, Quadrangle{2,T})
    @test Tables.rowaccess(dat)
    rows = Tables.rows(dat)
    @test Tables.schema(rows) == sch
    @test collect(rows) == [
      (a=1, b=5, geometry=dom[1]),
      (a=2, b=6, geometry=dom[2]),
      (a=3, b=7, geometry=dom[3]),
      (a=4, b=8, geometry=dom[4])
    ]
    @test collect(Tables.columns(dat)) == [[1, 2, 3, 4], [5, 6, 7, 8], [dom[1], dom[2], dom[3], dom[4]]]
    @test Tables.materializer(dat) <: DummyType

    # dataframe interface
    grid = CartesianGrid{T}(2, 2)
    data = dummy(grid, (a=[1, 2, 3, 4], b=[5, missing, 7, 8]))
    @test propertynames(data) == [:a, :b, :geometry]
    @test isequal(data.a, [1, 2, 3, 4])
    @test isequal(data.b, [5, missing, 7, 8])
    @test data.geometry == grid
    @test_throws ErrorException data.c
    for (a, b, geometry) in [(:a, :b, :geometry), ("a", "b", "geometry")]
      @test data[1:2, [a, b]] == dummy(view(grid, 1:2), (a=[1, 2], b=[5, missing]))
      @test data[1:2, [a, b, geometry]] == dummy(view(grid, 1:2), (a=[1, 2], b=[5, missing]))
      @test isequal(data[1:2, a], [1, 2])
      @test isequal(data[1:2, b], [5, missing])
      @test isequal(data[1:2, geometry], view(grid, 1:2))
      @test data[1:2, :] == dummy(view(grid, 1:2), (a=[1, 2], b=[5, missing]))
      @test isequal(data[1, [a, b]], (a=1, b=5, geometry=grid[1]))
      @test isequal(data[1, [a, b, geometry]], (a=1, b=5, geometry=grid[1]))
      @test isequal(data[1, a], 1)
      @test isequal(data[1, b], 5)
      @test isequal(data[1, geometry], grid[1])
      @test isequal(data[1, :], (a=1, b=5, geometry=grid[1]))
      @test data[:, [a, b]] == data
      @test data[:, [a, b, geometry]] == data
      @test isequal(data[:, a], [1, 2, 3, 4])
      @test isequal(data[:, b], [5, missing, 7, 8])
      @test isequal(data[:, geometry], grid)
    end
    # regex
    @test data[3, r"a"] == (a=3, geometry=grid[3])
    @test data[3:4, r"b"] == dummy(view(grid, 3:4), (b=[7, 8],))
    @test data[:, r"[ab]"] == data

    # hcat
    dom = PointSet(rand(P2, 10))
    data₁ = dummy(dom, (a=rand(10), b=rand(10)))
    data₂ = dummy(dom, (c=rand(10), d=rand(10)))
    data₃ = dummy(dom, (e=rand(10), f=rand(10)))
    @test hcat(data₁) == data₁
    hdata = hcat(data₁, data₂)
    @test hdata.a == data₁.a
    @test hdata.b == data₁.b
    @test hdata.c == data₂.c
    @test hdata.d == data₂.d
    @test hdata.geometry == dom
    hdata = hcat(data₁, data₂, data₃)
    @test hdata.a == data₁.a
    @test hdata.b == data₁.b
    @test hdata.c == data₂.c
    @test hdata.d == data₂.d
    @test hdata.e == data₃.e
    @test hdata.f == data₃.f
    @test hdata.geometry == dom
    # throws
    data₁ = dummy(dom, (a=rand(10), b=rand(10)))
    data₂ = dummy(dom, (a=rand(10), c=rand(10)))
    @test_throws ArgumentError hcat(data₁, data₂)
    data₁ = dummy(PointSet(rand(P2, 10)), (a=rand(10), b=rand(10)))
    data₂ = dummy(PointSet(rand(P2, 10)), (a=rand(10), b=rand(10)))
    @test_throws ArgumentError hcat(data₁, data₂)

    # vcat
    data₁ = dummy(PointSet(rand(P2, 10)), (a=rand(10), b=rand(10)))
    data₂ = dummy(PointSet(rand(P2, 10)), (a=rand(10), b=rand(10)))
    data₃ = dummy(PointSet(rand(P2, 10)), (a=rand(10), b=rand(10)))
    @test vcat(data₁) == data₁
    vdata = vcat(data₁, data₂)
    @test vdata.a == [data₁.a; data₂.a]
    @test vdata.b == [data₁.b; data₂.b]
    @test vdata.geometry == PointSet([collect(data₁.geometry); collect(data₂.geometry)])
    vdata = vcat(data₁, data₂, data₃)
    @test vdata.a == [data₁.a; data₂.a; data₃.a]
    @test vdata.b == [data₁.b; data₂.b; data₃.b]
    @test vdata.geometry == PointSet([collect(data₁.geometry); collect(data₂.geometry); collect(data₃.geometry)])
    # throws
    data₁ = dummy(PointSet(rand(P2, 10)), (a=rand(10), b=rand(10)))
    data₂ = dummy(PointSet(rand(P2, 10)), (a=rand(10), c=rand(10)))
    @test_throws ArgumentError vcat(data₁, data₂)

    # variables interface
    data = dummy(PointSet(rand(P2, 4)), (a=[1, 2, 3, 4], b=[5, 6, 7, 8]))
    @test asarray(data, :a) == asarray(data, "a") == [1, 2, 3, 4]
    @test asarray(data, :b) == asarray(data, "b") == [5, 6, 7, 8]
    data = dummy(CartesianGrid{T}(2, 2), (a=[1, 2, 3, 4], b=[5, 6, 7, 8]))
    @test asarray(data, :a) == asarray(data, "a") == [1 3; 2 4]
    @test asarray(data, :b) == asarray(data, "b") == [5 7; 6 8]

    data = dummy(CartesianGrid{T}(2, 2), (a=[1, 2, 3, 4], b=[5, 6, 7, 8]))
    @test sprint(show, data) == "4 $DummyType"
    if T == Float32
      @test sprint(show, MIME"text/plain"(), data) ==
            "2×2 CartesianGrid{2,Float32}\n  variables (rank 2)\n    └─a (Int64)\n    └─b (Int64)"
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), data) ==
            "2×2 CartesianGrid{2,Float64}\n  variables (rank 2)\n    └─a (Int64)\n    └─b (Int64)"
    end
  end
end
