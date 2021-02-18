@testset "Views" begin
  @testset "Domain" begin
    pset = PointSet(rand(P3, 100))
    inds = rand(1:100, 3)
    v = view(pset, inds)
    @test nelements(v) == 3
    for i in 1:3
      p = pset[inds[i]]
      @test v[i] == p
      @test coordinates(v, i) == coordinates(p)
    end

    grid = CartesianGrid{T}(10, 10)
    inds = rand(1:100, 3)
    v = view(grid, inds)
    @test nelements(v) == 3
    for i in 1:3
      e = grid[inds[i]]
      @test v[i] == e
      @test coordinates(v, i) == coordinates(centroid(e))
    end

    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = UnstructuredMesh(points, connec)
    inds = rand(1:4, 3)
    v = view(mesh, inds)
    @test nelements(v) == 3
    for i in 1:3
      e = mesh[inds[i]]
      @test v[i] == e
      @test coordinates(v, i) == coordinates(centroid(e))
    end

    if visualtests
      d = CartesianGrid{T}(10, 10)
      v = view(d, 1:50)
      @test_ref_plot "data/domain-view-$T.png" plot(v)
    end
  end

  @testset "Data" begin
    # dummy type implementing the Data trait
    struct DummyData{𝒟,𝒯} <: Data
      domain::𝒟
      table::𝒯
    end
    Meshes.domain(data::DummyData) = data.domain
    Meshes.values(data::DummyData) = data.table

    dom = CartesianGrid{T}(2,2)
    dat = DummyData(dom, (a=[1,2,3,4], b=[5,6,7,8]))

    v = view(dat, 2:4)
    @test domain(v) == view(dom, 2:4)
    @test values(v) == (a=[2,3,4], b=[6,7,8])
    @test coordinates(domain(v), 1:nelements(domain(v))) == T[1.5 0.5 1.5; 0.5 1.5 1.5]
    @test v[:a] == [2,3,4]
    @test v[:b] == [6,7,8]

    v = view(dat, [:a])
    @test domain(v) == view(dom, 1:4)
    @test values(v) == (a=[1,2,3,4],)
    @test coordinates(domain(v), 1:nelements(domain(v))) == T[0.5 1.5 0.5 1.5; 0.5 0.5 1.5 1.5]
    @test v[:a] == [1,2,3,4]
    @test_throws ErrorException v[:b]

    v = view(dat, 1:3, [:a])
    @test domain(v) == view(dom, 1:3)
    @test values(v) == (a=[1,2,3],)
    @test coordinates(domain(v), 1:nelements(domain(v))) == T[0.5 1.5 0.5; 0.5 0.5 1.5]
    @test v[:a] == [1,2,3]
    @test_throws ErrorException v[:b]
  end
end
