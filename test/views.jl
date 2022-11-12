@testset "Views" begin
  @testset "Domain" begin
    pset = PointSet(rand(P3, 100))
    inds = rand(1:100, 3)
    v = view(pset, inds)
    @test nelements(v) == 3
    for i in 1:3
      p = pset[inds[i]]
      @test v[i] == p
      @test centroid(v, i) == p
    end

    grid = CartesianGrid{T}(10, 10)
    inds = rand(1:100, 3)
    v = view(grid, inds)
    @test nelements(v) == 3
    for i in 1:3
      e = grid[inds[i]]
      @test v[i] == e
      @test centroid(v, i) == centroid(e)
    end

    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    connec = connect.([(1,2,5),(2,4,5),(4,3,5),(3,1,5)], Triangle)
    mesh = SimpleMesh(points, connec)
    inds = rand(1:4, 3)
    v = view(mesh, inds)
    @test nelements(v) == 3
    for i in 1:3
      e = mesh[inds[i]]
      @test v[i] == e
      @test centroid(v, i) == centroid(e)
    end

    # view of view stores the correct domain
    g = CartesianGrid{T}(10, 10)
    v = view(view(g, 11:20), 1:3)
    @test v isa DomainView{2,T,<:CartesianGrid}
    @test v[1] == g[11]
    @test v[2] == g[12]
    @test v[3] == g[13]

    # centroid of view of PointSet
    points = P2[(0,0), (1,0), (0,1), (1,1), (0.5,0.5)]
    pview  = view(PointSet(points), 1:4)
    @test centroid(pview) == P2(0.5, 0.5)
  end

  @testset "Data" begin
    dummydata(domain, table) = DummyData(domain, Dict(paramdim(domain) => table))
    dummymeta(domain, table) = meshdata(domain, Dict(paramdim(domain) => table))

    for dummy in [dummydata, dummymeta]
      dom = CartesianGrid{T}(2,2)
      dat = dummy(dom, (a=[1,2,3,4], b=[5,6,7,8]))
      v = view(dat, 2:4)
      @test domain(v) == view(dom, 2:4)
      @test Tables.columntable(values(v)) == (a=[2,3,4], b=[6,7,8])
      @test centroid(v, 1) == P2(1.5,0.5)
      @test centroid(v, 2) == P2(0.5,1.5)
      @test centroid(v, 3) == P2(1.5,1.5)
      @test v.a == v."a" == [2,3,4]
      @test v.b == v."b" == [6,7,8]
    end
  end
end
