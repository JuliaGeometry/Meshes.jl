@testset "Traits" begin
  @testset "Domain" begin
    # dummy type implementing the Domain trait
    struct DummyDomain{Dim,T} <: Domain{Dim,T}
      origin::Point{Dim,T}
    end
    function Base.getindex(domain::DummyDomain{Dim,T}, ind::Int) where {Dim,T}
      c = domain.origin + Vec(ntuple(i->T(ind), Dim))
      r = one(T)
      Ball(c, r)
    end
    Meshes.nelements(d::DummyDomain) = 3

    # basic properties
    dom = DummyDomain(P2(0,0))
    @test embeddim(dom) == 2
    @test coordtype(dom) == T

    # indexable/iterable interface
    dom = DummyDomain(P2(0,0))
    @test dom[begin] == Ball(P2(1,1), T(1))
    @test dom[end]   == Ball(P2(3,3), T(1))
    @test eltype(dom) <: Ball{2,T}
    @test length(dom) == 3
    @test collect(dom) == [Ball(P2(i,i), T(1)) for i in 1:3]

    # coordinates of centroids
    dom = DummyDomain(P2(1,1))
    pts = centroid.(Ref(dom), 1:3)
    @test pts == P2[(2,2), (3,3), (4,4)]

    dom = DummyDomain(P2(0,0))
    @test sprint(show, dom) == "3 DummyDomain{2,$T}"
    if T == Float32
      @test sprint(show, MIME"text/plain"(), dom) == "3 DummyDomain{2,Float32}\n  └─Ball{2,Float32}(Point(1.0f0, 1.0f0), 1.0))\n  └─Ball{2,Float32}(Point(2.0f0, 2.0f0), 1.0))\n  └─Ball{2,Float32}(Point(3.0f0, 3.0f0), 1.0))"
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), dom) == "3 DummyDomain{2,Float64}\n  └─Ball{2,Float64}(Point(1.0, 1.0), 1.0))\n  └─Ball{2,Float64}(Point(2.0, 2.0), 1.0))\n  └─Ball{2,Float64}(Point(3.0, 3.0), 1.0))"
    end

    if visualtests
      dom = DummyDomain(P2(0,0))
      @test_reference "data/domain-$T.png" plot(dom)
      @test_reference "data/domain-data-$T.png" plot(dom,1:3)
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
    Meshes.constructor(::Type{D}) where {D<:DummyData} = DummyData

    # fallback constructor with spatial table
    dom = CartesianGrid{T}(2,2)
    tab = DummyData(dom, (a=[1,2,3,4], b=[5,6,7,8]))
    dat = DummyData(tab)
    @test domain(dat) == domain(tab)
    @test values(dat) == values(tab)

    # equality of data sets
    data₁ = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,6,7,8]))
    data₂ = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,6,7,8]))
    data₃ = DummyData(PointSet(rand(P2,4)), (a=[1,2,3,4], b=[5,6,7,8]))
    @test data₁ == data₂
    @test data₁ != data₃
    @test data₂ != data₃

    # Tables interface
    dom = CartesianGrid{T}(2,2)
    dat = DummyData(dom, (a=[1,2,3,4], b=[5,6,7,8]))
    @test Tables.istable(dat)
    @test Tables.rowaccess(dat)
    rows = Tables.rows(dat)
    schema = Tables.schema(rows)
    @test schema.names == (:a,:b,:geometry)
    @test schema.types == (Int, Int, Quadrangle{2,T,Vector{P2}})
    @test collect(rows) == [
      (a=1, b=5, geometry=dom[1]),
      (a=2, b=6, geometry=dom[2]),
      (a=3, b=7, geometry=dom[3]),
      (a=4, b=8, geometry=dom[4])
    ]
    @test collect(Tables.columns(dat)) == [
      [1,2,3,4],
      [5,6,7,8],
      [dom[1],dom[2],dom[3],dom[4]]
    ]
    @test Tables.materializer(dat) <: DummyData

    # Query interface
    dom = CartesianGrid{T}(2,2)
    dat = DummyData(dom, (a=[1,2,3,4], b=[5,6,7,8]))
    new = dat |> @mutate(geometry=centroid(_.geometry)) |> DummyData
    @test domain(new) isa PointSet
    @test values(new) == values(dat)

    # variables interface
    data = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,missing,7,8]))
    @test variables(data) == (Variable(:a,Int), Variable(:b,Int))
    @test name.(variables(data)) == (:a,:b)
    @test mactype.(variables(data)) == (Int,Int)
    @test data[:a] == data["a"] == [1,2,3,4]
    @test isequal(data[:b], [5,missing,7,8])
    @test_throws ErrorException data[:c] 
    @test_throws ErrorException data[:geometry]

    # utility functions
    data = DummyData(PointSet(rand(P2,4)), (a=[1,2,3,4], b=[5,6,7,8]))
    @test asarray(data, :a) == asarray(data, "a") == [1,2,3,4]
    @test asarray(data, :b) == asarray(data, "b") == [5,6,7,8]
    data = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,6,7,8]))
    @test asarray(data, :a) == asarray(data, "a") == [1 3; 2 4]
    @test asarray(data, :b) == asarray(data, "b") == [5 7; 6 8]

    data = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,6,7,8]))
    @test sprint(show, data) == "4 DummyData{2,$T}"
    if T == Float32
      @test sprint(show, MIME"text/plain"(), data) == "4 DummyData{2,Float32}\n  variables\n    └─a (Int64)\n    └─b (Int64)\n  domain: 2×2 CartesianGrid{2,Float32}"
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), data) == "4 DummyData{2,Float64}\n  variables\n    └─a (Int64)\n    └─b (Int64)\n  domain: 2×2 CartesianGrid{2,Float64}"
    end

    if visualtests
      data = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,6,7,8]))
      @test_reference "data/data-$T.png" plot(data)
      data = DummyData(CartesianGrid{T}(2,2), (c=categorical([1,2,3,4]),))
      @test_reference "data/data-categorical-$T.png" plot(data)
    end
  end
end
