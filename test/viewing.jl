@testset "Viewing" begin
  @testset "Domain" begin
    g = CartesianGrid{T}(10,10)
    b = Box(P2(1,1), P2(5,5))
    v = view(g, b)
    @test v == CartesianGrid(P2(1,1), P2(5,5), dims=(4,4))

    p = PointSet(collect(vertices(g)))
    v = view(p, b)
    @test coordinates(v, 1) == T[1,1]
    @test coordinates(v, nelements(v)) == T[5,5]
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

    g = CartesianGrid{T}(10,10)
    t = (a=1:100, b=1:100)
    d = DummyData(g, t)
    b = Box(P2(1,1), P2(5,5))
    v = view(d, b)
    @test domain(v) == CartesianGrid(P2(1,1), P2(5,5), dims=(4,4))
    @test values(v) == (a=[12,13,14,15,22,23,24,25,32,33,34,35,42,43,44,45], b=[12,13,14,15,22,23,24,25,32,33,34,35,42,43,44,45])

    p = PointSet(collect(vertices(g)))
    d = DummyData(p, t)
    v = view(d, b)
    dd = domain(v)
    @test coordinates(dd, 1) == T[1,1]
    @test coordinates(dd, nelements(dd)) == T[5,5]
    tt = values(v)
    @test tt == (a=[13,14,15,16,17,24,25,26,27,28,35,36,37,38,39,46,47,48,49,50,57,58,59,60,61],
                 b=[13,14,15,16,17,24,25,26,27,28,35,36,37,38,39,46,47,48,49,50,57,58,59,60,61])
  end
end
