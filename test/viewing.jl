@testset "Viewing" begin
  @testset "Domain" begin
    g = CartesianGrid{T}(10,10)
    b = Box(P2(1,1), P2(5,5))
    v = view(g, b)
    @test v == CartesianGrid(P2(1,1), P2(5,5), dims=(4,4))

    p = PointSet(collect(vertices(g)))
    v = view(p, b)
    @test centroid(v, 1) == P2(1,1)
    @test centroid(v, nelements(v)) == P2(5,5)
  end

  @testset "Data" begin
    # dummy type implementing the Data trait
    struct DummyData{𝒟,𝒯} <: Data
      domain::𝒟
      table::𝒯
    end
    Meshes.domain(data::DummyData) = getfield(data, :domain)
    Meshes.values(data::DummyData) = getfield(data, :table)
    Meshes.constructor(::Type{D}) where {D<:DummyData} = DummyData

    g = CartesianGrid{T}(10,10)
    t = (a=1:100, b=1:100)
    d = DummyData(g, t)
    b = Box(P2(1,1), P2(5,5))
    v = view(d, b)
    @test domain(v) == CartesianGrid(P2(1,1), P2(5,5), dims=(4,4))
    @test Tables.columntable(values(v)) == (a=[12,13,14,15,22,23,24,25,32,33,34,35,42,43,44,45],
                                            b=[12,13,14,15,22,23,24,25,32,33,34,35,42,43,44,45])

    p = PointSet(collect(vertices(g)))
    d = DummyData(p, t)
    v = view(d, b)
    dd = domain(v)
    @test centroid(dd, 1) == P2(1,1)
    @test centroid(dd, nelements(dd)) == P2(5,5)
    tt = Tables.columntable(values(v))
    @test tt == (a=[13,14,15,16,17,24,25,26,27,28,35,36,37,38,39,46,47,48,49,50,57,58,59,60,61],
                 b=[13,14,15,16,17,24,25,26,27,28,35,36,37,38,39,46,47,48,49,50,57,58,59,60,61])
    
    g = CartesianGrid{T}(250,250)
    t = (a=rand(250*250), b=rand(250*250))
    d = DummyData(g, t)
    s1 = slice(d, T(50.5):T(100.2), T(41.7):T(81.3))
    pts1 = [centroid(s1, i) for i in 1:nelements(s1)]
    X1 = reduce(hcat, coordinates.(pts1))
    @test all(T[50.5,41.7] .≤ minimum(X1, dims=2))
    @test all(maximum(X1, dims=2) .≤ T[100.2,81.3])

    p = sample(d, 100)
    s2 = slice(p, T(50.5):T(150.7), T(175.2):T(250.3))
    pts2 = [centroid(s2, i) for i in 1:nelements(s2)]
    X2 = reduce(hcat, coordinates.(pts2))
    @test all(T[50.5,175.2] .≤ minimum(X2, dims=2))
    @test all(maximum(X2, dims=2) .≤ T[150.7,250.3])
  end
end
