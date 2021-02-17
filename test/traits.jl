@testset "Traits" begin
  @testset "Domain" begin
    # TODO
  end

  @testset "Data" begin
    # dummy type implementing the Data trait
    struct DummyData{𝒟,𝒯} <: Data
      domain::𝒟
      table::𝒯
    end
    Meshes.domain(data::DummyData) = data.domain
    Meshes.values(data::DummyData) = data.table

    # equality of data sets
    data₁ = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,6,7,8]))
    data₂ = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,6,7,8]))
    data₃ = DummyData(PointSet(rand(P2,4)), (a=[1,2,3,4], b=[5,6,7,8]))
    @test data₁ == data₂
    @test data₁ != data₃
    @test data₂ != data₃

    # Tables interface
    data = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,6,7,8]))
    @test Tables.istable(data)
    @test Tables.rowaccess(data)
    s = Tables.schema(data)
    @test s.names == (:a,:b,:geometry)
    @test s.types == (Int, Int, Quadrangle{2,T,Vector{P2}})

    # variables interface
    data = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,6,7,8]))
    @test data[:a] == [1,2,3,4]
    @test data[:b] == [5,6,7,8]
    @test_throws ErrorException data[:c] 
    @test_throws ErrorException data[:geometry]

    # utility functions
    data = DummyData(PointSet(rand(P2,4)), (a=[1,2,3,4], b=[5,6,7,8]))
    @test asarray(data, :a) == [1,2,3,4]
    @test asarray(data, :b) == [5,6,7,8]
    data = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,6,7,8]))
    @test asarray(data, :a) == [1 3; 2 4]
    @test asarray(data, :b) == [5 7; 6 8]

    data = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,6,7,8]))
    @test sprint(show, data) == "4 DummyData{2,$T}"
    if T == Float32
      @test sprint(show, MIME"text/plain"(), data) == "4 DummyData{2,Float32}\n  variables\n    └─a (Int64)\n    └─b (Int64)\n  domain: 2×2 CartesianGrid{2,Float32}"
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), data) == "4 DummyData{2,Float64}\n  variables\n    └─a (Int64)\n    └─b (Int64)\n  domain: 2×2 CartesianGrid{2,Float64}"
    end

    if visualtests
      data = DummyData(CartesianGrid{T}(2,2), (a=[1,2,3,4], b=[5,6,7,8]))
      @test_ref_plot "data/data-$T.png" plot(data)
    end
  end
end
