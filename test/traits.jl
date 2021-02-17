@testset "Traits" begin
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
  end
end
