@testset "Viewing" begin
  @testset "Domain" begin
    d = CartesianGrid{T}(10,10)
    b = Box(P2(1,1), P2(5,5))
    v = view(d, b)
    @test v == CartesianGrid(P2(1,1), P2(5,5), dims=(4,4))

    p = PointSet(collect(vertices(d)))
    v = view(p, b)
    @test coordinates(v, 1) == T[1,1]
    @test coordinates(v, nelements(v)) == T[5,5]
  end

  @testset "Data" begin
    # TODO
  end
end
