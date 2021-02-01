@testset "Discretization" begin
  @testset "FIST" begin
    𝒫 = Chain(P2[(0,0),(1,0),(1,1),(2,1),(2,2),(1,2),(0,0)])
    @test Meshes.ears(𝒫) == [5]

    𝒫 = Chain(P2[(0,0),(1,0),(1,1),(2,1),(1,2),(0,0)])
    @test Meshes.ears(𝒫) == [4]

    𝒫 = Chain(P2[(0,0),(1,0),(1,1),(1,2),(0,0)])
    @test Meshes.ears(𝒫) == [2,4]

    𝒫 = Chain(P2[(0,0),(1,1),(1,2),(0,0)])
    @test Meshes.ears(𝒫) == []

    points = P2[(0,0),(1,0),(1,1),(2,1),(2,2),(1,2),(0,0)]
    connec = connect.([(4,5,6),(3,4,6),(3,6,1),(1,2,3)], Triangle)
    target = UnstructuredMesh(points[1:end-1], connec)
    𝒫 = PolyArea(points)
    mesh = discretize(𝒫, FIST())
    @test mesh == target
  end
end
