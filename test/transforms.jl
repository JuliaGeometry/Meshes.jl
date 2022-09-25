@testset "Transforms" begin
  @testset "TaubinSmoothing" begin
    mesh  = readply(T, joinpath(datadir,"beethoven.ply"))
    trans = TaubinSmoothing(30)
    smesh = trans(mesh)
    # smoothing doesn't change the topology
    @test nvertices(smesh) == nvertices(mesh)
    @test nelements(smesh) == nelements(mesh)
    @test topology(smesh) == topology(mesh)
  end
end
