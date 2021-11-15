@testset "Smoothing" begin
  @testset "TaubinSmoothing" begin
    mesh  = readply(T, joinpath(datadir,"beethoven.ply"))
    smesh = smooth(mesh, TaubinSmoothing(30))
    # smoothing doesn't change the topology
    @test nvertices(smesh) == nvertices(mesh)
    @test nelements(smesh) == nelements(mesh)
    @test topology(smesh) == topology(mesh)
  end
end
