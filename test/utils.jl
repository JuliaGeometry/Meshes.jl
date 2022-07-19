@testset "Utilities" begin
  a, b, c = P2(0,0), P2(1,0), P2(0,1)
  @test signarea(a, b, c) == T(0.5)
  a, b, c = P2(0,0), P2(0,1), P2(1,0)
  @test signarea(a, b, c) == T(-0.5)

  p1, p2, p3 = P2(0,0), P2(1,1), P2(0.25,0.5)
  s = Segment(P2(0.5,0.0), P2(0.0,1.0))
  @test sideof(p1, s) == :LEFT
  @test sideof(p2, s) == :RIGHT
  @test sideof(p3, s) == :ON

  p1, p2, p3 = P2(0.5,0.5), P2(1.5,0.5), P2(1,1)
  c = Chain(P2[(0,0),(1,0),(1,1),(0,1),(0,0)])
  @test sideof(p1, c) == :INSIDE
  @test sideof(p2, c) == :OUTSIDE
  @test sideof(p3, c) == :INSIDE

  # drop units from unitful value and type
  @test Meshes.dropunits(1.0u"mm") == Float64
  @test Meshes.dropunits(typeof(1.0u"mm")) == Float64
  
  # return the same type in case of no units
  @test Meshes.dropunits(1.0) == Float64
  @test Meshes.dropunits(Float64) == Float64

  normals = [
    V3( 1, 0, 0),
    V3( 0, 1, 0),
    V3( 0, 0, 1),
    V3(-1, 0, 0),
    V3( 0,-1, 0),
    V3( 0, 0,-1),
    V3(rand(3).-0.5)
  ]
  for n in normals    
    u, v = householderbasis(n)
    @test u × v ≈ n ./ norm(n)
  end
end
