@testset "Merging" begin
  s = Sphere(P3(0,0,0), T(1))
  c = CylinderSurface(T(1))
  m = merge(s, c)
  @test m isa Multi
  @test eltype(collect(m)) <: Primitive

  s  = Sphere(P3(0,0,0), T(1))
  b  = Box(P3(0,0,0), P3(1,1,1))
  ms = Multi([s])
  mb = Multi([b])
  @test merge(ms, b) == merge(ms, mb) == merge(s, mb)
  m = merge(ms, mb)
  @test m isa Multi
  @test eltype(collect(m)) <: Primitive
end