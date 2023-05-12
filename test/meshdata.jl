@testset "MeshData" begin
  etable = (a=rand(100), b=rand(100))
  d = meshdata(CartesianGrid{T}(10, 10), Dict(2 => etable))
  @test domain(d) == CartesianGrid{T}(10, 10)
  @test values(d, 0) === nothing
  @test values(d, 1) === nothing
  @test values(d, 2) == etable

  vertices = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
  elements = connect.([(1, 2, 3), (3, 4, 1)])
  d = meshdata(vertices, elements,
    Dict(
      0 => (temperature=[1.0, 2.0, 3.0, 4.0], pressure=[4.0, 3.0, 2.0, 1.0]),
      2 => (quality=["A", "B"], state=[true, false])
    )
  )
  @test values(d, 0) == (temperature=[1.0, 2.0, 3.0, 4.0], pressure=[4.0, 3.0, 2.0, 1.0])
  @test values(d, 1) === nothing
  @test values(d, 2) == (quality=["A", "B"], state=[true, false])

  etable = (a=rand(100), b=rand(100))
  d = meshdata(CartesianGrid{T}(10, 10), etable=etable)
  @test domain(d) == CartesianGrid{T}(10, 10)
  @test values(d, 0) === nothing
  @test values(d, 1) === nothing
  @test values(d, 2) == etable

  vertices = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
  elements = connect.([(1, 2, 3), (3, 4, 1)])
  d = meshdata(vertices, elements,
    vtable=(temperature=[1.0, 2.0, 3.0, 4.0], pressure=[4.0, 3.0, 2.0, 1.0]),
    etable=(quality=["A", "B"], state=[true, false])
  )
  @test domain(d) == SimpleMesh(vertices, elements)
  @test values(d, 0) == (temperature=[1.0, 2.0, 3.0, 4.0], pressure=[4.0, 3.0, 2.0, 1.0])
  @test values(d, 1) === nothing
  @test values(d, 2) == (quality=["A", "B"], state=[true, false])

  etable = [(a=1.0, b=2.0), (a=3.0, b=4.0), (a=5.0, b=6.0), (a=7.0, b=8.0)]
  d = meshdata(CartesianGrid{T}(2, 2), etable=etable)
  @test d.a == [1.0, 3.0, 5.0, 7.0]
  @test d.b == [2.0, 4.0, 6.0, 8.0]
end
