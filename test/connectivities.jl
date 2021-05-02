@testset begin "Connectivities"
  # basic tests
  c = connect((1,2,3), Triangle)
  @test paramdim(c) == 2
  @test indices(c) == (1,2,3)
  @test materialize(c, P2[(0,0),(1,0),(0,1)]) == Triangle(P2(0,0), P2(1,0), P2(0,1))

  # tuple from other collections
  c = connect(Tuple([1,2,3]), Triangle)
  @test paramdim(c) == 2
  @test indices(c) == (1,2,3)
  @test materialize(c, P2[(0,0),(1,0),(0,1)]) == Triangle(P2(0,0), P2(1,0), P2(0,1))

  # incorrect number of vertices for polytope
  @test_throws AssertionError connect((1,2,3,4), Triangle)
  @test_throws AssertionError connect((1,2,3), Quadrangle)

  # heterogeneous collections
  c = connect.([(1,2,6,5),(2,4,6),(4,3,5,6),(1,5,3)], Ngon)
  @test c[1] isa Connectivity{Quadrangle}
  @test c[2] isa Connectivity{Triangle}
  @test c[3] isa Connectivity{Quadrangle}
  @test c[4] isa Connectivity{Triangle}
end
