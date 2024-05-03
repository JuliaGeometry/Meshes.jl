@testset "Connectivities" begin
  # basic tests
  c = connect((1, 2, 3), Triangle)
  @test pltype(c) == Triangle
  @test paramdim(c) == 2
  @test issimplex(c)
  @test indices(c) == (1, 2, 3)
  @test materialize(c, point.([(0, 0), (1, 0), (0, 1)])) == Triangle(point(0, 0), point(1, 0), point(0, 1))

  # tuple from other collections
  c = connect(Tuple([1, 2, 3]), Triangle)
  @test pltype(c) == Triangle
  @test paramdim(c) == 2
  @test issimplex(c)
  @test indices(c) == (1, 2, 3)
  @test materialize(c, point.([(0, 0), (1, 0), (0, 1)])) == Triangle(point(0, 0), point(1, 0), point(0, 1))

  # incorrect number of vertices for polytope
  @test_throws AssertionError connect((1, 2, 3, 4), Triangle)
  @test_throws AssertionError connect((1, 2, 3), Quadrangle)

  # Ngon requires 3 or more vertices
  @test_throws AssertionError connect((1, 2), Ngon)

  # heterogeneous collections
  c = connect.([(1, 2, 6, 5), (2, 4, 6), (4, 3, 5, 6), (1, 5, 3)], Ngon)
  @test c[1] isa Connectivity{Quadrangle}
  @test c[2] isa Connectivity{Triangle}
  @test c[3] isa Connectivity{Quadrangle}
  @test c[4] isa Connectivity{Triangle}
  @test pltype.(c) == [Quadrangle, Triangle, Quadrangle, Triangle]
  @test issimplex.(c) == [false, true, false, true]

  # ommitting polytope type means polygon or segment
  @test connect((1, 2)) isa Connectivity{Segment}
  @test connect((1, 2, 3)) isa Connectivity{Triangle}
  @test connect((1, 2, 3, 4)) isa Connectivity{Quadrangle}
  @test connect((1, 2, 3, 4, 5)) isa Connectivity{Pentagon}
  @test connect((1, 2, 3, 4, 5, 6)) isa Connectivity{Hexagon}

  # polyhedron connectivities
  c = connect((1, 2, 3, 4), Tetrahedron)
  @test pltype(c) == Tetrahedron
  @test paramdim(c) == 3
  @test issimplex(c)
  @test indices(c) == (1, 2, 3, 4)
  points = point.([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
  @test materialize(c, points) == Tetrahedron(points...)
end
