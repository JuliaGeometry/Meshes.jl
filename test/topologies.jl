@testitem "GridTopology" setup = [Setup] begin
  t = GridTopology(3)
  @test paramdim(t) == 1
  @test size(t) == (3,)
  @test elementtype(t) == Segment
  @test facettype(t) == Point
  @test elem2cart(t, 1) == (1,)
  @test elem2cart(t, 2) == (2,)
  @test elem2cart(t, 3) == (3,)
  @test cart2corner(t, 1) == 1
  @test cart2corner(t, 2) == 2
  @test cart2corner(t, 3) == 3
  @test elem2corner(t, 1) == 1
  @test elem2corner(t, 2) == 2
  @test elem2corner(t, 3) == 3
  @test corner2elem(t, 1) == 1
  @test corner2elem(t, 2) == 2
  @test corner2elem(t, 3) == 3
  @test nelements(t) == 3
  @test nfacets(t) == 4
  @test nvertices(t) == 4
  @test nfaces(t, 1) == 3
  @test nfaces(t, 0) == 4
  @test element(t, 1) == connect((1, 2))
  @test element(t, 2) == connect((2, 3))
  @test element(t, 3) == connect((3, 4))
  @test faces(t, 1) == elements(t)
  @test faces(t, 0) == vertices(t)
  @test vertices(t) == 1:4
  @test vertex(t, 1) == 1
  @test vertex(t, 4) == 4

  t = GridTopology(3, 4)
  @test paramdim(t) == 2
  @test size(t) == (3, 4)
  @test elementtype(t) == Quadrangle
  @test facettype(t) == Segment
  @test elem2cart(t, 1) == (1, 1)
  @test elem2cart(t, 2) == (2, 1)
  @test elem2cart(t, 3) == (3, 1)
  @test elem2cart(t, 4) == (1, 2)
  @test elem2cart(t, 5) == (2, 2)
  @test elem2cart(t, 6) == (3, 2)
  @test elem2cart(t, 7) == (1, 3)
  @test elem2cart(t, 8) == (2, 3)
  @test elem2cart(t, 9) == (3, 3)
  @test elem2cart(t, 10) == (1, 4)
  @test elem2cart(t, 11) == (2, 4)
  @test elem2cart(t, 12) == (3, 4)
  @test cart2corner(t, 1, 1) == 1
  @test cart2corner(t, 2, 1) == 2
  @test cart2corner(t, 3, 1) == 3
  @test cart2corner(t, 1, 2) == 5
  @test cart2corner(t, 2, 2) == 6
  @test cart2corner(t, 3, 2) == 7
  @test cart2corner(t, 1, 3) == 9
  @test cart2corner(t, 2, 3) == 10
  @test cart2corner(t, 3, 3) == 11
  @test cart2corner(t, 1, 4) == 13
  @test cart2corner(t, 2, 4) == 14
  @test cart2corner(t, 3, 4) == 15
  @test elem2corner(t, 1) == 1
  @test elem2corner(t, 2) == 2
  @test elem2corner(t, 3) == 3
  @test elem2corner(t, 4) == 5
  @test elem2corner(t, 5) == 6
  @test elem2corner(t, 6) == 7
  @test elem2corner(t, 7) == 9
  @test elem2corner(t, 8) == 10
  @test elem2corner(t, 9) == 11
  @test elem2corner(t, 10) == 13
  @test elem2corner(t, 11) == 14
  @test elem2corner(t, 12) == 15
  @test corner2elem(t, 1) == 1
  @test corner2elem(t, 2) == 2
  @test corner2elem(t, 3) == 3
  @test corner2elem(t, 5) == 4
  @test corner2elem(t, 6) == 5
  @test corner2elem(t, 7) == 6
  @test corner2elem(t, 9) == 7
  @test corner2elem(t, 10) == 8
  @test corner2elem(t, 11) == 9
  @test corner2elem(t, 13) == 10
  @test corner2elem(t, 14) == 11
  @test corner2elem(t, 15) == 12
  @test nelements(t) == 12
  @test nfacets(t) == 31
  @test nvertices(t) == 20
  @test nfaces(t, 2) == 12
  @test nfaces(t, 1) == 31
  @test nfaces(t, 0) == 20
  @test element(t, 1) == connect((1, 2, 6, 5))
  @test element(t, 5) == connect((6, 7, 11, 10))
  @test faces(t, 2) == elements(t)
  @test faces(t, 0) == vertices(t)
  @test vertices(t) == 1:20
  @test vertex(t, 1) == 1
  @test vertex(t, 20) == 20
  @test facet.(Ref(t), 1:31) ==
        connect.([
    (1, 5),
    (2, 6),
    (3, 7),
    (4, 8),
    (5, 9),
    (6, 10),
    (7, 11),
    (8, 12),
    (9, 13),
    (10, 14),
    (11, 15),
    (12, 16),
    (13, 17),
    (14, 18),
    (15, 19),
    (16, 20),
    (1, 2),
    (5, 6),
    (9, 10),
    (13, 14),
    (17, 18),
    (2, 3),
    (6, 7),
    (10, 11),
    (14, 15),
    (18, 19),
    (3, 4),
    (7, 8),
    (11, 12),
    (15, 16),
    (19, 20)
  ])

  t = GridTopology(3, 4, 2)
  @test paramdim(t) == 3
  @test size(t) == (3, 4, 2)
  @test elementtype(t) == Hexahedron
  @test facettype(t) == Quadrangle
  @test elem2cart(t, 1) == (1, 1, 1)
  @test elem2cart(t, 2) == (2, 1, 1)
  @test elem2cart(t, 3) == (3, 1, 1)
  @test elem2cart(t, 4) == (1, 2, 1)
  @test elem2cart(t, 5) == (2, 2, 1)
  @test elem2cart(t, 6) == (3, 2, 1)
  @test elem2cart(t, 7) == (1, 3, 1)
  @test elem2cart(t, 8) == (2, 3, 1)
  @test elem2cart(t, 9) == (3, 3, 1)
  @test elem2cart(t, 10) == (1, 4, 1)
  @test elem2cart(t, 11) == (2, 4, 1)
  @test elem2cart(t, 12) == (3, 4, 1)
  @test elem2cart(t, 13) == (1, 1, 2)
  @test elem2cart(t, 14) == (2, 1, 2)
  @test elem2cart(t, 15) == (3, 1, 2)
  @test elem2cart(t, 16) == (1, 2, 2)
  @test elem2cart(t, 17) == (2, 2, 2)
  @test elem2cart(t, 18) == (3, 2, 2)
  @test elem2cart(t, 19) == (1, 3, 2)
  @test elem2cart(t, 20) == (2, 3, 2)
  @test elem2cart(t, 21) == (3, 3, 2)
  @test elem2cart(t, 22) == (1, 4, 2)
  @test elem2cart(t, 23) == (2, 4, 2)
  @test elem2cart(t, 24) == (3, 4, 2)
  @test cart2corner(t, 1, 1, 1) == 1
  @test cart2corner(t, 2, 1, 1) == 2
  @test cart2corner(t, 3, 1, 1) == 3
  @test cart2corner(t, 1, 2, 1) == 5
  @test cart2corner(t, 2, 2, 1) == 6
  @test cart2corner(t, 3, 2, 1) == 7
  @test cart2corner(t, 1, 3, 1) == 9
  @test cart2corner(t, 2, 3, 1) == 10
  @test cart2corner(t, 3, 3, 1) == 11
  @test cart2corner(t, 1, 4, 1) == 13
  @test cart2corner(t, 2, 4, 1) == 14
  @test cart2corner(t, 3, 4, 1) == 15
  @test cart2corner(t, 1, 1, 2) == 21
  @test cart2corner(t, 2, 1, 2) == 22
  @test cart2corner(t, 3, 1, 2) == 23
  @test cart2corner(t, 1, 2, 2) == 25
  @test cart2corner(t, 2, 2, 2) == 26
  @test cart2corner(t, 3, 2, 2) == 27
  @test cart2corner(t, 1, 3, 2) == 29
  @test cart2corner(t, 2, 3, 2) == 30
  @test cart2corner(t, 3, 3, 2) == 31
  @test cart2corner(t, 1, 4, 2) == 33
  @test cart2corner(t, 2, 4, 2) == 34
  @test cart2corner(t, 3, 4, 2) == 35
  @test elem2corner(t, 1) == 1
  @test elem2corner(t, 2) == 2
  @test elem2corner(t, 3) == 3
  @test elem2corner(t, 4) == 5
  @test elem2corner(t, 5) == 6
  @test elem2corner(t, 6) == 7
  @test elem2corner(t, 7) == 9
  @test elem2corner(t, 8) == 10
  @test elem2corner(t, 9) == 11
  @test elem2corner(t, 10) == 13
  @test elem2corner(t, 11) == 14
  @test elem2corner(t, 12) == 15
  @test elem2corner(t, 13) == 21
  @test elem2corner(t, 14) == 22
  @test elem2corner(t, 15) == 23
  @test elem2corner(t, 16) == 25
  @test elem2corner(t, 17) == 26
  @test elem2corner(t, 18) == 27
  @test elem2corner(t, 19) == 29
  @test elem2corner(t, 20) == 30
  @test elem2corner(t, 21) == 31
  @test elem2corner(t, 22) == 33
  @test elem2corner(t, 23) == 34
  @test elem2corner(t, 24) == 35
  @test corner2elem(t, 1) == 1
  @test corner2elem(t, 2) == 2
  @test corner2elem(t, 3) == 3
  @test corner2elem(t, 5) == 4
  @test corner2elem(t, 6) == 5
  @test corner2elem(t, 7) == 6
  @test corner2elem(t, 9) == 7
  @test corner2elem(t, 10) == 8
  @test corner2elem(t, 11) == 9
  @test corner2elem(t, 13) == 10
  @test corner2elem(t, 14) == 11
  @test corner2elem(t, 15) == 12
  @test corner2elem(t, 21) == 13
  @test corner2elem(t, 22) == 14
  @test corner2elem(t, 23) == 15
  @test corner2elem(t, 25) == 16
  @test corner2elem(t, 26) == 17
  @test corner2elem(t, 27) == 18
  @test corner2elem(t, 29) == 19
  @test corner2elem(t, 30) == 20
  @test corner2elem(t, 31) == 21
  @test corner2elem(t, 33) == 22
  @test corner2elem(t, 34) == 23
  @test corner2elem(t, 35) == 24
  @test nelements(t) == 24
  @test nfacets(t) == 3 * 24 + 3 * 4 + 4 * 2 + 3 * 2
  @test nvertices(t) == 60
  @test nfaces(t, 3) == 24
  @test nfaces(t, 2) == 3 * 24 + 3 * 4 + 4 * 2 + 3 * 2
  @test nfaces(t, 0) == 60
  @test element(t, 1) == connect((1, 2, 6, 5, 21, 22, 26, 25), Hexahedron)
  @test element(t, 5) == connect((6, 7, 11, 10, 26, 27, 31, 30), Hexahedron)
  @test faces(t, 3) == elements(t)
  @test faces(t, 0) == vertices(t)
  @test vertices(t) == 1:60
  @test vertex(t, 1) == 1
  @test vertex(t, 60) == 60

  t = GridTopology((3,), (true,))
  @test isperiodic(t) == (true,)
  @test nvertices(t) == 3
  @test nelements(t) == 3
  @test nfacets(t) == 3
  @test element(t, 1) == connect((1, 2))
  @test element(t, 2) == connect((2, 3))
  @test element(t, 3) == connect((3, 1))

  t = GridTopology((2, 3), (true, true))
  @test isperiodic(t) == (true, true)
  @test nvertices(t) == 2 * 3
  @test nelements(t) == 6
  @test nfacets(t) == 12
  @test element(t, 1) == connect((1, 2, 4, 3))
  @test element(t, 2) == connect((2, 1, 3, 4))
  @test element(t, 3) == connect((3, 4, 6, 5))
  @test element(t, 4) == connect((4, 3, 5, 6))
  @test element(t, 5) == connect((5, 6, 2, 1))
  @test element(t, 6) == connect((6, 5, 1, 2))

  t = GridTopology((2, 3), (false, true))
  @test isperiodic(t) == (false, true)
  @test nvertices(t) == 3 * 3
  @test nelements(t) == 6
  @test nfacets(t) == 15
  @test element(t, 1) == connect((1, 2, 5, 4))
  @test element(t, 2) == connect((2, 3, 6, 5))
  @test element(t, 3) == connect((4, 5, 8, 7))
  @test element(t, 4) == connect((5, 6, 9, 8))
  @test element(t, 5) == connect((7, 8, 2, 1))
  @test element(t, 6) == connect((8, 9, 3, 2))

  t = GridTopology((2, 3), (true, false))
  @test isperiodic(t) == (true, false)
  @test nvertices(t) == 2 * 4
  @test nelements(t) == 6
  @test nfacets(t) == 14
  @test element(t, 1) == connect((1, 2, 4, 3))
  @test element(t, 2) == connect((2, 1, 3, 4))
  @test element(t, 3) == connect((3, 4, 6, 5))
  @test element(t, 4) == connect((4, 3, 5, 6))
  @test element(t, 5) == connect((5, 6, 8, 7))
  @test element(t, 6) == connect((6, 5, 7, 8))

  t = GridTopology((2, 3, 4), (true, true, true))
  @test isperiodic(t) == (true, true, true)
  @test nvertices(t) == 2 * 3 * 4
  @test nelements(t) == 2 * 3 * 4
  @test nfacets(t) == 3 * (2 * 3 * 4)
  @test element(t, 1) == connect((1, 2, 4, 3, 7, 8, 10, 9), Hexahedron)
  @test element(t, 2) == connect((2, 1, 3, 4, 8, 7, 9, 10), Hexahedron)
  @test element(t, 24) == connect((24, 23, 19, 20, 6, 5, 1, 2), Hexahedron)

  t = GridTopology((2, 3, 4), (false, true, true))
  @test isperiodic(t) == (false, true, true)
  @test nvertices(t) == 3 * 3 * 4
  @test nelements(t) == 2 * 3 * 4
  @test nfacets(t) == 3 * (2 * 3 * 4) + 3 * 4
  @test element(t, 1) == connect((1, 2, 5, 4, 10, 11, 14, 13), Hexahedron)
  @test element(t, 2) == connect((2, 3, 6, 5, 11, 12, 15, 14), Hexahedron)
  @test element(t, 24) == connect((35, 36, 30, 29, 8, 9, 3, 2), Hexahedron)

  t = GridTopology((2, 3, 4), (true, false, true))
  @test isperiodic(t) == (true, false, true)
  @test nvertices(t) == 2 * 4 * 4
  @test nelements(t) == 2 * 3 * 4
  @test nfacets(t) == 3 * (2 * 3 * 4) + 2 * 4
  @test element(t, 1) == connect((1, 2, 4, 3, 9, 10, 12, 11), Hexahedron)
  @test element(t, 2) == connect((2, 1, 3, 4, 10, 9, 11, 12), Hexahedron)
  @test element(t, 24) == connect((30, 29, 31, 32, 6, 5, 7, 8), Hexahedron)

  t = GridTopology((2, 3, 4), (true, true, false))
  @test isperiodic(t) == (true, true, false)
  @test nvertices(t) == 2 * 3 * 5
  @test nelements(t) == 2 * 3 * 4
  @test nfacets(t) == 3 * (2 * 3 * 4) + 2 * 3
  @test element(t, 1) == connect((1, 2, 4, 3, 7, 8, 10, 9), Hexahedron)
  @test element(t, 2) == connect((2, 1, 3, 4, 8, 7, 9, 10), Hexahedron)
  @test element(t, 24) == connect((24, 23, 19, 20, 30, 29, 25, 26), Hexahedron)

  t = GridTopology((2, 3, 4), (true, false, false))
  @test isperiodic(t) == (true, false, false)
  @test nvertices(t) == 2 * 4 * 5
  @test nelements(t) == 2 * 3 * 4
  @test nfacets(t) == 3 * (2 * 3 * 4) + 2 * 4 + 2 * 3
  @test element(t, 1) == connect((1, 2, 4, 3, 9, 10, 12, 11), Hexahedron)
  @test element(t, 2) == connect((2, 1, 3, 4, 10, 9, 11, 12), Hexahedron)
  @test element(t, 24) == connect((30, 29, 31, 32, 38, 37, 39, 40), Hexahedron)

  t = GridTopology((2, 3, 4), (false, true, false))
  @test isperiodic(t) == (false, true, false)
  @test nvertices(t) == 3 * 3 * 5
  @test nelements(t) == 2 * 3 * 4
  @test nfacets(t) == 3 * (2 * 3 * 4) + 3 * 4 + 2 * 3
  @test element(t, 1) == connect((1, 2, 5, 4, 10, 11, 14, 13), Hexahedron)
  @test element(t, 2) == connect((2, 3, 6, 5, 11, 12, 15, 14), Hexahedron)
  @test element(t, 24) == connect((35, 36, 30, 29, 44, 45, 39, 38), Hexahedron)

  t = GridTopology((2, 3, 4), (false, false, true))
  @test isperiodic(t) == (false, false, true)
  @test nvertices(t) == 3 * 4 * 4
  @test nelements(t) == 2 * 3 * 4
  @test nfacets(t) == 3 * (2 * 3 * 4) + 3 * 4 + 2 * 4
  @test element(t, 1) == connect((1, 2, 5, 4, 13, 14, 17, 16), Hexahedron)
  @test element(t, 2) == connect((2, 3, 6, 5, 14, 15, 18, 17), Hexahedron)
  @test element(t, 24) == connect((44, 45, 48, 47, 8, 9, 12, 11), Hexahedron)

  # indexable api
  t = GridTopology(10, 10)
  @test t[begin] == connect((1, 2, 13, 12), Quadrangle)
  @test t[end] == connect((109, 110, 121, 120), Quadrangle)
  @test t[10] == connect((10, 11, 22, 21), Quadrangle)
  @test length(t) == 100
  @test eltype(t) == Connectivity{Quadrangle,4}
  for e in t
    @test e isa Connectivity{Quadrangle,4}
  end
end

@testitem "HalfEdgeTopology" setup = [Setup] begin
  function test_halfedge(elems, topology)
    @test nelements(topology) == length(elems)
    for e in 1:nelements(topology)
      he = half4elem(topology, e)
      inds = indices(elems[e])
      for _ in inds
        @test he.elem == e
        @test he.head ∈ inds
        @test he.next.elem == e
        @test he.prev.elem == e
        @test he.next.prev == he
        @test he.prev.next == he
        he = he.next
      end
    end
  end

  # 2 triangles as a list of half-edges
  h1 = HalfEdge(1, 1)
  h2 = HalfEdge(2, nothing)
  h3 = HalfEdge(2, 1)
  h4 = HalfEdge(3, 2)
  h5 = HalfEdge(3, 1)
  h6 = HalfEdge(1, nothing)
  h7 = HalfEdge(2, 2)
  h8 = HalfEdge(4, nothing)
  h9 = HalfEdge(4, 2)
  h10 = HalfEdge(3, nothing)
  h1.half = h2
  h2.half = h1
  h3.half = h4
  h4.half = h3
  h5.half = h6
  h6.half = h5
  h7.half = h8
  h8.half = h7
  h9.half = h10
  h10.half = h9
  h1.prev = h5
  h1.next = h3
  h3.prev = h1
  h3.next = h5
  h4.prev = h9
  h4.next = h7
  h5.prev = h3
  h5.next = h1
  h7.prev = h4
  h7.next = h9
  h9.prev = h7
  h9.next = h4
  halves = [(h1, h2), (h3, h4), (h5, h6), (h7, h8), (h9, h10)]
  struc = HalfEdgeTopology(halves)
  @test half4elem(struc, 1) == h1
  @test half4elem(struc, 2) == h4
  @test half4vert(struc, 1) == h1
  @test half4vert(struc, 2) == h3
  @test half4vert(struc, 3) == h4
  @test half4vert(struc, 4) == h9
  @test edge4pair(struc, (1, 2)) == 1
  @test edge4pair(struc, (2, 1)) == 1
  @test edge4pair(struc, (2, 3)) == 2
  @test edge4pair(struc, (3, 2)) == 2
  @test edge4pair(struc, (3, 1)) == 3
  @test edge4pair(struc, (1, 3)) == 3
  @test edge4pair(struc, (2, 4)) == 4
  @test edge4pair(struc, (4, 2)) == 4
  @test edge4pair(struc, (4, 3)) == 5
  @test edge4pair(struc, (3, 4)) == 5

  # 2 triangles
  elems = connect.([(1, 2, 3), (4, 3, 2)])
  t = HalfEdgeTopology(elems)
  @test paramdim(t) == 2
  @test nelements(t) == 2
  @test nfacets(t) == 5
  @test nvertices(t) == 4
  @test nfaces(t, 2) == 2
  @test nfaces(t, 1) == 5
  @test nfaces(t, 0) == 4
  test_halfedge(elems, t)

  # 2 triangles + 2 quadrangles
  elems = connect.([(1, 2, 6, 5), (2, 4, 6), (4, 3, 5, 6), (1, 5, 3)])
  t = HalfEdgeTopology(elems)
  @test paramdim(t) == 2
  @test nelements(t) == 4
  @test nfacets(t) == 9
  @test nvertices(t) == 6
  @test nfaces(t, 2) == 4
  @test nfaces(t, 1) == 9
  @test nfaces(t, 0) == 6
  test_halfedge(elems, t)

  # 1 triangle + 3 quadrangles + 1 triangle hole
  elems = connect.([(1, 2, 6, 5), (2, 4, 7, 6), (4, 3, 7), (3, 1, 5, 7)])
  t = HalfEdgeTopology(elems)
  @test paramdim(t) == 2
  @test nelements(t) == 4
  @test nfacets(t) == 11
  @test nvertices(t) == 7
  @test nfaces(t, 2) == 4
  @test nfaces(t, 1) == 11
  @test nfaces(t, 0) == 7
  @test vertices(t) == 1:7
  @test vertex(t, 1) == 1
  @test vertex(t, 7) == 7
  test_halfedge(elems, t)

  # no need to sort elements with consistent orientation
  elems = connect.([(1, 2, 6, 5), (2, 4, 7, 6), (4, 3, 7), (3, 1, 5, 7)])
  t = HalfEdgeTopology(elems, sort=false)
  @test paramdim(t) == 2
  @test nelements(t) == 4
  @test nfacets(t) == 11
  @test nvertices(t) == 7
  @test nfaces(t, 2) == 4
  @test nfaces(t, 1) == 11
  @test nfaces(t, 0) == 7
  test_halfedge(elems, t)

  # correct construction from inconsistent orientation
  e = connect.([(1, 2, 3), (3, 4, 2), (4, 3, 5), (6, 3, 1)])
  t = HalfEdgeTopology(e)
  test_halfedge(e, t)
  n = collect(elements(t))
  @test n[1] == e[1]
  @test n[2] != e[2]
  @test n[3] != e[3]
  @test n[4] != e[4]

  # more challenging case with inconsistent orientation
  e = connect.([(4, 1, 5), (2, 6, 4), (3, 5, 6), (4, 5, 6)])
  t = HalfEdgeTopology(e)
  test_halfedge(e, t)
  n = collect(elements(t))
  @test n == connect.([(4, 1, 5), (4, 6, 2), (6, 5, 3), (4, 5, 6)])

  e = connect.([(1, 2, 3), (1, 3, 4), (2, 5, 3), (5, 4, 6), (3, 5, 4)])
  t = HalfEdgeTopology(e)
  test_halfedge(e, t)

  # indexable api
  g = GridTopology(10, 10)
  t = convert(HalfEdgeTopology, g)
  @test t[begin] == connect((1, 2, 13, 12), Quadrangle)
  @test t[end] == connect((120, 109, 110, 121), Quadrangle)
  @test t[10] == connect((22, 21, 10, 11), Quadrangle)
  @test length(t) == 100
  @test eltype(t) == Connectivity{Quadrangle,4}
  for e in t
    @test e isa Connectivity{Quadrangle,4}
  end
end

@testitem "SimpleTopology" setup = [Setup] begin
  # 2 triangles
  elems = connect.([(1, 2, 3), (4, 3, 2)])
  t = SimpleTopology(elems)
  @test paramdim(t) == 2
  @test connec4elem(t, 1) == (1, 2, 3)
  @test connec4elem(t, 2) == (4, 3, 2)
  @test nvertices(t) == 4
  @test nelements(t) == 2
  @test vertices(t) == 1:4
  @test vertex(t, 1) == 1
  @test vertex(t, 4) == 4
  @test nfaces(t, 2) == 2
  @test nfaces(t, 1) == 0
  @test nfaces(t, 0) == 4

  # 2 triangles + 2 quadrangles
  elems = connect.([(1, 2, 6, 5), (2, 4, 6), (4, 3, 5, 6), (1, 5, 3)])
  t = SimpleTopology(elems)
  @test connec4elem(t, 1) == (1, 2, 6, 5)
  @test connec4elem(t, 2) == (2, 4, 6)
  @test connec4elem(t, 3) == (4, 3, 5, 6)
  @test connec4elem(t, 4) == (1, 5, 3)
  @test nelements(t) == 4
  @test nfacets(t) == 0
  @test nvertices(t) == 6
  @test nfaces(t, 2) == 4
  @test nfaces(t, 1) == 0
  @test nfaces(t, 0) == 6

  # 1 triangle + 3 quadrangles + 1 triangle hole
  elems = connect.([(1, 2, 6, 5), (2, 4, 7, 6), (4, 3, 7), (3, 1, 5, 7)])
  t = SimpleTopology(elems)
  @test connec4elem(t, 1) == (1, 2, 6, 5)
  @test connec4elem(t, 2) == (2, 4, 7, 6)
  @test connec4elem(t, 3) == (4, 3, 7)
  @test connec4elem(t, 4) == (3, 1, 5, 7)
  @test nelements(t) == 4
  @test nfacets(t) == 0
  @test nvertices(t) == 7
  @test nfaces(t, 2) == 4
  @test nfaces(t, 1) == 0
  @test nfaces(t, 0) == 7

  # convert from other topologies
  g = GridTopology(2, 2)
  t = convert(SimpleTopology, g)
  @test nelements(t) == 4
  @test nfacets(t) == 12
  @test nvertices(t) == 9
  @test nfaces(t, 2) == 4
  @test nfaces(t, 1) == 12
  @test nfaces(t, 0) == 9

  # indexable api
  g = GridTopology(10, 10)
  t = convert(SimpleTopology, g)
  @test t[begin] == connect((1, 2, 13, 12), Quadrangle)
  @test t[end] == connect((109, 110, 121, 120), Quadrangle)
  @test t[10] == connect((10, 11, 22, 21), Quadrangle)
  @test length(t) == 100
  @test eltype(t) == Connectivity{Quadrangle,4}
  for e in t
    @test e isa Connectivity{Quadrangle,4}
  end
end
