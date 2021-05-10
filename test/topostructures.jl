@testset "TopologicalStructure" begin
  @testset "FullStructure" begin
    # 2 triangles
    elems = connect.([(1,2,3),(4,3,2)])
    struc = FullStructure(elems)
    @test nvertices(struc) == 4
    @test nelements(struc) == 2

    # 2 triangles + 2 quadrangles
    elems = connect.([(1,2,6,5),(2,4,6),(4,3,5,6),(1,5,3)])
    struc = FullStructure(elems)
    @test nvertices(struc) == 6
    @test nelements(struc) == 4

    # 1 triangle + 3 quadrangles + 1 triangle hole
    elems = connect.([(1,2,6,5),(2,4,7,6),(4,3,7),(3,1,5,7)])
    struc = FullStructure(elems)
    @test nvertices(struc) == 7
    @test nelements(struc) == 4
  end

  @testset "HalfEdgeStructure" begin
    function test_halfedge(elems, structure)
      @test nelements(structure) == length(elems)
      for e in 1:nelements(structure)
        he = half4elem(e, structure)
        inds = indices(elems[e])
        @test he.elem == e
        @test he.head ∈ inds
      end
    end

    # 2 triangles as a list of half-edges
    h1  = HalfEdge(1, 1)
    h2  = HalfEdge(2, nothing)
    h3  = HalfEdge(2, 1)
    h4  = HalfEdge(3, 2)
    h5  = HalfEdge(3, 1)
    h6  = HalfEdge(1, nothing)
    h7  = HalfEdge(2, 2)
    h8  = HalfEdge(4, nothing)
    h9  = HalfEdge(4, 2)
    h10 = HalfEdge(3, nothing)
    h1.half =  h2; h2.half  = h1
    h3.half =  h4; h4.half  = h3
    h5.half =  h6; h6.half  = h5
    h7.half =  h8; h8.half  = h7
    h9.half = h10; h10.half = h9
    h1.prev = h5;  h1.next = h3
    h3.prev = h1;  h3.next = h5
    h4.prev = h9;  h4.next = h7
    h5.prev = h3;  h5.next = h1
    h7.prev = h4;  h7.next = h9
    h9.prev = h7;  h9.next = h4
    halves  = [(h1,h2),(h3,h4),(h5,h6),(h7,h8),(h9,h10)]
    struc = HalfEdgeStructure(halves)
    @test half4elem(1, struc) == h1
    @test half4elem(2, struc) == h4
    @test half4vert(1, struc) == h1
    @test half4vert(2, struc) == h3
    @test half4vert(3, struc) == h4
    @test half4vert(4, struc) == h9
    @test edge4pair((1,2), struc) == 1
    @test edge4pair((2,1), struc) == 1
    @test edge4pair((2,3), struc) == 2
    @test edge4pair((3,2), struc) == 2
    @test edge4pair((3,1), struc) == 3
    @test edge4pair((1,3), struc) == 3
    @test edge4pair((2,4), struc) == 4
    @test edge4pair((4,2), struc) == 4
    @test edge4pair((4,3), struc) == 5
    @test edge4pair((3,4), struc) == 5

    # 2 triangles
    elems = connect.([(1,2,3),(4,3,2)])
    struc = HalfEdgeStructure(elems)
    @test nvertices(struc) == 4
    @test nelements(struc) == 2
    @test nfacets(struc) == 5
    test_halfedge(elems, struc)

    # 2 triangles + 2 quadrangles
    elems = connect.([(1,2,6,5),(2,4,6),(4,3,5,6),(1,5,3)])
    struc = HalfEdgeStructure(elems)
    @test nvertices(struc) == 6
    @test nelements(struc) == 4
    @test nfacets(struc) == 9
    test_halfedge(elems, struc)

    # 1 triangle + 3 quadrangles + 1 triangle hole
    elems = connect.([(1,2,6,5),(2,4,7,6),(4,3,7),(3,1,5,7)])
    struc = HalfEdgeStructure(elems)
    @test nvertices(struc) == 7
    @test nelements(struc) == 4
    @test nfacets(struc) == 11
    test_halfedge(elems, struc)
  end
end
