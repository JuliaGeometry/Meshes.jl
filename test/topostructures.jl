@testset "Topostructures" begin
  function test_halfedge(elems, structure)
    @test nelements(structure) == length(elems)
    for e in 1:nelements(structure)
      he = edgeonelem(e, structure)
      inds = indices(elems[e])
      @test he.elem == e
      @test he.head ∈ inds
    end
  end

  # 2 triangles
  elems = connect.([(1,2,3),(4,3,2)], Triangle)
  struc = HalfEdgeStructure(elems)
  test_halfedge(elems, struc)

  # 2 triangles + 2 quadrangles
  elems = connect.([(1,2,6,5),(2,4,6),(4,3,5,6),(1,5,3)], Ngon)
  struc = HalfEdgeStructure(elems)
  test_halfedge(elems, struc)

  # 1 triangle + 3 quadrangles + 1 triangle hole
  elems = connect.([(1,2,6,5),(2,4,7,6),(4,3,7),(3,1,5,7)], Ngon)
  struc = HalfEdgeStructure(elems)
  test_halfedge(elems, struc)
end
