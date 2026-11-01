@testitem "BVHTree" setup = [Setup] begin
  # leaf size greater than or equal to the number of elements
  boxes = [Box(cart(0, 0), cart(1, 1)), Box(cart(2, 2), cart(3, 3)), Box(cart(4, 4), cart(5, 5))]
  nboxes = length(boxes)
  for leafsize in (nboxes, nboxes + 1, 2nboxes)
    ltree = Meshes.BVHTree(boxes; leafsize)
    lroot = ltree.nodes[1]
    @test length(ltree.nodes) == 1
    @test Meshes.isleaf(lroot)
    @test lroot.left == -1
    @test lroot.right == nboxes
    @test sort(ltree.perm) == collect(1:nboxes)
  end

  # construction and traversal retain original element indices
  boxes = reverse([Box(cart(i, 0), cart(i + 1, 1)) for i in 1:40])
  @test_throws ArgumentError Meshes.BVHTree(boxes, leafsize=0)
  @test_throws ArgumentError Meshes.BVHTree(boxes, leafsize=-1)
  for leafsize in (1, 2, 8, 64)
    ltree = @inferred Meshes.BVHTree(boxes; leafsize)
    stack = Int[]
    for query in (Box(cart(5.5, 0), cart(12.5, 1)), boundingbox(boxes), Box(cart(100, 0), cart(101, 1)))
      candidates = Int[]
      Meshes.foreachcandidate(i -> push!(candidates, i), query, ltree, stack)
      @test sort(candidates) == findall(box -> intersects(box, query), boxes)
    end
  end

  # root bounding box
  boxes = reverse([Box(cart(i, 0), cart(i + 1, 1)) for i in 1:40])
  tree = Meshes.BVHTree(boxes)
  @test tree.nodes[1].box ≈ boundingbox(boxes)

  # internal node invariants
  tree = Meshes.BVHTree(boxes, leafsize=1)
  for node in tree.nodes
    if !Meshes.isleaf(node)
      leftbox = tree.nodes[node.left].box
      rightbox = tree.nodes[node.right].box
      @test node.box ≈ Meshes.boundingbox((leftbox, rightbox))
      @test node.left > 0
    end
  end

  # leaf invariants
  tree = Meshes.BVHTree(boxes, leafsize=2)
  for node in tree.nodes
    if Meshes.isleaf(node)
      @test 1 ≤ abs(node.left) ≤ node.right ≤ length(tree.perm)
      @test node.right - abs(node.left) + 1 ≤ tree.leafsize
    end
  end

  # leaf partition completeness
  leafinds = Int[]
  for node in tree.nodes
    if Meshes.isleaf(node)
      append!(leafinds, tree.perm[abs(node.left):node.right])
    end
  end
  @test sort(leafinds) == collect(1:length(boxes))
  @test length(unique(leafinds)) == length(boxes)
end
