@testitem "BVH" setup = [Setup] begin
  box = [Box(cart(0, 0), cart(1, 1)), Box(cart(2, 2), cart(3, 3)), Box(cart(4, 4), cart(5, 5))]
  domain = GeometrySet(box)

  # basic query
  bvh = BVH(domain; leafsize=1)
  query = Box(cart(0.5, 0.5), cart(2.5, 2.5))
  answer = findall(box -> intersects(box, query), box)
  result = sort(candidates(query, bvh))
  @test result == answer

  # no candidates
  bvh = BVH(domain)
  query = Box(cart(10, 10), cart(11, 11))
  @test isempty(candidates(query, bvh))

  # different leaf sizes
  query = Box(cart(0.5, 0.5), cart(4.5, 4.5))
  answer = findall(box -> intersects(box, query), box)
  for leafsize in (1, 2, 3, 8)
    leafbvh = BVH(domain; leafsize)
    @test sort(candidates(query, leafbvh)) == answer
  end

  # leaf size greater than or equal to the number of elements
  n = nelements(domain)
  for leafsize in (n, n + 1, 2n)
    leafbvh = BVH(domain; leafsize)
    root = leafbvh.nodes[1]
    @test length(leafbvh.nodes) == 1
    @test Meshes._isleaf(root)
    @test root.first == 1
    @test root.last == n
    @test sort(leafbvh.perm) == collect(1:n)
    @test sort(candidates(boundingbox(domain), leafbvh)) == collect(1:n)
  end

  # preallocated output
  bvh = BVH(domain)
  query = Box(cart(0.5, 0.5), cart(1.5, 1.5))
  inds = [100, 200]
  result = candidates!(inds, query, bvh)
  @test result === inds
  @test inds == [1]

  # invalid leaf size
  @test_throws ArgumentError BVH(domain; leafsize=0)
  @test_throws ArgumentError BVH(domain; leafsize=-1)

  # randomized brute-force equivalence
  rng = StableRNG(1234)
  randb = [
    let
      xmin = rand(rng) * 100
      ymin = rand(rng) * 100
      width = rand(rng) * 10
      height = rand(rng) * 10
      Box(cart(xmin, ymin), cart(xmin + width, ymin + height))
    end for _ in 1:200
  ]
  randd = GeometrySet(randb)
  for leafsize in (1, 2, 4, 8, 16)
    randbvh = BVH(randd; leafsize)
    for _ in 1:100
      local query, answer, result
      xmin = rand(rng) * 100
      ymin = rand(rng) * 100
      width = rand(rng) * 20
      height = rand(rng) * 20
      query = Box(cart(xmin, ymin), cart(xmin + width, ymin + height))
      answer = findall(box -> intersects(box, query), randb)
      result = sort(candidates(query, randbvh))
      @test result == answer
    end
  end

  # root bounding box
  bvh = BVH(randd)
  @test bvh.nodes[1].box ≈ boundingbox(randd)

  # internal node invariants
  bvh = BVH(randd; leafsize=1)
  for node in bvh.nodes
    if !Meshes._isleaf(node)
      leftbox = bvh.nodes[node.left].box
      rightbox = bvh.nodes[node.right].box
      @test node.box ≈ Meshes._bboxes((leftbox, rightbox))
      @test node.first == 0
      @test node.last == 0
    end
  end

  # leaf invariants
  bvh = BVH(randd; leafsize=2)
  for node in bvh.nodes
    if Meshes._isleaf(node)
      @test node.left == 0
      @test node.right == 0
      @test 1 ≤ node.first ≤ node.last ≤ length(bvh.perm)
      @test node.last - node.first + 1 ≤ bvh.leafsize
    end
  end

  # leaf partition completeness
  leafinds = Int[]
  for node in bvh.nodes
    if Meshes._isleaf(node)
      append!(leafinds, bvh.perm[node.first:node.last])
    end
  end
  @test sort(leafinds) == collect(1:nelements(randd))
  @test length(unique(leafinds)) == nelements(randd)

  # type stability tests
  bvh = BVH(domain)
  query = Box(cart(0.5, 0.5), cart(2.5, 2.5))
  @inferred candidates(query, bvh)
  inds = Int[]
  @inferred candidates!(inds, query, bvh)
end
