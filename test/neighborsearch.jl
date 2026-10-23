@testitem "BallSearch" setup = [Setup] begin
  𝒟 = CartesianGrid(T.((-0.5, -0.5)), T.((1.0, 1.0)), GridTopology(10, 10))

  s = BallSearch(𝒟, MetricBall(T(1)))
  n = search(cart(0, 0), s)
  @test Set(n) == Set([1, 2, 11])
  n = search(cart(9, 0), s)
  @test Set(n) == Set([9, 10, 20])
  n = search(cart(0, 9), s)
  @test Set(n) == Set([91, 81, 92])
  n = search(cart(9, 9), s)
  @test Set(n) == Set([100, 99, 90])

  s = BallSearch(𝒟, MetricBall(T(√2 + eps(T))))
  n = search(cart(0, 0), s)
  @test Set(n) == Set([1, 2, 11, 12])
  n = search(cart(9, 0), s)
  @test Set(n) == Set([9, 10, 19, 20])
  n = search(cart(0, 9), s)
  @test Set(n) == Set([81, 82, 91, 92])
  n = search(cart(9, 9), s)
  @test Set(n) == Set([89, 90, 99, 100])

  # different units
  s = BallSearch(𝒟, MetricBall(T(10) * u"dm"))
  n = search(Point(T(900) * u"cm", T(900) * u"cm"), s)
  @test Set(n) == Set([100, 99, 90])
  n = search(Point(T(9000) * u"mm", T(9000) * u"mm"), s)
  @test Set(n) == Set([100, 99, 90])

  # construct from vector of geometries
  s = BallSearch([cart(1, 1), cart(2, 2), cart(3, 3)], MetricBall(T(1)))
  @test s isa BallSearch
end

@testitem "KNearestSearch" setup = [Setup] begin
  𝒟 = CartesianGrid(T.((-0.5, -0.5)), T.((1.0, 1.0)), GridTopology(10, 10))
  s = KNearestSearch(𝒟, 3)
  n = search(cart(0, 0), s)
  @test Set(n) == Set([1, 2, 11])
  n = search(cart(9, 0), s)
  @test Set(n) == Set([9, 10, 20])
  n = search(cart(0, 9), s)
  @test Set(n) == Set([91, 81, 92])
  n = search(cart(9, 9), s)
  @test Set(n) == Set([100, 99, 90])
  n, d = searchdists(cart(9, 9), s)
  @test Set(n) == Set([100, 99, 90])
  @test length(d) == 3
  n = Vector{Int}(undef, maxneighbors(s))
  nn = search!(n, cart(9, 9), s)
  @test nn == 3
  @test Set(n[1:nn]) == Set([100, 99, 90])
  n = Vector{Int}(undef, maxneighbors(s))
  d = Vector{ℳ}(undef, maxneighbors(s))
  nn = searchdists!(n, d, cart(9, 9), s)
  @test nn == 3
  @test Set(n[1:nn]) == Set([100, 99, 90])

  # different units
  s = KNearestSearch(𝒟, 3)
  n = search(Point(T(900) * u"cm", T(900) * u"cm"), s)
  @test Set(n) == Set([100, 99, 90])
  n = search(Point(T(9000) * u"mm", T(9000) * u"mm"), s)
  @test Set(n) == Set([100, 99, 90])

  # construct from vector of geometries
  s = KNearestSearch([cart(1, 1), cart(2, 2), cart(3, 3)], 3)
  @test s isa KNearestSearch
end

@testitem "KBallSearch" setup = [Setup] begin
  𝒟 = CartesianGrid(T.((-0.5, -0.5)), T.((1.0, 1.0)), GridTopology(10, 10))

  s = KBallSearch(𝒟, 10, MetricBall(T(100)))
  n = search(cart(5, 5), s)
  @test length(n) == 10

  s = KBallSearch(𝒟, 10, MetricBall(T.((100, 100))))
  n = search(cart(5, 5), s)
  @test length(n) == 10

  s = KBallSearch(𝒟, 10, MetricBall(T(1)))
  n = search(cart(5, 5), s)
  @test length(n) == 5
  @test n[1] == 56

  s = KBallSearch(𝒟, 10, MetricBall(T(1)))
  n, d = searchdists(cart(5, 5), s)
  @test length(n) == 5
  @test length(d) == 5

  s = KBallSearch(𝒟, 10, MetricBall(T(1)))
  n = Vector{Int}(undef, maxneighbors(s))
  nn = search!(n, cart(5, 5), s)
  @test nn == 5

  s = KBallSearch(𝒟, 10, MetricBall(T(1)))
  n = Vector{Int}(undef, maxneighbors(s))
  d = Vector{ℳ}(undef, maxneighbors(s))
  nn = searchdists!(n, d, cart(5, 5), s)
  @test nn == 5

  mask = trues(nelements(𝒟))
  mask[56] = false
  n = search(cart(5, 5), s, mask=mask)
  @test length(n) == 4
  n = search(cart(-0.2, -0.2), s)
  @test length(n) == 1
  n = search(cart(-10, -10), s)
  @test length(n) == 0
  n, d = searchdists(cart(5, 5), s, mask=mask)
  @test length(n) == 4
  @test length(d) == 4

  # different units
  s = KBallSearch(𝒟, 10, MetricBall(T(10) * u"dm"))
  n = search(Point(T(500) * u"cm", T(500) * u"cm"), s)
  @test Set(n) == Set([56, 66, 55, 57, 46])
  n = search(Point(T(5000) * u"mm", T(5000) * u"mm"), s)
  @test Set(n) == Set([56, 66, 55, 57, 46])

  # construct from vector of geometries
  s = KBallSearch([cart(1, 1), cart(2, 2), cart(3, 3)], 3, MetricBall(T(1)))
  @test s isa KBallSearch
end

@testitem "BBoxSearch" setup = [Setup] begin
  # basic query
  boxes = [Box(cart(0, 0), cart(1, 1)), Box(cart(2, 2), cart(3, 3)), Box(cart(4, 4), cart(5, 5))]
  searcher = BBoxSearch(boxes, 1)
  query = Box(cart(0.5, 0.5), cart(2.5, 2.5))
  answer = findall(box -> intersects(box, query), boxes)
  result = sort(search(query, searcher))
  @test result == answer

  # constructor paths with default and explicit leaf sizes
  geoms = [Triangle(cart(i, i), cart(i + 1, i), cart(i, i + 1)) for i in (0, 2, 4)]
  boxes = boundingbox.(geoms)
  query = Box(cart(0.5, 0.5), cart(2.5, 2.5))
  for input in (boxes, GeometrySet(boxes), geoms, GeometrySet(geoms))
    @test sort(search(query, BBoxSearch(input))) == [1, 2]
    @test sort(search(query, BBoxSearch(input, 1))) == [1, 2]
  end

  # searchers retain their original boxes when the caller changes its vector
  for leafsize in (1, 8)
    boxes = [Box(cart(0, 0), cart(1, 1)), Box(cart(2, 2), cart(3, 3)), Box(cart(4, 4), cart(5, 5))]
    searchers = (BBoxSearch(boxes, leafsize), BBoxSearch(GeometrySet(boxes), leafsize))
    boxes[1] = Box(cart(10, 10), cart(11, 11))
    for searcher in searchers
      @test sort(search(query, searcher)) == [1, 2]
    end
    reverse!(boxes)
    push!(boxes, Box(cart(0, 0), cart(1, 1)))
    for searcher in searchers
      @test sort(search(query, searcher)) == [1, 2]
    end
  end

  # no candidates
  boxes = [Box(cart(0, 0), cart(1, 1)), Box(cart(2, 2), cart(3, 3)), Box(cart(4, 4), cart(5, 5))]
  searcher = BBoxSearch(boxes)
  query = Box(cart(10, 10), cart(11, 11))
  @test isempty(search(query, searcher))

  # different leaf sizes
  boxes = [Box(cart(0, 0), cart(1, 1)), Box(cart(2, 2), cart(3, 3)), Box(cart(4, 4), cart(5, 5))]
  query = Box(cart(0.5, 0.5), cart(4.5, 4.5))
  answer = findall(box -> intersects(box, query), boxes)
  for leafsize in (1, 2, 3, 8)
    leafsearcher = BBoxSearch(boxes, leafsize)
    @test sort(search(query, leafsearcher)) == answer
  end

  # leaf size greater than or equal to the number of elements
  n = length(boxes)
  for leafsize in (n, n + 1, 2n)
    @test sort(search(boundingbox(boxes), BBoxSearch(boxes, leafsize))) == collect(eachindex(boxes))
  end

  # masks preserve original domain indices
  boxes = [Box(cart(0, 0), cart(1, 1)), Box(cart(2, 2), cart(3, 3)), Box(cart(4, 4), cart(5, 5))]
  searcher = BBoxSearch(boxes)
  query = boundingbox(boxes)
  mask = Bool[false, true, false]
  @test search(query, searcher; mask=mask) == [2]
  @test isempty(search(query, searcher; mask=falses(length(boxes))))

  # callback adapter preserves traversal order
  candidates = Int[]
  Meshes._foreachcandidate(i -> push!(candidates, i), query, searcher, Int[])
  @test candidates == search(query, searcher)

  # invalid leaf size
  boxes = [Box(cart(0, 0), cart(1, 1)), Box(cart(2, 2), cart(3, 3)), Box(cart(4, 4), cart(5, 5))]
  @test_throws ArgumentError BBoxSearch(boxes, 0)
  @test_throws ArgumentError BBoxSearch(boxes, -1)

  # randomized brute-force equivalence
  rng = StableRNG(1234)
  boxes = [
    let
      xmin = rand(rng) * 100
      ymin = rand(rng) * 100
      width = rand(rng) * 10
      height = rand(rng) * 10
      Box(cart(xmin, ymin), cart(xmin + width, ymin + height))
    end for _ in 1:200
  ]
  for leafsize in (1, 2, 4, 8, 16)
    randbvh = BBoxSearch(boxes, leafsize)
    for _ in 1:100
      local query, answer, result
      xmin = rand(rng) * 100
      ymin = rand(rng) * 100
      width = rand(rng) * 20
      height = rand(rng) * 20
      query = Box(cart(xmin, ymin), cart(xmin + width, ymin + height))
      answer = findall(box -> intersects(box, query), boxes)
      result = sort(search(query, randbvh))
      @test result == answer
    end
  end

  # type stability tests
  boxes = [Box(cart(0, 0), cart(1, 1)), Box(cart(2, 2), cart(3, 3)), Box(cart(4, 4), cart(5, 5))]
  searcher = BBoxSearch(boxes)
  query = Box(cart(0.5, 0.5), cart(2.5, 2.5))
  @inferred search(query, searcher)
end

@testitem "BVHTree" setup = [Setup] begin
  # leaf size greater than or equal to the number of elements
  boxes = [Box(cart(0, 0), cart(1, 1)), Box(cart(2, 2), cart(3, 3)), Box(cart(4, 4), cart(5, 5))]
  n = length(boxes)
  for leafsize in (n, n + 1, 2n)
    tree = Meshes.BVHTree(boxes, leafsize)
    root = tree.nodes[1]
    @test length(tree.nodes) == 1
    @test Meshes._isleaf(root)
    @test root.left == -1
    @test root.right == n
    @test sort(tree.perm) == collect(1:n)
  end

  # construction and traversal retain original element indices
  boxes = [Box(cart(i, 0), cart(i + 1, 1)) for i in 1:40]
  reverse!(boxes)
  @test_throws ArgumentError Meshes.BVHTree(boxes, 0)
  @test_throws ArgumentError Meshes.BVHTree(boxes, -1)
  for leafsize in (1, 2, 8, 64)
    tree = @inferred Meshes.BVHTree(boxes, leafsize)
    stack = Int[]
    for query in (Box(cart(5.5, 0), cart(12.5, 1)), boundingbox(boxes), Box(cart(100, 0), cart(101, 1)))
      candidates = Int[]
      Meshes._foreachcandidate(i -> push!(candidates, i), query, tree, stack)
      @test sort(candidates) == findall(box -> intersects(box, query), boxes)
    end
  end

  # root bounding box
  tree = Meshes.BVHTree(boxes)
  @test tree.nodes[1].box ≈ boundingbox(boxes)

  # internal node invariants
  tree = Meshes.BVHTree(boxes, 1)
  for node in tree.nodes
    if !Meshes._isleaf(node)
      leftbox = tree.nodes[node.left].box
      rightbox = tree.nodes[node.right].box
      @test node.box ≈ Meshes._bboxes((leftbox, rightbox))
      @test node.left > 0
    end
  end

  # leaf invariants
  tree = Meshes.BVHTree(boxes, 2)
  for node in tree.nodes
    if Meshes._isleaf(node)
      @test 1 ≤ abs(node.left) ≤ node.right ≤ length(tree.perm)
      @test node.right - abs(node.left) + 1 ≤ tree.leafsize
    end
  end

  # leaf partition completeness
  leafinds = Int[]
  for node in tree.nodes
    if Meshes._isleaf(node)
      append!(leafinds, tree.perm[abs(node.left):node.right])
    end
  end
  @test sort(leafinds) == collect(1:length(boxes))
  @test length(unique(leafinds)) == length(boxes)
end
