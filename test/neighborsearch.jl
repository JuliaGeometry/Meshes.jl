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

  # non MinkowskiMetric example
  𝒟 = CartesianGrid(T.((0.0, -90.0)), T.((1.0, 1.0)), GridTopology(360, 180))
  s = BallSearch(𝒟, MetricBall(T(150), Haversine(T(6371))))
  n = search(cart(0, 0), s)
  @test Set(n) == Set([32041, 32400, 32401, 32760])

  # construct from vector of geometries
  s = BallSearch([cart(1, 1), cart(2, 2), cart(3, 3)], MetricBall(T(1)))
  @test s isa BallSearch

  # latlon coodinates
  𝒟 = RegularGrid(latlon(0, 0), T.((1.0, 1.0)), GridTopology(10, 10))
  s = BallSearch(𝒟, MetricBall(T(3e5), Haversine()))
  n = search(latlon(0, 0), s)
  @test Set(n) == Set([1, 2, 3, 11, 12, 21])
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

  # latlon coodinates
  𝒟 = RegularGrid(latlon(0, 0), T.((1.0, 1.0)), GridTopology(10, 10))
  s = KNearestSearch(𝒟, 3, metric=Haversine())
  n = search(latlon(0, 0), s)
  @test Set(n) == Set([1, 2, 11])
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

  # latlon coodinates
  𝒟 = RegularGrid(latlon(0, 0), T.((1.0, 1.0)), GridTopology(10, 10))
  s = KBallSearch(𝒟, 10, MetricBall(T(3e5), Haversine()))
  n = search(latlon(5, 5), s)
  @test length(n) == 10
end

@testitem "BoundingBoxSearch" setup = [Setup] begin
  box = [Box(cart(0, 0), cart(1, 1)), Box(cart(2, 2), cart(3, 3)), Box(cart(4, 4), cart(5, 5))]
  domain = GeometrySet(box)

  # basic query
  bvh = BoundingBoxSearch(domain; leafsize=1)
  query = Box(cart(0.5, 0.5), cart(2.5, 2.5))
  answer = findall(box -> intersects(box, query), box)
  result = sort(search(query, bvh))
  @test result == answer

  # no candidates
  bvh = BoundingBoxSearch(domain)
  query = Box(cart(10, 10), cart(11, 11))
  @test isempty(search(query, bvh))

  # different leaf sizes
  query = Box(cart(0.5, 0.5), cart(4.5, 4.5))
  answer = findall(box -> intersects(box, query), box)
  for leafsize in (1, 2, 3, 8)
    leafbvh = BoundingBoxSearch(domain; leafsize)
    @test sort(search(query, leafbvh)) == answer
  end

  # leaf size greater than or equal to the number of elements
  n = nelements(domain)
  for leafsize in (n, n + 1, 2n)
    leafbvh = BoundingBoxSearch(domain; leafsize)
    root = leafbvh.nodes[1]
    @test length(leafbvh.nodes) == 1
    @test Meshes._isleaf(root)
    @test root.first == 1
    @test root.last == n
    @test sort(leafbvh.perm) == collect(1:n)
    @test sort(search(boundingbox(domain), leafbvh)) == collect(1:n)
  end

  # preallocated output
  bvh = BoundingBoxSearch(domain)
  query = Box(cart(0.5, 0.5), cart(1.5, 1.5))
  inds = [100, 200]
  result = search!(inds, query, bvh)
  @test result === inds
  @test inds == [1]

  # invalid leaf size
  @test_throws ArgumentError BoundingBoxSearch(domain; leafsize=0)
  @test_throws ArgumentError BoundingBoxSearch(domain; leafsize=-1)

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
    randbvh = BoundingBoxSearch(randd; leafsize)
    for _ in 1:100
      local query, answer, result
      xmin = rand(rng) * 100
      ymin = rand(rng) * 100
      width = rand(rng) * 20
      height = rand(rng) * 20
      query = Box(cart(xmin, ymin), cart(xmin + width, ymin + height))
      answer = findall(box -> intersects(box, query), randb)
      result = sort(search(query, randbvh))
      @test result == answer
    end
  end

  # root bounding box
  bvh = BoundingBoxSearch(randd)
  @test bvh.nodes[1].box ≈ boundingbox(randd)

  # internal node invariants
  bvh = BoundingBoxSearch(randd; leafsize=1)
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
  bvh = BoundingBoxSearch(randd; leafsize=2)
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
  bvh = BoundingBoxSearch(domain)
  query = Box(cart(0.5, 0.5), cart(2.5, 2.5))
  @inferred search(query, bvh)
  inds = Int[]
  @inferred search!(inds, query, bvh)
end
