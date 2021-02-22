@testset "Polytopes" begin
  @testset "Segment" begin
    s1 = Segment(P2(0,0), P2(1,0))
    s2 = Segment(P2(0.5,0.0), P2(2,0))
    @test s1 ∩ s2 == Segment(P2(0.5,0.0), P2(1,0))
    @test s2 ∩ s1 == Segment(P2(0.5,0.0), P2(1,0))

    s1 = Segment(P2(0,0), P2(1,0))
    s2 = Segment(P2(0,0), P2(0,1))
    @test s1 ∩ s2 == P2(0,0)
    @test s2 ∩ s1 == P2(0,0)

    s1 = Segment(P2(0,0), P2(1,0))
    s2 = Segment(P2(0,0), P2(-1,0))
    @test s1 ∩ s2 == P2(0,0)
    @test s2 ∩ s1 == P2(0,0)

    s1 = Segment(P2(0,0), P2(0,1))
    s2 = Segment(P2(0,0), P2(0,-1))
    @test s1 ∩ s2 == P2(0,0)
    @test s2 ∩ s1 == P2(0,0)

    s1 = Segment(P2(1,1), P2(1,2))
    s2 = Segment(P2(1,1), P2(1,0))
    @test s1 ∩ s2 == P2(1,1)
    @test s2 ∩ s1 == P2(1,1)

    s1 = Segment(P2(1,1), P2(2,1))
    s2 = Segment(P2(1,0), P2(3,0))
    @test s1 ∩ s2 === nothing
    @test s2 ∩ s1 === nothing

    s1 = Segment(P2(0.181429364026879, 0.546811355144474),
                 P2(0.38282226144778, 0.107781953228536))
    s2 = Segment(P2(0.412498700935005, 0.212081819871479),
                 P2(0.395936725690311, 0.252041094122474))
    @test s1 ∩ s2 === nothing
    @test s2 ∩ s1 === nothing

    s1 = Segment(P2(1,2), P2(1,0))
    s2 = Segment(P2(1,0), P2(1,1))
    @test s1 ∩ s2 == Segment(P2(1,1), P2(1,0))

    s1 = Segment(P2(0,0), P2(2,0))
    s2 = Segment(P2(-2,0), P2(-1,0))
    s3 = Segment(P2(-1,0), P2(-2,0))
    @test s1 ∩ s2 === s2 ∩ s1 === nothing
    @test s1 ∩ s3 === s3 ∩ s1 === nothing

    s1 = Segment(P2(-1,0), P2(0,0))
    s2 = Segment(P2(0,0), P2(2,0))
    @test s1 ∩ s2 == s2 ∩ s1 == P2(0,0)

    s1 = Segment(P2(-1,0), P2(1,0))
    s2 = Segment(P2(0,0), P2(3,0))
    @test s1 ∩ s2 == s2 ∩ s1 == Segment(P2(0,0), P2(1,0))

    s1 = Segment(P2(0,0), P2(1,0))
    s2 = Segment(P2(0,0), P2(2,0))
    @test s1 ∩ s2 == s2 ∩ s1 == Segment(P2(0,0), P2(1,0))

    s1 = Segment(P2(0,0), P2(3,0))
    s2 = Segment(P2(1,0), P2(2,0))
    @test s1 ∩ s2 == s2 ∩ s1 == s2

    s1 = Segment(P2(0,0), P2(2,0))
    s2 = Segment(P2(1,0), P2(2,0))
    @test s1 ∩ s2 == s2 ∩ s1 == s2

    s1 = Segment(P2(0,0), P2(2,0))
    s2 = Segment(P2(1,0), P2(3,0))
    @test s1 ∩ s2 == s2 ∩ s1 == Segment(P2(1,0), P2(2,0))

    s1 = Segment(P2(0,0), P2(2,0))
    s2 = Segment(P2(2,0), P2(3,0))
    @test s1 ∩ s2 == s2 ∩ s1 == P2(2,0)

    s1 = Segment(P2(0,0), P2(2,0))
    s2 = Segment(P2(3,0), P2(4,0))
    @test s1 ∩ s2 === s2 ∩ s1 === nothing

    s1 = Segment(P2(2,1), P2(1,2))
    s2 = Segment(P2(1,0), P2(1,1))
    @test s1 ∩ s2 === s2 ∩ s1 === nothing
    
    s1 = Segment(P2(1.5,1.5), P2(3.0,1.5))
    s2 = Segment(P2(3.0,1.0), P2(2.0,2.0))
    @test s1 ∩ s2 == s2 ∩ s1 == P2(2.5,1.5)
  end

  @testset "Triangles" begin
    t = Triangle(P2(0,0), P2(1,0), P2(0,1))
    @test signarea(t) == T(0.5)
    @test area(t) == T(0.5)

    t = Triangle(P2(0,0), P2(0,1), P2(1,0))
    @test signarea(t) == T(-0.5)
    @test area(t) == T(0.5)

    t = Triangle(P2(0,0), P2(1,0), P2(1,1))
    for p in P2[(0,0),(1,0),(1,1),(0.5,0.),(1.,0.5),(0.5,0.5)]
      @test p ∈ t
    end
    for p in P2[(-1,0),(0,-1),(0.5,1.)]
      @test p ∉ t
    end
  end

  @testset "Quadrangles" begin
    q = Quadrangle(P2(0,0), P2(1,0), P2(1,1), P2(0,1))
    @test area(q) == T(1)

    q = Quadrangle(P2(0,0), P2(1,0), P2(1.5,1.0), P2(0.5,1.0))
    @test area(q) == T(1)

    q = Quadrangle(P2(0,0), P2(1,0), P2(1.5,1.0), P2(0.5,1.0))
    for p in P2[(0,0),(1,0),(1.5,1.0),(0.5,1.0),(0.5,0.5)]
      @test p ∈ q
    end
    for p in P2[(0,1),(1.5,0.0)]
      @test p ∉ q
    end
  end

  @testset "Chains" begin
    # unique vertices
    c = Chain(P2[(1,1),(2,2),(2,2),(3,3)])
    @test unique(c) == Chain(P2[(1,1),(2,2),(3,3)])
    @test c == Chain(P2[(1,1),(2,2),(2,2),(3,3)])
    unique!(c)
    @test c == Chain(P2[(1,1),(2,2),(3,3)])

    # closing/opening chains
    c = Chain(P2[(1,1),(2,2),(3,3)])
    close!(c)
    @test c == Chain(P2[(1,1),(2,2),(3,3),(1,1)])
    open!(c)
    @test c == Chain(P2[(1,1),(2,2),(3,3)])

    # reversing chains
    c = Chain(P2[(1,1),(2,2),(3,3)])
    reverse!(c)
    @test c == Chain(P2[(3,3),(2,2),(1,1)])
    c = Chain(P2[(1,1),(2,2),(3,3)])
    @test reverse(c) == Chain(P2[(3,3),(2,2),(1,1)])

    # angles and inner angles
    c = Chain(P2[(0,0),(1,0),(1,1),(0,1),(0,0)])
    @test angles(c) ≈ [-π/2, -π/2, -π/2, -π/2]
    c = Chain(P2[(0,0),(1,0),(1,1),(0,1)])
    @test angles(c) ≈ [-π/2, -π/2]
    c = Chain(P2[(0,0),(1,0),(1,1),(2,1),(2,2),(1,2),(0,0)])
    @test angles(c) ≈ [-atan(2), -π/2, +π/2, -π/2, -π/2, -(π-atan(2))]
    @test innerangles(c) ≈ [atan(2), π/2, 3π/2, π/2, π/2, π-atan(2)]

    # reconstruct chain from vertices
    c1 = Chain(P2[(0,0),(1,0),(1,1),(0,1),(0,0)])
    c2 = Chain(vertices(c1))
    @test c1 == c2

    # centroid
    c = Chain(P2[(0,0),(1,0),(1,1),(0,1),(0,0)])
    @test centroid(c) == P2(0.5, 0.5)
  end

  @testset "PolyAreas" begin
    # COMMAND USED TO GENERATE TEST FILES (VARY --seed = 1, 2, ..., 5)
    # rpg --cluster 30 --algo 2opt --format line --seed 1 --output poly1
    fnames = ["poly$i.line" for i in 1:5]
    polys1 = [readpoly(joinpath(datadir, fname)) for fname in fnames]
    for poly in polys1
      @test nvertices(first(chains(poly))) == 30
      @test !hasholes(poly)
      @test issimple(poly)
      @test orientation(poly) == [:CCW]
      @test unique(poly) == poly
    end

    # COMMAND USED TO GENERATE TEST FILES (VARY --seed = 1, 2, ..., 5)
    # rpg --cluster 30 --algo 2opt --format line --seed 1 --output smooth1 --smooth 2
    fnames = ["smooth$i.line" for i in 1:5]
    polys2 = [readpoly(joinpath(datadir, fname)) for fname in fnames]
    for poly in polys2
      @test nvertices(first(chains(poly))) == 120
      @test !hasholes(poly)
      @test issimple(poly)
      @test orientation(poly) == [:CCW]
      @test unique(poly) == poly
    end

    # COMMAND USED TO GENERATE TEST FILES (VARY --seed = 1, 2, ..., 5)
    # rpg --cluster 30 --algo 2opt --format line --seed 1 --output hole1 --holes 2
    fnames = ["hole$i.line" for i in 1:5]
    polys3 = [readpoly(joinpath(datadir, fname)) for fname in fnames]
    for poly in polys3
      rings = chains(poly)
      @test nvertices(first(rings)) < 30
      @test all(nvertices.(rings[2:end]) .< 18)
      @test hasholes(poly)
      @test !issimple(poly)
      orients = orientation(poly)
      @test orients[1] == :CCW
      @test all(orients[2:end] .== :CW)
      @test unique(poly) == poly
    end

    # test bridges
    for poly in [polys1; polys2; polys3]
      b  = bridge(poly)
      nb = nvertices(b)
      np = nvertices.(chains(poly))
      @test orientation(b) == :CCW
      @test nb ≥ sum(np)
    end

    # bridges between holes
    outer = P2[(0,0),(1,0),(1,1),(0,1),(0,0)]
    hole1 = P2[T.((0.2,0.2)),T.((0.4,0.2)),T.((0.4,0.4)),T.((0.2,0.4)),T.((0.2,0.2))]
    hole2 = P2[T.((0.6,0.2)),T.((0.8,0.2)),T.((0.8,0.4)),T.((0.6,0.4)),T.((0.6,0.2))]
    poly  = PolyArea(outer, [hole1, hole2])
    chain = bridge(poly)
    target = T[
      0.0  0.2  0.2  0.4  0.4  0.6  0.6  0.8  0.8  0.6  0.4  0.2  0.0  1.0  1.0  0.0
      0.0  0.2  0.4  0.4  0.2  0.2  0.4  0.4  0.2  0.2  0.2  0.2  0.0  0.0  1.0  1.0
    ]
    @test vertices(chain) == Point.(eachcol(target))

    # test uniqueness
    points = P2[(1,1),(2,2),(2,2),(3,3),(1,1)]
    poly   = PolyArea(points)
    unique!(poly)
    @test first(chains(poly)) == Chain(points)

    # unique and bridges
    poly = PolyArea(P2[(0,0),(1,0),(1,0),(1,1),(1,2),(0,2),(0,1),(0,1),(0,0)])
    chain = poly |> unique |> bridge
    @test chain == Chain(P2[(0,0),(1,0),(1,1),(1,2),(0,2),(0,1),(0,0)])

    # centroid
    poly = PolyArea(P2[(0,0),(1,0),(1,1),(0,1),(0,0)])
    @test centroid(poly) == P2(0.5, 0.5)
  end
end
