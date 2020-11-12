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
  end

  @testset "Chains" begin
    c = Chain(P2[(1,1),(2,2),(2,2),(3,3)])
    @test unique(c) == Chain(P2[(1,1),(2,2),(3,3)])
    @test c == Chain(P2[(1,1),(2,2),(2,2),(3,3)])
    unique!(c)
    @test c == Chain(P2[(1,1),(2,2),(3,3)])

    c = Chain(P2[(1,1),(2,2),(3,3)])
    close!(c)
    @test c == Chain(P2[(1,1),(2,2),(3,3),(1,1)])
    open!(c)
    @test c == Chain(P2[(1,1),(2,2),(3,3)])

    c = Chain(P2[(1,1),(2,2),(3,3)])
    reverse!(c)
    @test c == Chain(P2[(3,3),(2,2),(1,1)])
  end

  @testset "PolySurfaces" begin
    # helper function to read *.line files containing polygons
    # generated with RPG (https://github.com/cgalab/genpoly-rpg)
    function readpoly(fname)
      open(fname, "r") do f
        # read outer chain
        n = parse(Int, readline(f))
        outer = map(1:n) do _
          coords = readline(f)
          x, y = parse.(Float64, split(coords))
          Point(x, y)
        end

        # read inner chains
        inners = []
        while !eof(f)
          n = parse(Int, readline(f))
          inner = map(1:n) do _
            coords = readline(f)
            x, y = parse.(Float64, split(coords))
            Point(x, y)
          end
          push!(inners, inner)
        end

        PolySurface(outer, inners)
      end
    end

    # COMMAND USED TO GENERATE TEST FILES (VARY --seed = 1, 2, ..., 5)
    # rpg --cluster 30 --algo 2opt --format line --seed 1 --output poly1
    fnames = ["poly$i.line" for i in 1:5]
    polys1 = [readpoly(joinpath(datadir, fname)) for fname in fnames]
    for poly in polys1
      @test nvertices(first(chains(poly))) == 31
      @test !hasholes(poly)
      @test issimple(poly)
      @test orientation(poly) == :CCW
      @test unique(poly) == poly
    end

    # COMMAND USED TO GENERATE TEST FILES (VARY --seed = 1, 2, ..., 5)
    # rpg --cluster 30 --algo 2opt --format line --seed 1 --output smooth1 --smooth 2
    fnames = ["smooth$i.line" for i in 1:5]
    polys2 = [readpoly(joinpath(datadir, fname)) for fname in fnames]
    for poly in polys2
      @test nvertices(first(chains(poly))) == 121
      @test !hasholes(poly)
      @test issimple(poly)
      @test orientation(poly) == :CCW
      @test unique(poly) == poly
    end

    # COMMAND USED TO GENERATE TEST FILES (VARY --seed = 1, 2, ..., 5)
    # rpg --cluster 30 --algo 2opt --format line --seed 1 --output hole1 --holes 2
    fnames = ["hole$i.line" for i in 1:5]
    polys3 = [readpoly(joinpath(datadir, fname)) for fname in fnames]
    for poly in polys3
      outer, inners = chains(poly)
      @test nvertices(outer) < 31
      @test all(nvertices.(inners) .< 18)
      @test hasholes(poly)
      @test !issimple(poly)
      oorient, iorients = orientation(poly)
      @test oorient == :CCW
      @test all(iorients .== :CW)
      @test unique(poly) == poly
    end

    c = Chain(P2[(1,1),(2,2),(2,2),(3,3),(1,1)])
    p = PolySurface(c)
    unique!(p)
    outer, inners = chains(p)
    @test outer == Chain(P2[(1,1),(2,2),(3,3),(1,1)])
  end
end
