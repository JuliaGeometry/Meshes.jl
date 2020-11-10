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
    @test nvertices(first(rings(poly))) == 31
    @test !hasholes(poly)
    @test issimple(poly)
    @test orientation(poly) == :CCW
  end

  # COMMAND USED TO GENERATE TEST FILES (VARY --seed = 1, 2, ..., 5)
  # rpg --cluster 30 --algo 2opt --format line --seed 1 --output smooth1 --smooth 2
  fnames = ["smooth$i.line" for i in 1:5]
  polys2 = [readpoly(joinpath(datadir, fname)) for fname in fnames]
  for poly in polys2
    @test nvertices(first(rings(poly))) == 121
    @test !hasholes(poly)
    @test issimple(poly)
    @test orientation(poly) == :CCW
  end

  # COMMAND USED TO GENERATE TEST FILES (VARY --seed = 1, 2, ..., 5)
  # rpg --cluster 30 --algo 2opt --format line --seed 1 --output hole1 --holes 2
  fnames = ["hole$i.line" for i in 1:5]
  polys3 = [readpoly(joinpath(datadir, fname)) for fname in fnames]
  for poly in polys3
    outer, inners = rings(poly)
    @test nvertices(outer) < 31
    @test all(nvertices.(inners) .< 18)
    @test hasholes(poly)
    @test !issimple(poly)
    oorient, iorients = orientation(poly)
    @test oorient == :CCW
    @test all(iorients .== :CW)
  end
end
