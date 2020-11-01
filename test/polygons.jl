@testset "Polygons" begin
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

      Polygon(outer, inners)
    end
  end

  # rpg --cluster 30 --algo 2opt --format line --output rpg
  files1 = ["poly$i.line" for i in 1:5]
  for file in files1
    P = readpoly(joinpath(datadir, file))
  end

  # rpg --cluster 30 --algo 2opt --format line --output rpg --holes 2
  files2 = ["hole$i.line" for i in 1:5]
  for file in files2
    P = readpoly(joinpath(datadir, file))
  end
end
