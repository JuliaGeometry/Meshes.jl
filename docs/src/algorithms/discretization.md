# Discretization

```@docs
discretize
DiscretizationMethod
FIST
```
# FIST  
Examples for the FIST algorithm. 

### Usage  
  
```@example FIST
using Meshes
using Plots; gr(size = (500,300))

poly = PolyArea(Point2f[(0.22926679f0, 0.47329807f0), (0.45326096f0, 0.2665109f0), (0.51871038f0, 0.22148979f0), (0.3152017f0, 0.31646582f0), (0.5953296f0, 0.62900037f0), (0.5951828f0, 0.6113712f0), (0.74374446f0, 0.72971905f0), (0.29018405f0, 0.5538437f0), (0.25439468f0, 0.5678829f0), (0.22926679f0, 0.47329807f0)], [])
mesh = discretize(poly, FIST())
plot(plot(poly), plot(mesh))
```
###### Examples of polygons and their corresponding FIST-discretized triangular meshes:  
```@eval 
using Plots
using Meshes
using Markdown
gr(size = (500,300))
function readpoly(T, fname) 
  open(fname, "r") do f
    # read outer chain
    n = parse(Int, readline(f))
    outer = map(1:n) do _
      coords = readline(f)
      x, y = parse.(T, split(coords))
      Point(x, y)
    end

    # read inner chains
    inners = []
    while !eof(f)
      n = parse(Int, readline(f))
      inner = map(1:n) do _
        coords = readline(f)
        x, y = parse.(T, split(coords))
        Point(x, y)
      end
      push!(inners, inner)
    end

    # return polygonal area
    PolyArea(outer, inners)
  end
end
poly = readpoly(Float32, "../../../test/data/smooth1.line")
mesh = discretize(poly, FIST())
plot(plot(poly), plot(mesh))
savefig("plot1.svg")
Markdown.parse("![Plot](plot1.svg)")
```
```@eval
using Plots
using Meshes
using Markdown
gr(size = (500,300))
function readpoly(T, fname) 
  open(fname, "r") do f
    # read outer chain
    n = parse(Int, readline(f))
    outer = map(1:n) do _
      coords = readline(f)
      x, y = parse.(T, split(coords))
      Point(x, y)
    end

    # read inner chains
    inners = []
    while !eof(f)
      n = parse(Int, readline(f))
      inner = map(1:n) do _
        coords = readline(f)
        x, y = parse.(T, split(coords))
        Point(x, y)
      end
      push!(inners, inner)
    end

    # return polygonal area
    PolyArea(outer, inners)
  end
end
poly = readpoly(Float32, "../../../test/data/poly1.line")
mesh = discretize(poly, FIST())
plot(plot(poly), plot(mesh))
savefig("plot2.svg")
Markdown.parse("![Plot](plot2.svg)")
```
```@eval
using Plots
using Meshes
using Markdown
gr(size = (500,300))
function readpoly(T, fname) 
  open(fname, "r") do f
    # read outer chain
    n = parse(Int, readline(f))
    outer = map(1:n) do _
      coords = readline(f)
      x, y = parse.(T, split(coords))
      Point(x, y)
    end

    # read inner chains
    inners = []
    while !eof(f)
      n = parse(Int, readline(f))
      inner = map(1:n) do _
        coords = readline(f)
        x, y = parse.(T, split(coords))
        Point(x, y)
      end
      push!(inners, inner)
    end

    # return polygonal area
    PolyArea(outer, inners)
  end
end
poly = readpoly(Float32, "../../../test/data/hole1.line")
mesh = discretize(poly, FIST())
plot(plot(poly), plot(mesh))
savefig("plot3.svg")
Markdown.parse("![Plot](plot3.svg)")
```
```@eval
using Plots
using Meshes
using Markdown
gr(size = (500,300))
function readpoly(T, fname) 
  open(fname, "r") do f
    # read outer chain
    n = parse(Int, readline(f))
    outer = map(1:n) do _
      coords = readline(f)
      x, y = parse.(T, split(coords))
      Point(x, y)
    end

    # read inner chains
    inners = []
    while !eof(f)
      n = parse(Int, readline(f))
      inner = map(1:n) do _
        coords = readline(f)
        x, y = parse.(T, split(coords))
        Point(x, y)
      end
      push!(inners, inner)
    end

    # return polygonal area
    PolyArea(outer, inners)
  end
end
poly = readpoly(Float32, "../../../test/data/poly3.line")
mesh = discretize(poly, FIST())
plot(plot(poly), plot(mesh))
savefig("plot4.svg")
Markdown.parse("![Plot](plot4.svg)")
```