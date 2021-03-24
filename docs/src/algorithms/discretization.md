# Discretization

```@docs
discretize
DiscretizationMethod
FIST
```
# FIST  
Examples for the FIST algorithm. 

### Usage  
  
```@example overview
#Import required Libraries
using GeoStats
using Plots; gr(size=(600,600))

#Let's write a function to read out polygon
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

#Construct the polygon from the points in 'smooth1.line' 
poly = readpoly(Float32, joinpath("../Meshes/test/data, "smooth1.line"))

#Plot the Polygon
plot(poly)

#Discretize the polygon into Triangles using FIST
mesh = discretize(poly, FIST())

#Plot the Mesh
plot(mesh)
```  





