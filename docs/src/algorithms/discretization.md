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
# Import required Libraries
using GeoStats
using Plots; gr(size=(300,300))

# Function to read Polygon from given values
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
```  
Now, we will read point values from an example file `smooth1.line`, and construct a polygon from those files using the `readpoly` function in the following manner:  
```@example overview
# Construct the polygon from the points in 'smooth1.line' 
poly = readpoly(Float32, "/home/atreyamaj/.julia/dev/Meshes/test/data/smooth1.line")
plot(poly)
```  
We will now use the `FIST` algorithm to discretize the polygon we previously obtained into triangles, and plot the discretized mesh in the following manner:  
```@example overview
# Discretize the polygon into Triangles using FIST
mesh = discretize(poly, FIST())
plot(mesh)
```
