"""
    upperhull(pset)

Upper part of the convex hull of a given set of 2D points `pset`. 
The upperhull is computed according to the Andrew's Monotone Chain Algorithm.
See https://en.wikipedia.org/wiki/Convex_hull_algorithms.

The input is sorted on x-coordinate, and the hull is computed in two steps producing the upper and the lower parts of the hull respectively. 
This is a modification of the Graham's Scan algorithm, devised by A. M. Andrew, known as Andrew's Monotone Chain Algorithm. It has the same basic properties as Graham's scan. 
See https://en.wikipedia.org/wiki/Graham_scan
"""
function upperhull(pset::PointSet)
  xs = sort(coordinates.(pset), by = first)
  upperHull = [xs[1],xs[2]]
  for i in 3:length(xs)
    push!(upperHull,xs[i])
    A = Point(upperHull[end-2])
    B = Point(upperHull[end-1])
    C = Point(upperHull[end])
    # remove the second to last point in upperHull while the last 3 points make a "left turn"
    while length(upperHull)>2 && ∠(A,B,C)<0
      splice!(upperHull,lastindex(upperHull)-1)
      if length(upperHull)>2 
        A = Point(upperHull[end-2])
        B = Point(upperHull[end-1])
        C = Point(upperHull[end])
      end
    end
  end 
  upperHull
end 