"""
Compute the upper part of the convex hull of a given set of 2D points 'pset' using the "Graham Scan" algorithm.

*) Sort the points lexicographic

*) Assess the last 3 points in 'upperHull'. If the the clockwise angle is <180° remove the second to last point.

"""

function upperhull(pset::PointSet) 

  xs=sort(coordinates.(pset), by=first)

  upperHull=[xs[1],xs[2]] 
  
  for i in 3:length(xs)
    push!(upperHull,xs[i])
    
    # remove the second to last point in upperHull while the last 3 points make a "left turn"
    while length(upperHull)>2 && ∠(Point(upperHull[lastindex(upperHull)-2]),Point(upperHull[lastindex(upperHull)-1]),Point(upperHull[lastindex(upperHull)]))<0 
      
      splice!(upperHull,lastindex(upperHull)-1)
    
    end 

  end

  upperHull

end