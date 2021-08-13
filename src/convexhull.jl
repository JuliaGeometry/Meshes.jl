function upperhull(pset::PointSet) 

  xs=Vector{Float64}[] 
  xs=sort(coordinates.(pset), by=first)

  convHull=[xs[1],xs[2]] 
  
  for i in 3:size(xs,1)
    push!(convHull,xs[i])

    # remove the second to last point in convHull while the last 3 points make a "left turn"
    while size(convHull)[1]>2 && ∠(Point(xs[i-2]),Point(xs[i-1]),Point(xs[i]))<0 
      splice!(convHull,i-1)
    end 

  end

  convHull

end