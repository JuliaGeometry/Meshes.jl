function UpperConvexHull(mySet::PointSet) 
  
    x=Array{Float64, 1}[] # Vector of Array{Float64,2}
   
    for p in mySet 
     push!(x,coordinates(p))   
    end
   
    x=sort(x, by=first)
    
    convHull=[x[1],x[2]] 
    for i in 3:size(x,1)
     push!(convHull,x[i])
         while size(convHull)[1]>2 && ∠(Point(x[i-2]),Point(x[i-1]),Point(x[i]))<0
 
              # os ultimos 3 pontos fazem volta pra esquerda i,i-1,i-2
         #remover i-1 de convHull
        # print(i)
         #print(convHull)
         splice!(convHull,i-1)
 
         end 
    
     end
     convHull
     #=
     conn=Connectivity{Ngon,size(convHull,1)}
     print(size(convHull,1))
     print(conn)
     SimpleMesh(PointSet(convHull),conn)
 =#
 end