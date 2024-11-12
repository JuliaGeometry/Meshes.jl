struct VertexItr{T}
  el::T
end

Base.iterate(itr::VertexItr{T}, i=1) where {T<:Polytope} =
  (@inline; (i - 1) % UInt < nvertices(itr.el) % UInt ? (@inbounds vertex(itr.el, i), i + 1) : nothing)

Base.iterate(itr::VertexItr{T}, state=(1, 1)) where {T<:Multi} = begin
  ig, ivg = state # current geometry index, current vertex index in the current geometry
  ig > length(itr.el.geoms) && return nothing

  # iterate through the current geometry
  is = iterate(itr.el.geoms[ig], ivg)

  # start next geometry if current one is done
  isnothing(is) && return iterate(itr, (ig + 1, 1))

  v, ivg = is
  return (v, (ig, ivg))
end

Base.IteratorSize(itr::VertexItr{T}) where {T<:Polytope} = Base.HasLength()
Base.IteratorSize(itr::VertexItr{T}) where {T<:Multi} = Base.HasLength()

Base.length(itr::VertexItr{T}) where {T<:Polytope} = nvertices(T)
Base.length(itr::VertexItr{T}) where {T<:Multi} = sum(nvertices, itr.el.geoms)
Base.eltype(itr::VertexItr{T}) where {T} = Point

verticesiter(e::T) where {T<:Union{Geometry,Multi}} = VertexIter(e)