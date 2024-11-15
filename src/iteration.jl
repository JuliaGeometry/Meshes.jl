struct VertexItr{T}
  el::T
end

eachvertex(e::T) where {T} = error("Vertex iterator not implemented for type $T")

eachvertex(e::Union{Mesh,Polytope,MultiPolytope}) = VertexItr(e)

_v_iterate(el::Polytope, i) =
  (@inline; (i - 1) % UInt < nvertices(el) % UInt ? (@inbounds vertex(el, i), i + 1) : nothing)

_v_iterate(el::Mesh, i) = (@inline; (i - 1) % UInt < nvertices(el) % UInt ? (@inbounds vertex(el, i), i + 1) : nothing)

Base.iterate(itr::VertexItr{<:Mesh}, i=1) = _v_iterate(itr.el, i)

Base.iterate(itr::VertexItr{<:Polytope}, i=1) = _v_iterate(itr.el, i)

@propagate_inbounds Base.iterate(itr::VertexItr{<:Meshes.MultiPolytope}, state=(1, 1)) = begin
  ig, ivg = state
  ig > length(itr.el.geoms) && return nothing

  is = _v_iterate(itr.el.geoms[ig], ivg)
  is === nothing && return Base.iterate(itr, (ig + 1, 1))

  v, ivg = is
  return (v, (ig, ivg))
end

Base.length(itr::VertexItr{<:Mesh}) = nvertices(T)
Base.length(itr::VertexItr{<:Polytope}) = nvertices(T)
Base.length(itr::VertexItr{<:Meshes.MultiPolytope}) = sum(nvertices, itr.el.geoms)

Base.IteratorSize(itr::VertexItr) = Base.HasLength()
Base.eltype(::VertexItr) = Point
