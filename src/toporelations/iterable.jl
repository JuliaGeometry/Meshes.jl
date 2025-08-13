"""
    <:AbstractRelation{Q,P}

Read as "\$typename Q's of P"

See also: [`Bounding`](@ref), [`Shared`](@ref), [`Adjacent`](@ref)
"""
abstract type AbstractRelation{Q,P} end

(T::Type{<:AbstractRelation{Q,P} where {Q,P}})(msh::HalfEdgeMesh) = T(topology(msh))

struct TopologyIterator{T<:AbstractRelation}
    relation::T
    e::HalfEdge
end

Base.IteratorSize(::Type{<:TopologyIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:TopologyIterator}) = Int

struct ParametricDimension{N} end

const VERTEX = ParametricDimension{0}
const SEGMENT = ParametricDimension{1}
const FACE = ParametricDimension{2}
const CELL = ParametricDimension{3}

_dim(::Type{ParametricDimension{N}}) where N = N

struct Bounding{Q<:ParametricDimension,P<:ParametricDimension,T<:HalfEdgeTopology} <: AbstractRelation{Q,P}
    t::T

    function Bounding{Q,P}(t::T) where {Q<:ParametricDimension,P<:ParametricDimension,T<:HalfEdgeTopology}
        _dim(P) ≤ paramdim(t) || throw(DomainError((P,T), "topology with rank $(ParametricDimension{paramdim(t)}) has no $P. `P` must be ≤ 2"))
        _dim(Q) < _dim(P) || throw(DomainError((Q,P), "cannot calculate bounding $Q's of $P"))
        P === VERTEX && throw(DomainError(P, "the boundary of a VERTEX is not well-defined"))

        return new{Q,P,T}(t)
    end
end

function (B::Bounding{Q,FACE})(elem::Int) where {Q}
    return TopologyIterator(B, half4elem(B.t, elem))
end

function (B::Bounding{Q,SEGMENT})(edge::Int) where {Q}
    return TopologyIterator(B, half4edge(B.t, edge))
end

function Base.iterate(itr::TopologyIterator{<:Bounding{SEGMENT}})
    if isnothing(itr.e.elem)
        return nothing
    else
        return edge4pair(itr.relation.t, (itr.e.head, itr.e.half.head)), itr.e.next
    end
end
function Base.iterate(itr::TopologyIterator{<:Bounding{SEGMENT}}, state::HalfEdge)
    if state.next === itr.e.next
        return nothing
    else
        return edge4pair(itr.relation.t, (state.head, state.half.head)), state.next
    end
end

function Base.iterate(itr::TopologyIterator{<:Bounding{VERTEX,FACE}})
    if isnothing(itr.e.elem)
        return nothing
    else
        return itr.e.head, itr.e
    end
end

function Base.iterate(itr::TopologyIterator{<:Bounding{VERTEX,FACE}}, state::HalfEdge)
    n = state.next
    if n === itr.e
        return nothing
    else
        return n.head, n
    end
end

Base.IteratorSize(::Type{<:TopologyIterator{<:Bounding{VERTEX,SEGMENT}}}) = Base.HasLength()
Base.length(::TopologyIterator{<:Bounding{VERTEX,SEGMENT}}) = 2

Base.iterate(itr::TopologyIterator{<:Bounding{VERTEX,SEGMENT}}) = itr.e.head, itr.e.half
function Base.iterate(itr::TopologyIterator{<:Bounding{VERTEX,SEGMENT}}, state::HalfEdge)
    if state === itr.e
        return nothing
    else
        return state.head, state.half
    end
end

struct Shared{Q<:ParametricDimension,P<:ParametricDimension,T<:HalfEdgeTopology} <: AbstractRelation{Q,P}
    t::T

    function Shared{Q,P}(t::T) where {Q<:ParametricDimension,P<:ParametricDimension,T<:HalfEdgeTopology}
        _dim(Q) ≤ paramdim(t) || throw(DomainError((Q,T), "topology with rank $(ParametricDimension{paramdim(t)}) has no $Q. `Q` must be ≤ 2"))
        _dim(Q) > _dim(P) || throw(DomainError((Q,P), "cannot calculate shared $Q's of $P"))

        return new{Q,P,T}(t)
    end
end

function (S::Shared{Q,FACE})(elem::Int) where {Q}
    return TopologyIterator(S, half4elem(S.t, elem))
end

function (S::Shared{Q,SEGMENT})(edge::Int) where {Q}
    return TopologyIterator(S, half4edge(S.t, edge))
end

function (S::Shared{Q,VERTEX})(vert::Int) where {Q}
    return TopologyIterator(S, half4vert(S.t, vert))
end

function Base.iterate(itr::TopologyIterator{<:Shared{SEGMENT,VERTEX}})
    e = itr.e
    si = edge4pair(itr.relation.t, (e.head, e.half.head))
    return si, (e, :prev)
end

function Base.iterate(itr::TopologyIterator{<:Shared{SEGMENT,VERTEX}}, (e, dir)::Tuple{HalfEdge,Symbol})
    if dir === :prev
        if !isnothing(e.elem)
            h = e.prev.half
            if h.elem != itr.e.elem
                si = edge4pair(itr.relation.t, (h.head, h.half.head))
                return si, (h, :prev)
            else
                return nothing
            end
        elseif !isnothing(itr.e.half.elem)
            h = itr.e.half.next
            si = edge4pair(itr.relation.t, (h.head, h.half.head))
            return si, (h, :next)
        else
            return nothing
        end
    else # dir === :next
        h = e.half
        if !isnothing(h.elem)
            n = h.next
            si = edge4pair(itr.relation.t, (n.head, n.half.head))
            return si, (n, :next)
        else
            return nothing
        end
    end
end

function Base.iterate(itr::TopologyIterator{<:Shared{FACE,VERTEX}})
    if isnothing(itr.e.elem)
        h = itr.e.half
        return h.elem, (h, :next)
    else
        h = itr.e
        return h.elem, (h, :prev)
    end
end

function Base.iterate(itr::TopologyIterator{<:Shared{FACE,VERTEX}}, (e, dir)::Tuple{HalfEdge,Symbol})
    if dir === :prev
        h = e.prev.half
        if !isnothing(h.elem)
            if h.elem != itr.e.elem
                return h.elem, (h, :prev)
            else
                return nothing
            end
        elseif !isnothing(itr.e.half.elem)
            return itr.e.half.elem, (itr.e.half, :next)
        else
            return nothing
        end
    else # dir === :next
        h = e.next.half
        if !isnothing(h.elem)
            return h.elem, (h, :next)
        else
            return nothing
        end
    end
end

function Base.iterate(itr::TopologyIterator{<:Shared{FACE,SEGMENT}})
    e = itr.e
    if isnothing(e.elem)
        return e.half.elem, e.half
    else
        return e.elem, e
    end
end

function Base.iterate(itr::TopologyIterator{<:Shared{FACE,SEGMENT}}, e::HalfEdge)
    if e === itr.e.half || isnothing(e.half.elem)
        return nothing
    else
        return e.half.elem, e.half
    end
end

struct Adjacent{Q<:ParametricDimension,P<:ParametricDimension,T<:HalfEdgeTopology} <: AbstractRelation{Q,P}
    t::T

    function Adjacent{Q,P}(t::T) where {Q<:ParametricDimension,P<:ParametricDimension,T<:HalfEdgeTopology}
        _dim(Q) ≤ paramdim(t) || throw(DomainError((Q,T), "topology with rank $(ParametricDimension{paramdim(t)}) has no $Q. `Q` and `P` must be ≤ 2"))
        _dim(Q) == _dim(P) || throw(DomainError((Q,P), "cannot calculate adjacent $Q's of $P"))

        return new{Q,P,T}(t)
    end
end

# Convenience constructor to avoid repeating `Q`
Adjacent{Q}(t::HalfEdgeTopology) where Q = Adjacent{Q,Q}(t)
Adjacent{Q}(msh::HalfEdgeMesh) where Q = Adjacent{Q,Q}(topology(msh))

function (A::Adjacent{VERTEX})(vert::Int)
    return TopologyIterator(A, half4vert(A.t, vert))
end

function (A::Adjacent{FACE})(elem::Int)
    return TopologyIterator(A, half4elem(A.t, elem))
end

function Base.iterate(itr::TopologyIterator{<:Adjacent{VERTEX}})
    e = itr.e
    return e.half.head, (e, :prev)
end

function Base.iterate(itr::TopologyIterator{<:Adjacent{VERTEX}}, (e, dir)::Tuple{HalfEdge,Symbol})
    if dir === :prev
        if !isnothing(e.elem)
            h = e.prev.half
            if h.elem != itr.e.elem
                return h.half.head, (h, :prev)
            else
                return nothing
            end
        elseif !isnothing(itr.e.half.elem)
            h = itr.e.half.next
            return h.half.head, (h, :next)
        else
            return nothing
        end
    else # dir === :next
        h = e.half
        if !isnothing(h.elem)
            n = h.next
            return n.half.head, (n, :next)
        else
            return nothing
        end
    end
end

function Base.iterate(itr::TopologyIterator{<:Adjacent{FACE}})
    if isnothing(itr.e.elem)
        return nothing
    else
        n = itr.e
        while isnothing(n.half.elem)
            n = n.next
        end
        return n.half.elem, n.next
    end
end
function Base.iterate(itr::TopologyIterator{<:Adjacent{FACE}}, state::HalfEdge)
    while isnothing(state.half.elem) && state !== itr.e
        state = state.next
    end
    if state === itr.e
        return nothing
    else
        return state.half.elem, state.next
    end
end

