# Based on MIT licensed code from DataStructures.jl and DataStructures.jl#873 PR
# from @laborg
mutable struct ListNode{T}
    data::T
    prev::ListNode{T}
    next::ListNode{T}
    function ListNode{T}() where {T}
        node = new{T}()
        node.next = node.prev = node
        return node
    end
    function ListNode{T}(prev, data, next) where {T}
        node = new{T}(data, prev, next)
        prev.next = next.prev = node
        return node
    end
end

mutable struct MutableLinkedList{T}
    len::Int
    node::ListNode{T}
    function MutableLinkedList{T}() where {T}
        return new{T}(0, ListNode{T}())
    end
end

MutableLinkedList() = MutableLinkedList{Any}()
MutableLinkedList(rg::AbstractRange{T}) where {T} = append!(MutableLinkedList{T}(), rg)
MutableLinkedList(elts...) = MutableLinkedList{eltype(elts)}(elts...)
MutableLinkedList{T}(elts...) where {T} = append!(MutableLinkedList{T}(), elts)

Base.iterate(l::MutableLinkedList) = l.len == 0 ? nothing : (l.node.next.data, l.node.next.next)
Base.iterate(l::MutableLinkedList, n::ListNode) = n === l.node ? nothing : (n.data, n.next)

Base.isempty(l::MutableLinkedList) = l.len == 0
Base.length(l::MutableLinkedList) = l.len

# mend two nodes together
nconnect(a::ListNode{T}, b::ListNode{T}) where {T} = a.next, b.prev = b, a

# remove node from the list, the node itself is unchanged
function nremove(l::MutableLinkedList, node::ListNode)
    nconnect(node.prev, node.next)
    l.len -= 1
    return node
end

# create a new node, insert it after the provided `node` and return it
function ninsert(l::MutableLinkedList{T}, node::ListNode{T}, data) where {T}
    ins = ListNode{T}(node, data, node.next)
    l.len += 1
    return ins
end

function Base.append!(l::MutableLinkedList, collections...)
    node = l.node.prev
    for c in collections
        for e in c
            node = ninsert(l, node, e)
        end
    end
    return l
end

function Base.popfirst!(l::MutableLinkedList)
    isempty(l) && throw(ArgumentError("List must be non-empty"))
    node = l.node.next
    nremove(l, node)
    return node.data
end

function popprev!(l::MutableLinkedList, node)
    nremove(l, node.prev)
    return node
end

