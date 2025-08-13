@testset "HalfEdgeTopology topological iteration" begin
#
#           3 ------------ 4
#         /  \            /
#        /    \    2     /
#       /      \        /
#      /        \      /
#     /    1     \    /
#    /            \  /
#   1 ------------ 2
#   This is a simple 2 triangle topology
    simple = HalfEdgeTopology(connect.([(1, 2, 3), (4, 3, 2)]))

#                  6
#                 /  \
#                /    \
#               /      \
#              /   3    \
#             /          \
#            /            \
#           1 ------------ 3
#          /  \           /  \
#         /    \    1    /    \
#        /      \       /      \
#       /        \     /        \
#      /    4     \   /    2     \
#     /            \ /            \
#    5 ------------ 2 ------------ 4
#   This topology has one triangle (1) which is fully surrounded (neighbors on every edge)
    triforce = HalfEdgeTopology(connect.([(1, 2, 3), (4, 3, 2), (1, 3, 6), (1, 2, 5)]))

#                   3
#                 / | \
#                /  |  \
#               /   |   \
#              /    |    \
#             /  1  |  2  \
#            /      |      \
#           1 ------2------ 4
#            \      |      /
#             \  4  |  3  /
#              \    |    /
#               \   |   /
#                \  |  /
#                 \ | /
#                   5
#   This topology has one vertex (2) which is fully surrounded (the connected edges are a cycle)
    diamond = HalfEdgeTopology(connect.([(1, 2, 3), (4, 3, 2), (2, 4, 5), (1, 2, 5)]))

    nonnothing_elem = x -> !isnothing(x.elem)
    @testset "Bounding" begin
        @test_throws DomainError Bounding{FACE,CELL}(simple)
        @test_throws DomainError Bounding{FACE,FACE}(simple)
        @test_throws DomainError Bounding{ParametricDimension{-1},VERTEX}(simple)

        @testset "Bounding{SEGMENT,FACE}" begin
            ∂ = Boundary{2,1}(simple)
            B = Bounding{SEGMENT,FACE}(simple)
            for i in 1:nfaces(simple,2)
                @test issetequal(B(i), ∂(i))
            end
            @test isnothing(iterate(TopologyIterator(B, half4pair(simple, (1,2)).half)))
            for e in filter(nonnothing_elem, simple.halfedges)
                @test issetequal(TopologyIterator(B, e), ∂(e.elem))
            end

            ∂ = Boundary{2,1}(triforce)
            B = Bounding{SEGMENT,FACE}(triforce)
            for i in 1:nfaces(triforce,2)
                @test issetequal(B(i), ∂(i))
            end
            for e in filter(nonnothing_elem, triforce.halfedges)
                @test issetequal(TopologyIterator(B, e), ∂(e.elem))
            end

            ∂ = Boundary{2,1}(diamond)
            B = Bounding{SEGMENT,FACE}(diamond)
            for i in 1:nfaces(diamond,2)
                @test issetequal(B(i), ∂(i))
            end
            for e in filter(nonnothing_elem, diamond.halfedges)
                @test issetequal(TopologyIterator(B, e), ∂(e.elem))
            end
        end

        @testset "Bounding{VERTEX,FACE}" begin
            ∂ = Boundary{2,0}(simple)
            B = Bounding{VERTEX,FACE}(simple)
            for i in 1:nfaces(simple,2)
                @test issetequal(B(i), ∂(i))
            end
            @test isnothing(iterate(TopologyIterator(B, half4pair(simple, (1,2)).half)))
            for e in filter(nonnothing_elem, simple.halfedges)
                @test issetequal(TopologyIterator(B, e), ∂(e.elem))
            end

            ∂ = Boundary{2,0}(triforce)
            B = Bounding{VERTEX,FACE}(triforce)
            for i in 1:nfaces(triforce,2)
                @test issetequal(B(i), ∂(i))
            end
            for e in filter(nonnothing_elem, triforce.halfedges)
                @test issetequal(TopologyIterator(B, e), ∂(e.elem))
            end

            ∂ = Boundary{2,0}(diamond)
            B = Bounding{VERTEX,FACE}(diamond)
            for i in 1:nfaces(diamond,2)
                @test issetequal(B(i), ∂(i))
            end
            for e in filter(nonnothing_elem, diamond.halfedges)
                @test issetequal(TopologyIterator(B, e), ∂(e.elem))
            end
        end

        @testset "Bounding{VERTEX,SEGMENT}" begin
            ∂ = Boundary{1,0}(simple)
            B = Bounding{VERTEX,SEGMENT}(simple)
            for i in 1:nfacets(simple)
                @test issetequal(B(i), ∂(i))
            end
            for e in simple.halfedges
                @test issetequal(TopologyIterator(B, e), ∂(edge4pair(simple, (e.head, e.half.head))))
            end

            ∂ = Boundary{1,0}(triforce)
            B = Bounding{VERTEX,SEGMENT}(triforce)
            for i in 1:nfacets(triforce)
                @test issetequal(B(i), ∂(i))
            end
            for e in triforce.halfedges
                @test issetequal(TopologyIterator(B, e), ∂(edge4pair(triforce, (e.head, e.half.head))))
            end

            ∂ = Boundary{1,0}(diamond)
            B = Bounding{VERTEX,SEGMENT}(diamond)
            for i in 1:nfacets(diamond)
                @test issetequal(B(i), ∂(i))
            end
            for e in diamond.halfedges
                @test issetequal(TopologyIterator(B, e), ∂(edge4pair(diamond, (e.head, e.half.head))))
            end
        end
    end

    @testset "Shared" begin
        @test_throws DomainError Bounding{CELL,FACE}(simple)
        @test_throws DomainError Bounding{FACE,FACE}(simple)

        @testset "Shared{SEGMENT,VERTEX}" begin
            𝒞 = Coboundary{0,1}(simple)
            S = Shared{SEGMENT,VERTEX}(simple)
            for i in 1:nvertices(simple)
                @test issetequal(S(i), 𝒞(i))
            end
            for e in simple.halfedges
                @test issetequal(TopologyIterator(S, e), 𝒞(e.head))
            end

            𝒞 = Coboundary{0,1}(triforce)
            S = Shared{SEGMENT,VERTEX}(triforce)
            for i in 1:nvertices(triforce)
                @test issetequal(S(i), 𝒞(i))
            end
            for e in triforce.halfedges
                @test issetequal(TopologyIterator(S, e), 𝒞(e.head))
            end

            𝒞 = Coboundary{0,1}(diamond)
            S = Shared{SEGMENT,VERTEX}(diamond)
            for i in 1:nvertices(diamond)
                @test issetequal(S(i), 𝒞(i))
            end
            for e in diamond.halfedges
                @test issetequal(TopologyIterator(S, e), 𝒞(e.head))
            end
        end

        @testset "Shared{FACE,VERTEX}" begin
            𝒞 = Coboundary{0,2}(simple)
            S = Shared{FACE,VERTEX}(simple)
            for i in 1:nvertices(simple)
                @test issetequal(S(i), 𝒞(i))
            end
            for e in simple.halfedges
                @test issetequal(TopologyIterator(S, e), 𝒞(e.head))
            end

            𝒞 = Coboundary{0,2}(triforce)
            S = Shared{FACE,VERTEX}(triforce)
            for i in 1:nvertices(triforce)
                @test issetequal(S(i), 𝒞(i))
            end
            for e in triforce.halfedges
                @test issetequal(TopologyIterator(S, e), 𝒞(e.head))
            end

            𝒞 = Coboundary{0,2}(diamond)
            S = Shared{FACE,VERTEX}(diamond)
            for i in 1:nvertices(diamond)
                @test issetequal(S(i), 𝒞(i))
            end
            for e in diamond.halfedges
                @test issetequal(TopologyIterator(S, e), 𝒞(e.head))
            end
        end

        @testset "Shared{FACE,SEGMENT}" begin
            𝒞 = Coboundary{1,2}(simple)
            S = Shared{FACE,SEGMENT}(simple)
            for i in 1:nfacets(simple)
                @test issetequal(S(i), 𝒞(i))
            end
            for (e,i) in zip(simple.halfedges, repeat(1:nfacets(simple), inner=2))
                @test issetequal(TopologyIterator(S, e), 𝒞(i))
            end

            𝒞 = Coboundary{1,2}(triforce)
            S = Shared{FACE,SEGMENT}(triforce)
            for i in 1:nfacets(triforce)
                @test issetequal(S(i), 𝒞(i))
            end
            for (e,i) in zip(triforce.halfedges, repeat(1:nfacets(triforce), inner=2))
                @test issetequal(TopologyIterator(S, e), 𝒞(i))
            end

            𝒞 = Coboundary{1,2}(diamond)
            S = Shared{FACE,SEGMENT}(diamond)
            for i in 1:nfacets(diamond)
                @test issetequal(S(i), 𝒞(i))
            end
            for (e,i) in zip(diamond.halfedges, repeat(1:nfacets(diamond), inner=2))
                @test issetequal(TopologyIterator(S, e), 𝒞(i))
            end
        end
    end

    @testset "Adjacent" begin
        @test_throws DomainError Adjacent{CELL,FACE}(simple)
        @test_throws DomainError Adjacent{VERTEX,SEGMENT}(simple)

        @testset "Adjacent{VERTEX,VERTEX}" begin
            𝒜 = Adjacency{0}(simple)
            A = Adjacent{VERTEX}(simple)
            for i in 1:nvertices(simple)
                @test issetequal(A(i), 𝒜(i))
            end
            for e in simple.halfedges
                @test issetequal(TopologyIterator(A, e), 𝒜(e.head))
            end

            𝒜 = Adjacency{0}(triforce)
            A = Adjacent{VERTEX}(triforce)
            for i in 1:nvertices(triforce)
                @test issetequal(A(i), 𝒜(i))
            end
            for e in triforce.halfedges
                @test issetequal(TopologyIterator(A, e), 𝒜(e.head))
            end

            𝒜 = Adjacency{0}(diamond)
            A = Adjacent{VERTEX}(diamond)
            for i in 1:nvertices(diamond)
                @test issetequal(A(i), 𝒜(i))
            end
            for e in diamond.halfedges
                @test issetequal(TopologyIterator(A, e), 𝒜(e.head))
            end
        end

        @testset "Adjacent{FACE,FACE}" begin
            𝒜 = Adjacency{2}(simple)
            A = Adjacent{FACE}(simple)
            for i in 1:nfaces(simple,2)
                @test issetequal(A(i), 𝒜(i))
            end
            @test isnothing(iterate(TopologyIterator(A, half4pair(simple, (1,2)).half)))
            for e in filter(nonnothing_elem, simple.halfedges)
                @test issetequal(TopologyIterator(A, e), 𝒜(e.elem))
            end

            𝒜 = Adjacency{2}(triforce)
            A = Adjacent{FACE}(triforce)
            for i in 1:nfaces(triforce,2)
                @test issetequal(A(i), 𝒜(i))
            end
            for e in filter(nonnothing_elem, triforce.halfedges)
                @test issetequal(TopologyIterator(A, e), 𝒜(e.elem))
            end

            𝒜 = Adjacency{2}(diamond)
            A = Adjacent{FACE}(diamond)
            for i in 1:nfaces(diamond,2)
                @test issetequal(A(i), 𝒜(i))
            end
            for e in filter(nonnothing_elem, diamond.halfedges)
                @test issetequal(TopologyIterator(A, e), 𝒜(e.elem))
            end
        end
    end
end

