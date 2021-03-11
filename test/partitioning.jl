@testset "Partitioning" begin
  setify(lists) = Set(Set.(lists))

  d = CartesianGrid{T}(10,10)
  p = partition(d, RandomPartition(100))
  @test sprint(show, p) == "100 Partition{2,$T}"
  @test sprint(show, MIME"text/plain"(), p) == "100 Partition{2,$T}\n  └─1 View{10×10 CartesianGrid{2,$T}}\n  └─1 View{10×10 CartesianGrid{2,$T}}\n  └─1 View{10×10 CartesianGrid{2,$T}}\n  └─1 View{10×10 CartesianGrid{2,$T}}\n  └─1 View{10×10 CartesianGrid{2,$T}}\n  ⋮\n  └─1 View{10×10 CartesianGrid{2,$T}}\n  └─1 View{10×10 CartesianGrid{2,$T}}\n  └─1 View{10×10 CartesianGrid{2,$T}}\n  └─1 View{10×10 CartesianGrid{2,$T}}\n  └─1 View{10×10 CartesianGrid{2,$T}}"

  @testset "RandomPartition" begin
    Random.seed!(123)

    grid = CartesianGrid{T}(3,3)
    p = partition(grid, RandomPartition(3, false))
    @test setify(subsets(p)) == setify([[1,2,3], [4,5,6], [7,8,9]])
    p = partition(grid, RandomPartition(3))
    @test setify(subsets(p)) == setify([[8,6,9], [4,1,7], [2,3,5]])

    grid = CartesianGrid{T}(2,3)
    p = partition(grid, RandomPartition(3, false))
    @test setify(subsets(p)) == setify([[1,2], [3,4], [5,6]])
  end

  @testset "DirectionPartition" begin
    grid = CartesianGrid{T}(3,3)

    # basic checks on small regular grid data
    p = partition(grid, DirectionPartition(T.((1,0))))
    @test setify(subsets(p)) == setify([[1,2,3], [4,5,6], [7,8,9]])

    p = partition(grid, DirectionPartition(T.((0,1))))
    @test setify(subsets(p)) == setify([[1,4,7], [2,5,8], [3,6,9]])

    p = partition(grid, DirectionPartition(T.((1,1))))
    @test setify(subsets(p)) == setify([[1,5,9], [2,6], [3], [4,8], [7]])

    p = partition(grid, DirectionPartition(T.((1,-1))))
    @test setify(subsets(p)) == setify([[1], [2,4], [3,5,7], [6,8], [9]])

    # opposite directions produce same partition
    dir1 = (rand(T), rand(T)); dir2 = .-dir1
    p1 = partition(grid, DirectionPartition(dir1))
    p2 = partition(grid, DirectionPartition(dir2))
    @test setify(subsets(p1)) == setify(subsets(p2))

    # partition of arbitrarily large regular grid always
    # returns the "lines" and "columns" of the grid
    for n in [10,100,200]
      grid = CartesianGrid{T}(n,n)

      p = partition(grid, DirectionPartition(T.((1,0))))
      @test setify(subsets(p)) == setify([collect((i-1)*n+1:i*n) for i in 1:n])
      ns = [nelements(d) for d in p]
      @test all(ns .== n)

      p = partition(grid, DirectionPartition(T.((0,1))))
      @test setify(subsets(p)) == setify([collect(i:n:n*n) for i in 1:n])
      ns = [nelements(d) for d in p]
      @test all(ns .== n)
    end
  end

  @testset "FractionPartition" begin
    grid = CartesianGrid{T}(10,10)

    p = partition(grid, FractionPartition(T(0.5)))
    @test nelements(p[1]) == nelements(p[2]) == 50
    @test length(p) == 2

    p = partition(grid, FractionPartition(T(0.7)))
    @test nelements(p[1]) == 70
    @test nelements(p[2]) == 30

    p = partition(grid, FractionPartition(T(0.3)))
    @test nelements(p[1]) == 30
    @test nelements(p[2]) == 70
  end

  @testset "BlockPartition" begin
    grid = CartesianGrid{T}(10,10)

    p = partition(grid, BlockPartition(T(5),T(5)))
    @test length(p) == 4
    @test all(nelements.(p) .== 25)

    p = partition(grid, BlockPartition(T(5),T(2)))
    @test length(p) == 12
    @test Set(nelements.(p)) == Set([5,10])

    grid = CartesianGrid{T}(50, 50, 50)

    p = partition(grid, BlockPartition(T(1.), T(1.), T(1.), neighbors = false))
    @test length(p) == 125000
    @test Set(nelements.(p)) == Set(1)
    @test metadata(p) == Dict{Any,Any}()

    p = partition(grid, BlockPartition(T(5.), T(5.), T(5.), neighbors = true))
    @test length(p) == 1000
    @test Set(nelements.(p)) == Set(125)
    n = metadata(p)[:neighbors]
    @test length(n) == length(p)
    @test all(0 .< length.(n) .< 27) 
  end

  @testset "BisectPointPartition" begin
    grid = CartesianGrid((10,10), T.((-0.5,-0.5)), T.((1.0, 1.0)))

    p = partition(grid, BisectPointPartition(T.((0.0,1.0)), T.((5.0,5.1))))
    p1, p2 = p[1], p[2]
    @test nelements(p1) == 60
    @test nelements(p2) == 40

    # all points in p1 are below those in p2
    pts1 = [centroid(p1, i) for i in 1:nelements(p1)]
    pts2 = [centroid(p2, i) for i in 1:nelements(p2)]
    X1 = reduce(hcat, coordinates.(pts1))
    X2 = reduce(hcat, coordinates.(pts2))
    M1 = maximum(X1, dims=2)
    m2 = minimum(X2, dims=2)
    @test all(X1[2,j] < m2[2] for j in 1:size(X1,2))
    @test all(X2[2,j] > M1[2] for j in 1:size(X2,2))

    # flipping normal direction is equivalent to swapping subsets
    p₁ = partition(grid, BisectPointPartition(T.(( 1.0,0.0)), T.((5.1,5.0))))
    p₂ = partition(grid, BisectPointPartition(T.((-1.0,0.0)), T.((5.1,5.0))))
    @test nelements(p₁[1]) == nelements(p₂[2]) == 60
    @test nelements(p₁[2]) == nelements(p₂[1]) == 40
  end

  @testset "BisectFractionPartition" begin
    grid = CartesianGrid((10,10), T.((-0.5,-0.5)), T.((1.0,1.0)))

    p = partition(grid, BisectFractionPartition(T.((1.0,0.0)), T(0.2)))
    p1, p2 = p[1], p[2]
    @test nelements(p1) == 20
    @test nelements(p2) == 80

    # all points in p1 are to the left of p2
    pts1 = [centroid(p1, i) for i in 1:nelements(p1)]
    pts2 = [centroid(p2, i) for i in 1:nelements(p2)]
    X1 = reduce(hcat, coordinates.(pts1))
    X2 = reduce(hcat, coordinates.(pts2))
    M1 = maximum(X1, dims=2)
    m2 = minimum(X2, dims=2)
    @test all(X1[1,j] < m2[1] for j in 1:size(X1,2))
    @test all(X2[1,j] > M1[1] for j in 1:size(X2,2))

    # flipping normal direction is equivalent to swapping subsets
    p₁ = partition(grid, BisectFractionPartition(T.(( 1.0,0.0)), T(0.2)))
    p₂ = partition(grid, BisectFractionPartition(T.((-1.0,0.0)), T(0.8)))
    @test nelements(p₁[1]) == nelements(p₂[2]) == 20
    @test nelements(p₁[2]) == nelements(p₂[1]) == 80
  end

  @testset "BallPartition" begin
    pset = PointSet(T[
      0 1 1 0 0.2
      0 0 1 1 0.2
    ])

    # 3 balls with 1 point, and 1 ball with 2 points
    p = partition(pset, BallPartition(T(0.5)))
    n = nelements.(p)
    @test length(p) == 4
    @test count(i->i==1, n) == 3
    @test count(i->i==2, n) == 1
    @test setify(subsets(p)) == setify([[1,5],[2],[3],[4]])

    # 5 balls with 1 point each
    p = partition(pset, BallPartition(T(0.2)))
    @test length(p) == 5
    @test all(nelements.(p) .== 1)
    @test setify(subsets(p)) == setify([[1],[2],[3],[4],[5]])
  end

  @testset "PlanePartition" begin
    grid = CartesianGrid((3,3), T.((-0.5,-0.5)), T.((1.0,1.0)))
    p = partition(grid, PlanePartition(T.((0,1))))
    @test setify(subsets(p)) == setify([[1,2,3],[4,5,6],[7,8,9]])

    grid = CartesianGrid((4,4), T.((-0.5,-0.5)), T.((1.0,1.0)))
    p = partition(grid, PlanePartition(T.((0,1))))
    @test setify(subsets(p)) == setify([1:4,5:8,9:12,13:16])
  end

  @testset "PredicatePartition" begin
    grid = CartesianGrid((3,3), T.((-0.5,-0.5)), T.((1.0,1.0)))

    # partition even from odd locations
    pred(i,j) = iseven(i+j)
    p = partition(grid, PredicatePartition(pred))
    @test setify(subsets(p)) == setify([1:2:9,2:2:8])
  end

  @testset "SpatialPredicatePartition" begin
    g = CartesianGrid((10,10), T.((-0.5,-0.5)), T.((1.0,1.0)))

    # check if there are 100 partitions, each one having only 1 point
    sp = SpatialPredicatePartition((x,y) -> norm(x-y) < T(1))
    s = subsets(partition(g, sp))
    @test length(s) == 100
    nelms = [nelements(d) for d in partition(g, sp)]
    @test all(nelms .== 1)

    # defining a predicate to check if points x and y belong to the square [0.,5.]x[0.,5.]
    pred(x, y) = all(T[0,0] .<= x .<= T[5,5]) && all(T[0,0] .<= y .<= T[5,5])
    sp = SpatialPredicatePartition(pred)
    p = partition(g, sp)
    s = subsets(p)
    n = nelements.(p)

    # There will be 65 subsets:
    # 1 subset with 36 points (inside square [0.,5.]x[0.,5.])
    # 64 subsets with only 1 point inside each of them
    @test length(s) == 65
    @test maximum(length.(s)) == 36
    @test count(i->i==1, n) == 64
    @test count(i->i==36, n) == 1
  end

  @testset "ProductPartition" begin
    g = CartesianGrid((100,100), T.((-0.5,-0.5)), T.((1.0,1.0)))
    bm = BlockPartition(T(10), T(10))
    bn = BlockPartition(T(5), T(5))

    # Bm x Bn = Bn with m > n
    s1 = subsets(partition(g, bm*bn))
    s2 = subsets(partition(g, bn))
    @test setify(s1) == setify(s2)

    # pXp=p (for deterministic p)
    for p in [BlockPartition(T(10), T(10)),
              BisectFractionPartition(T.((0.1,0.1)))]
      s1 = subsets(partition(g, p*p))
      s2 = subsets(partition(g, p))
      @test setify(s1) == setify(s2)
    end
  end

  @testset "HierarchicalPartition" begin
    g = CartesianGrid((100,100), T.((-0.5,-0.5)), T.((1.0,1.0)))
    bm = BlockPartition(T(10), T(10))
    bn = BlockPartition(T(5), T(5))

    # Bn -> Bm = Bm with m > n
    s1 = subsets(partition(g, bm → bn))
    s2 = subsets(partition(g, bn))
    @test setify(s1) == setify(s2)
  end

  @testset "Mixed Tests" begin
    g = CartesianGrid((100,100), T.((-0.5,-0.5)), T.((1.0,1.0)))
    bm = BlockPartition(T(10), T(10))
    bn = BlockPartition(T(5), T(5))

    # Bm*Bn = Bm->Bn
    s1 = subsets(partition(g, bm * bn))
    s2 = subsets(partition(g, bm → bn))
    @test setify(s1) == setify(s2)
  end

  @testset "Utilities" begin
    d = CartesianGrid{T}(10,10)
    l, r = split(d, T(0.5))
    @test nelements(l) == 50
    @test nelements(r) == 50
    l, r = split(d, T(0.5), T.((1,0)))
    @test nelements(l) == 50
    @test nelements(r) == 50
    lpts = [centroid(l, i) for i in 1:nelements(l)]
    rpts = [centroid(r, i) for i in 1:nelements(r)]
    cl = mean(coordinates.(lpts))
    cr = mean(coordinates.(rpts))
    @test cl[1] < cr[1]
    @test cl[2] == cr[2]
    l, r = split(d, T(0.5), T.((0,1)))
    @test nelements(l) == 50
    @test nelements(r) == 50
    lpts = [centroid(l, i) for i in 1:nelements(l)]
    rpts = [centroid(r, i) for i in 1:nelements(r)]
    cl = mean(coordinates.(lpts))
    cr = mean(coordinates.(rpts))
    @test cl[1] == cr[1]
    @test cl[2] < cr[2]

    d = CartesianGrid{T}(10,10)
    s = slice(d, T(0.5):T(11), T(0.5):T(11))
    @test s == CartesianGrid(P2(1,1), P2(10,10), dims=(9,9))
  end
end
