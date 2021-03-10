@testset "Discretization" begin
  function is_fully_connected(mesh::UnstructuredMesh)
    inds = collect(1:length(vertices(mesh)))
    connected = collect(Iterators.flatten(getproperty.(mesh.connec, :list)))
    all(inds .== sort(unique(inds)))
  end
  function has_same_vertices(𝒫::PolyArea, mesh::UnstructuredMesh)
    Set(vcat(vertices.(chains(𝒫))...)) == Set(vertices(mesh))
  end

  @testset "FIST" begin
    𝒫 = Chain(P2[(0,0),(1,0),(1,1),(2,1),(2,2),(1,2),(0,0)])
    @test Meshes.ears(𝒫) == [2,4,5]

    𝒫 = Chain(P2[(0,0),(1,0),(1,1),(2,1),(1,2),(0,0)])
    @test Meshes.ears(𝒫) == [2,4]

    𝒫 = Chain(P2[(0,0),(1,0),(1,1),(1,2),(0,0)])
    @test Meshes.ears(𝒫) == [2,4]

    𝒫 = Chain(P2[(0,0),(1,1),(1,2),(0,0)])
    @test Meshes.ears(𝒫) == []

    𝒫 = Chain(P2[(0.443339268495331, 0.283757618605357),
                 (0.497822414616971, 0.398142813114205),
                 (0.770343126156527, 0.201815462842808),
                 (0.761236456732531, 0.330085709922366),
                 (0.985658085510286, 0.221530395507904),
                 (0.877899962498139, 0.325516131702896),
                 (0.561404274882782, 0.540334008885703),
                 (0.949459768187313, 0.396227653478068),
                 (0.594962560615951, 0.584927547374551),
                 (0.324208409133154, 0.607290684450708),
                 (0.424085089823892, 0.493532112641353),
                 (0.209843417261654, 0.590030658255966),
                 (0.27993878548962, 0.525162463476181),
                 (0.385557753911967, 0.322338556632868),
                 (0.443339268495331, 0.283757618605357)])
    @test Meshes.ears(𝒫) == [1,3,5,6,8,10,12,14]

    points = P2[(0,0),(1,0),(1,1),(2,1),(2,2),(1,2),(0,0)]
    connec = connect.([(4,5,6),(3,4,6),(3,6,1),(1,2,3)], Triangle)
    target = UnstructuredMesh(points[1:end-1], connec)
    𝒫 = PolyArea(points)
    mesh = discretize(𝒫, FIST())
    @test mesh == target
    @test is_fully_connected(mesh)
    @test has_same_vertices(𝒫, mesh)

    poly = readpoly(joinpath(datadir, "poly1.line"))
    points = P2.(coordinates.(vertices(first(chains(poly)))))
    connec = connect.([(26,27,28), (26,28,29), (23,24,25),
                       (23,25,26), (23,26,29), (23,29,30),
                       (21,22,23), (20,21,23), (16,17,18),
                       (16,18,19), (19,20,23), (19,23,30),
                       (13,14,15), (15,16,19), (11,12,13),
                       (13,15,19), (13,19,30), (8,9,10),
                       (7,8,10), (6,7,10), (4,5,6),
                       (4,6,10), (4,10,11), (11,13,30),
                       (11,30,1), (3,4,11), (3,11,1),
                       (1,2,3)], Triangle)
    target = UnstructuredMesh(points, connec)
    𝒫 = PolyArea([points; points[1]])
    mesh = discretize(𝒫, FIST())
    @test mesh == target
    @test is_fully_connected(mesh)
    @test has_same_vertices(𝒫, mesh)

    poly = readpoly(joinpath(datadir, "poly2.line"))
    points = P2.(coordinates.(vertices(first(chains(poly)))))
    connec = connect.([(28, 29, 30), (24, 25, 26), (23, 24, 26),
                       (22, 23, 26), (21, 22, 26), (21, 26, 27),
                       (21, 27, 28), (19, 20, 21), (19, 21, 28),
                       (19, 28, 30), (18, 19, 30), (15, 16, 17),
                       (17, 18, 30), (17, 30, 1), (17, 1, 2),
                       (14, 15, 17), (14, 17, 2), (11, 12, 13),
                       (13, 14, 2), (10, 11, 13), (9, 10, 13),
                       (8, 9, 13), (8, 13, 2), (7, 8, 2),
                       (6, 7, 2), (3, 4, 5), (5, 6, 2),
                       (2, 3, 5)], Triangle)
    target = UnstructuredMesh(points, connec)
    𝒫 = PolyArea([points; points[1]])
    mesh = discretize(𝒫, FIST())
    @test mesh == target
    @test is_fully_connected(mesh)
    @test has_same_vertices(𝒫, mesh)

    poly = readpoly(joinpath(datadir, "poly3.line"))
    points = P2.(coordinates.(vertices(first(chains(poly)))))
    connec = connect.([(28, 29, 30), (26, 27, 28), (24, 25, 26),
                       (23, 24, 26), (23, 26, 28), (21, 22, 23),
                       (21, 23, 28), (21, 28, 30), (19, 20, 21),
                       (19, 21, 30), (19, 30, 1), (15, 16, 17),
                       (12, 13, 14), (14, 15, 17), (12, 14, 17),
                       (10, 11, 12), (9, 10, 12), (8, 9, 12),
                       (5, 6, 7), (4, 5, 7), (3, 4, 7),
                       (7, 8, 12), (18, 19, 1), (3, 7, 12),
                       (12, 17, 18), (18, 1, 2), (18, 2, 3),
                       (3, 12, 18)], Triangle)
    target = UnstructuredMesh(points, connec)
    𝒫 = PolyArea([points; points[1]])
    mesh = discretize(𝒫, FIST())
    @test mesh == target
    @test is_fully_connected(mesh)
    @test has_same_vertices(𝒫, mesh)

    poly = readpoly(joinpath(datadir, "poly4.line"))
    points = P2.(coordinates.(vertices(first(chains(poly)))))
    connec = connect.([(24, 25, 26), (23, 24, 26), (23, 26, 27),
                       (23, 27, 28), (20, 21, 22), (22, 23, 28),
                       (22, 28, 29), (22, 29, 30), (17, 18, 19),
                       (16, 17, 19), (15, 16, 19), (13, 14, 15),
                       (13, 15, 19), (13, 19, 20), (20, 22, 30),
                       (11, 12, 13), (11, 13, 20), (11, 20, 30),
                       (11, 30, 1), (6, 7, 8), (6, 8, 9),
                       (6, 9, 10), (10, 11, 1), (10, 1, 2),
                       (4, 5, 6), (4, 6, 10), (4, 10, 2),
                       (2, 3, 4)], Triangle)
    target = UnstructuredMesh(points, connec)
    𝒫 = PolyArea([points; points[1]])
    mesh = discretize(𝒫, FIST())
    @test mesh == target
    @test is_fully_connected(mesh)
    @test has_same_vertices(𝒫, mesh)

    poly = readpoly(joinpath(datadir, "poly5.line"))
    points = P2.(coordinates.(vertices(first(chains(poly)))))
    connec = connect.([(25, 26, 27), (24, 25, 27), (23, 24, 27),
                       (23, 27, 28), (23, 28, 29), (19, 20, 21),
                       (18, 19, 21), (17, 18, 21), (15, 16, 17),
                       (12, 13, 14), (7, 8, 9), (6, 7, 9),
                       (6, 9, 10), (4, 5, 6), (2, 3, 4),
                       (2, 4, 6), (1, 2, 6), (30, 1, 6),
                       (22, 23, 29), (22, 29, 30), (11, 12, 14),
                       (10, 11, 14), (10, 14, 15), (10, 15, 17),
                       (6, 10, 17), (6, 17, 21), (6, 21, 22),
                       (6, 22, 30)], Triangle)
    target = UnstructuredMesh(points, connec)
    𝒫 = PolyArea([points; points[1]])
    mesh = discretize(𝒫, FIST())
    @test mesh == target
    @test is_fully_connected(mesh)
    @test has_same_vertices(𝒫, mesh)
  end
end
