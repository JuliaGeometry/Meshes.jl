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
    connec = connect.([(4,5,6),(1,2,3),(3,4,6),(1,3,6)], Triangle)
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
                       (20,21,22), (20,22,23), (16,17,18),
                       (16,18,19), (19,20,23), (19,23,30),
                       (13,14,15), (15,16,19), (8,9,10),
                       (4,5,6), (6,7,8), (6,8,10),
                       (1,2,3), (3,4,6), (3,6,10),
                       (3,10,11), (11,12,13), (13,15,19),
                       (13,19,30), (30,1,3), (30,3,11),
                       (11,13,30)], Triangle)
    target = UnstructuredMesh(points, connec)
    𝒫 = PolyArea([points; points[1]])
    mesh = discretize(𝒫, FIST())
    @test mesh == target
    @test is_fully_connected(mesh)
    @test has_same_vertices(𝒫, mesh)

    poly = readpoly(joinpath(datadir, "poly2.line"))
    points = P2.(coordinates.(vertices(first(chains(poly)))))
    connec = connect.([(28, 29, 30), (24, 25, 26), (21, 22, 23),
                       (23, 24, 26), (18, 19, 20), (15, 16, 17),
                       (11, 12, 13), (8, 9, 10), (10, 11, 13),
                       (3, 4, 5), (21, 23, 26), (21, 26, 27),
                       (21, 27, 28), (14, 15, 17), (8, 10, 13),
                       (2, 3, 5), (2, 5, 6), (6, 7, 8),
                       (6, 8, 13), (6, 13, 14), (20, 21, 28),
                       (20, 28, 30), (2, 6, 14), (2, 14, 17),
                       (18, 20, 30), (1, 2, 17), (17, 18, 30),
                       (1, 17, 30)], Triangle)
    target = UnstructuredMesh(points, connec)
    𝒫 = PolyArea([points; points[1]])
    mesh = discretize(𝒫, FIST())
    @test mesh == target
    @test is_fully_connected(mesh)
    @test has_same_vertices(𝒫, mesh)

    poly = readpoly(joinpath(datadir, "poly3.line"))
    points = P2.(coordinates.(vertices(first(chains(poly)))))
    connec = connect.([(28, 29, 30), (24, 25, 26), (26, 27, 28),
                       (21, 22, 23), (23, 24, 26), (23, 26, 28),
                       (15, 16, 17), (12, 13, 14), (14, 15, 17),
                       (9, 10, 11), (9, 11, 12), (9, 12, 14),
                       (5, 6, 7), (21, 23, 28), (21, 28, 30),
                       (8, 9, 14), (8, 14, 17), (4, 5, 7),
                       (7, 8, 17), (20, 21, 30), (3, 4, 7),
                       (3, 7, 17), (3, 17, 18), (19, 20, 30),
                       (19, 30, 1), (2, 3, 18), (18, 19, 1),
                       (1, 2, 18)], Triangle)
    target = UnstructuredMesh(points, connec)
    𝒫 = PolyArea([points; points[1]])
    mesh = discretize(𝒫, FIST())
    @test mesh == target
    @test is_fully_connected(mesh)
    @test has_same_vertices(𝒫, mesh)

    poly = readpoly(joinpath(datadir, "poly4.line"))
    points = P2.(coordinates.(vertices(first(chains(poly)))))
    connec = connect.([(24, 25, 26), (20, 21, 22), (17, 18, 19),
                       (13, 14, 15), (15, 16, 17), (15, 17, 19),
                       (6, 7, 8), (6, 8, 9), (6, 9, 10),
                       (2, 3, 4), (4, 5, 6), (4, 6, 10),
                       (23, 24, 26), (23, 26, 27), (23, 27, 28),
                       (13, 15, 19), (13, 19, 20), (2, 4, 10),
                       (22, 23, 28), (22, 28, 29), (22, 29, 30),
                       (12, 13, 20), (20, 22, 30), (1, 2, 10),
                       (1, 10, 11), (11, 12, 20), (11, 20, 30),
                       (1, 11, 30)], Triangle)
    target = UnstructuredMesh(points, connec)
    𝒫 = PolyArea([points; points[1]])
    mesh = discretize(𝒫, FIST())
    @test mesh == target
    @test is_fully_connected(mesh)
    @test has_same_vertices(𝒫, mesh)

    poly = readpoly(joinpath(datadir, "poly5.line"))
    points = P2.(coordinates.(vertices(first(chains(poly)))))
    connec = connect.([(25, 26, 27), (19, 20, 21), (15, 16, 17),
                       (17, 18, 19), (17, 19, 21), (12, 13, 14),
                       (7, 8, 9), (4, 5, 6), (6, 7, 9),
                       (6, 9, 10), (1, 2, 3), (1, 3, 4),
                       (1, 4, 6), (24, 25, 27), (11, 12, 14),
                       (23, 24, 27), (23, 27, 28), (23, 28, 29),
                       (10, 11, 14), (10, 14, 15), (10, 15, 17),
                       (22, 23, 29), (22, 29, 30), (6, 10, 17),
                       (6, 17, 21), (21, 22, 30), (21, 30, 1),
                       (1, 6, 21)], Triangle)
    target = UnstructuredMesh(points, connec)
    𝒫 = PolyArea([points; points[1]])
    mesh = discretize(𝒫, FIST())
    @test mesh == target
    @test is_fully_connected(mesh)
    @test has_same_vertices(𝒫, mesh)

    poly = readpoly(joinpath(datadir, "smooth1.line"))
    points = P2.(coordinates.(vertices(first(chains(poly)))))
    connec = connect.([(119, 120, 1), (106, 107, 108), (106, 108, 109),
                       (103, 104, 105), (105, 106, 109), (105, 109, 110),
                       (105, 110, 111), (94, 95, 96), (94, 96, 97),
                       (94, 97, 98), (94, 98, 99), (94, 99, 100),
                       (91, 92, 93), (93, 94, 100), (86, 87, 88),
                       (83, 84, 85), (85, 86, 88), (85, 88, 89),
                       (85, 89, 90), (80, 81, 82), (82, 83, 85),
                       (82, 85, 90), (82, 90, 91), (91, 93, 100),
                       (91, 100, 101), (66, 67, 68), (63, 64, 65),
                       (65, 66, 68), (60, 61, 62), (62, 63, 65),
                       (62, 65, 68), (62, 68, 69), (62, 69, 70),
                       (62, 70, 71), (62, 71, 72), (54, 55, 56),
                       (51, 52, 53), (53, 54, 56), (53, 56, 57),
                       (46, 47, 48), (43, 44, 45), (45, 46, 48),
                       (30, 31, 32), (30, 32, 33), (30, 33, 34),
                       (30, 34, 35), (30, 35, 36), (27, 28, 29),
                       (29, 30, 36), (24, 25, 26), (26, 27, 29),
                       (26, 29, 36), (18, 19, 20), (15, 16, 17),
                       (17, 18, 20), (12, 13, 14), (14, 15, 17),
                       (14, 17, 20), (14, 20, 21), (9, 10, 11),
                       (11, 12, 14), (11, 14, 21), (11, 21, 22),
                       (6, 7, 8), (8, 9, 11), (8, 11, 22),
                       (8, 22, 23), (23, 24, 26), (23, 26, 36),
                       (3, 4, 5), (5, 6, 8), (8, 23, 36),
                       (8, 36, 37), (8, 37, 38), (8, 38, 39),
                       (8, 39, 40), (119, 1, 2), (2, 3, 5),
                       (5, 8, 40), (5, 40, 41), (5, 41, 42),
                       (118, 119, 2), (103, 105, 111), (103, 111, 112),
                       (103, 112, 113), (103, 113, 114), (103, 114, 115),
                       (79, 80, 82), (79, 82, 91), (60, 62, 72),
                       (60, 72, 73), (51, 53, 57), (51, 57, 58),
                       (43, 45, 48), (2, 5, 42), (42, 43, 48),
                       (117, 118, 2), (102, 103, 115), (78, 79, 91),
                       (59, 60, 73), (50, 51, 58), (58, 59, 73),
                       (117, 2, 42), (116, 117, 42), (116, 42, 48),
                       (101, 102, 115), (115, 116, 48), (77, 78, 91),
                       (91, 101, 115), (91, 115, 48), (91, 48, 49),
                       (49, 50, 58), (49, 58, 73), (49, 73, 74),
                       (76, 77, 91), (91, 49, 74), (91, 74, 75),
                       (75, 76, 91)], Triangle)
    target = UnstructuredMesh(points, connec)
    𝒫 = PolyArea([points; points[1]])
    mesh = discretize(𝒫, FIST())
    @test mesh == target
    @test is_fully_connected(mesh)
    @test has_same_vertices(𝒫, mesh)
  end
end
