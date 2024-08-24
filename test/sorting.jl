@testitem "DirectionSort" begin
  g = cartgrid(3, 3)
  s = sort(g, DirectionSort((T(1), T(1))))
  @test centroid.(s) ==
        cart.([
    (0.5, 0.5),
    (1.5, 0.5),
    (0.5, 1.5),
    (2.5, 0.5),
    (1.5, 1.5),
    (0.5, 2.5),
    (2.5, 1.5),
    (1.5, 2.5),
    (2.5, 2.5)
  ])
end
