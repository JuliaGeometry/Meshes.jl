# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    isless(A::Point, B::Point)

The product order of points `A` and `B`.

`isless(A, B)` if `aᵢ < bᵢ` for all coordinates `aᵢ` and `bᵢ`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
Base.isless(A::Point, B::Point) = A < B

"""
    <(A::Point, B::Point)

The product order of points `A` and `B`.

`A < B` if `aᵢ < bᵢ` for all coordinates `aᵢ` and `bᵢ`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
<(A::Point, B::Point) = all(x -> x > zero(x), B - A)
<(A::Point{🌐}, B::Point{🌐}) = _lat(A) < _lat(B)

"""
    >(A::Point, B::Point)

The product order of points `A` and `B`.

`A > B` if `aᵢ > bᵢ` for all coordinates `aᵢ` and `bᵢ`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
>(A::Point, B::Point) = all(x -> x > zero(x), A - B)
>(A::Point{🌐}, B::Point{🌐}) = _lat(A) > _lat(B)

"""
    ≤(A::Point, B::Point)

The product order of points `A` and `B`.

`A ≤ B` if `aᵢ ≤ bᵢ` for all coordinates `aᵢ` and `bᵢ`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
≤(A::Point, B::Point) = all(x -> x ≥ zero(x), B - A)
≤(A::Point{🌐}, B::Point{🌐}) = _lat(A) ≤ _lat(B)

"""
    ≥(A::Point, B::Point)

The product order of points `A` and `B`.

`A ≥ B` if `aᵢ ≥ bᵢ` for all coordinates `aᵢ` and `bᵢ`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
≥(A::Point, B::Point) = all(x -> x ≥ zero(x), A - B)
≥(A::Point{🌐}, B::Point{🌐}) = _lat(A) ≥ _lat(B)
