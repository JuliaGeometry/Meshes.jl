# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ----------------------
# LEXICOGRAPHICAL ORDER
# ----------------------

"""
    <(A::Point, B::Point)

The lexicographical order of points `A` and `B` (`<`).

`A < B` if the tuples of coordinates satisfy `(a₁, a₂, ...) < (b₁, b₂, ...)`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
<(A::Point, B::Point) = CoordRefSystems.values(coords(A)) < CoordRefSystems.values(coords(B))

"""
    >(A::Point, B::Point)

The lexicographical order of points `A` and `B` (`>`).

`A > B` if the tuples of coordinates satisfy `(a₁, a₂, ...) > (b₁, b₂, ...)`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
>(A::Point, B::Point) = CoordRefSystems.values(coords(A)) > CoordRefSystems.values(coords(B))

"""
    ≤(A::Point, B::Point)

The lexicographical order of points `A` and `B` (`\\le`).

`A ≤ B` if the tuples of coordinates satisfy `(a₁, a₂, ...) ≤ (b₁, b₂, ...)`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
≤(A::Point, B::Point) = CoordRefSystems.values(coords(A)) ≤ CoordRefSystems.values(coords(B))

"""
    ≥(A::Point, B::Point)

The lexicographical order of points `A` and `B` (`\\ge`).

`A ≥ B` if the tuples of coordinates satisfy `(a₁, a₂, ...) ≥ (b₁, b₂, ...)`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
≥(A::Point, B::Point) = CoordRefSystems.values(coords(A)) ≥ CoordRefSystems.values(coords(B))

# --------------
# PRODUCT ORDER
# --------------

"""
    ≺(A::Point, B::Point)

The product order of points `A` and `B` (`\\prec`).

`A ≺ B` if `aᵢ < bᵢ` for all coordinates `aᵢ` and `bᵢ`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
≺(A::Point, B::Point) = all(x -> x > zero(x), B - A)
≺(A::Point{🌐}, B::Point{🌐}) = _lat(A) < _lat(B)

"""
    ≻(A::Point, B::Point)

The product order of points `A` and `B` (`\\succ`).

`A ≻ B` if `aᵢ > bᵢ` for all coordinates `aᵢ` and `bᵢ`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
≻(A::Point, B::Point) = all(x -> x > zero(x), A - B)
≻(A::Point{🌐}, B::Point{🌐}) = _lat(A) > _lat(B)

"""
    ⪯(A::Point, B::Point)

The product order of points `A` and `B` (`\\preceq`).

`A ⪯ B` if `aᵢ ≤ bᵢ` for all coordinates `aᵢ` and `bᵢ`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
⪯(A::Point, B::Point) = all(x -> x ≥ zero(x), B - A)
⪯(A::Point{🌐}, B::Point{🌐}) = _lat(A) ≤ _lat(B)

"""
    ⪰(A::Point, B::Point)

The product order of points `A` and `B` (`\\succeq`).

`A ⪰ B` if `aᵢ ≥ bᵢ` for all coordinates `aᵢ` and `bᵢ`.

See <https://en.wikipedia.org/wiki/Partially_ordered_set#Orders_on_the_Cartesian_product_of_partially_ordered_sets>
"""
⪰(A::Point, B::Point) = all(x -> x ≥ zero(x), A - B)
⪰(A::Point{🌐}, B::Point{🌐}) = _lat(A) ≥ _lat(B)
