# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ∠(A, B, C)

Angle ∠ABC between rays BA and BC.
See https://en.wikipedia.org/wiki/Angle.

## Example

```julia
∠(Point(0,1), Point(0,0), (1,0)) == π/2
```

## References

* Balbes, R. and Siegel, J. 1990. [A robust method for calculating
  the simplicity and orientation of planar polygons]
  (https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
function ∠(A::P, B::P, C::P) where {P<:Point{2}}
  BA = A - B
  BC = C - B

  cross = BA × BC
  inner = BA ⋅ BC

  R = norm(cross) / (norm(BA) * norm(BC))

  # Table 1 from Balbes, R. and Siegel, J. 1990.
  if cross ≥ 0 && inner ≥ 0
    asin(R)
  elseif cross ≥ 0 && inner < 0
    π - asin(R)
  elseif cross < 0 && inner ≥ 0
    -asin(R)
  elseif cross < 0 && inner < 0
    asin(R) - π
  end
end
