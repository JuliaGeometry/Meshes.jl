# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    issimple(polygon)

Tells whether or not the `polygon` is simple.
See [https://en.wikipedia.org/wiki/Simple_polygon]
(https://en.wikipedia.org/wiki/Simple_polygon).
"""
issimple(p::Polygon) = issimple(typeof(p))

issimple(p::PolyArea) = !hasholes(p) && issimple(first(rings(p)))

issimple(::Type{<:Ngon}) = true

"""
   issimple(chain)

Tells whether or not the `chain` is simple.

A chain is simple when all its segments only
intersect at end points.
"""
function issimple(c::Chain)
  λ(I) = !(type(I) == CornerTouching || type(I) == NotIntersecting)
  ss = collect(segments(c))
  result = Threads.Atomic{Bool}(true)
  Threads.@threads for i in 1:length(ss)
    for j in (i + 1):length(ss)
      if intersection(λ, ss[i], ss[j])
        result[] || break
        result[] = false
      end
    end
  end
  result[]
end
