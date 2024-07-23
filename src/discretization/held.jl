# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    HeldTriangulation([rng]; shuffle=true)

Fast Industrial-Strength Triangulation (FIST) of polygons.

This triangulation method is the method behind the famous Mapbox's
Earcut library. It is based on a ear clipping algorithm adapted
for complex n-gons with holes. It has O(n²) time complexity where
n is the number of vertices. In practice it is very efficient due
to heuristics implemented in the algorithm.

The option `shuffle` is used to shuffle the order in which ears
are clipped. It improves the quality of the triangles, which can
be very sliver otherwise. Optionally, specify the random number
generator `rng`.

## References

* Held, M. 1998. [FIST: Fast Industrial-Strength Triangulation of Polygons]
  (https://link.springer.com/article/10.1007/s00453-001-0028-4)
* Eder et al. 2018. [Parallelized ear clipping for the triangulation and
  constrained Delaunay triangulation of polygons]
  (https://www.sciencedirect.com/science/article/pii/S092577211830004X)
"""
struct HeldTriangulation{RNG<:AbstractRNG} <: BoundaryDiscretizationMethod
  rng::RNG
  shuffle::Bool
end

HeldTriangulation(rng=Random.default_rng(); shuffle=true) = HeldTriangulation(rng, shuffle)

function _discretizewithin(ring::Ring, ::Val{2}, method::HeldTriangulation)
  # helper function to shuffle ears
  earshuffle!(𝒬) = method.shuffle && shuffle!(method.rng, 𝒬)

  # input ring
  O = orientation(ring)
  ℛ = O == CCW ? ring : reverse(ring)

  # standardize coordinates
  𝒫 = ℛ |> StdCoords()

  # points of resulting mesh
  points = collect(vertices(ℛ))

  # standardized points for algorithm
  stdpts = collect(vertices(𝒫))

  # keep track of global indices
  I = CircularVector(1:nvertices(𝒫))

  # perform ear clipping
  𝒬 = earsccw(𝒫)
  earshuffle!(𝒬)
  n = nvertices(𝒫)
  𝒯 = Connectivity{Triangle,3}[]
  clipped = false
  while n > 3
    if !isempty(𝒬) # clip an ear
      # 0. select candidate ear
      i = pop!(𝒬)
      𝒬[𝒬 .> i] .-= 1
      # 1. push a new triangle to 𝒯
      push!(𝒯, connect((I[i - 1], I[i], I[i + 1])))
      # 2. remove the vertex from 𝒫
      I = I[setdiff(1:n, mod1(i, n))]
      𝒫 = Ring(stdpts[I])
      n = nvertices(𝒫)
      # 3. update 𝒬 near clipped ear
      for j in (i - 1, i)
        if isearccw(𝒫, j)
          𝒬 = 𝒬 ∪ [mod1(j, n)]
        else
          setdiff!(𝒬, [mod1(j, n)])
        end
      end
      clipped = true
    elseif clipped # recompute all ears
      𝒬 = earsccw(𝒫)
      earshuffle!(𝒬)
      clipped = false
    else # recovery process
      # check if consecutive edges vᵢ-1 -- vᵢ and vᵢ+1 -- vᵢ+2
      # intersect and fix the issue by clipping ear (vᵢ, vᵢ+1, vᵢ+2)
      v = vertices(𝒫)
      for i in 1:n
        s1 = Segment(v[i - 1], v[i])
        s2 = Segment(v[i + 1], v[i + 2])
        λ(I) = type(I) == Crossing
        if intersection(λ, s1, s2)
          # 1. push a new triangle to 𝒯
          push!(𝒯, connect((I[i], I[i + 1], I[i + 2])))
          # 2. remove the vertex from 𝒫
          I = I[setdiff(1:n, mod1(i + 1, n))]
          𝒫 = Ring(stdpts[I])
          n = nvertices(𝒫)
          clipped = true
          break
        end
      end

      # consecutive vertices vᵢ-1,  vᵢ, vᵢ+1 form a valid ear
      # if vᵢ-1 lies on the edge vᵢ+1 -- vᵢ+2
      v = vertices(𝒫)
      for i in 1:n
        if v[i - 1] ∈ Segment(v[i + 1], v[i + 2])
          # 1. push a new triangle to 𝒯
          push!(𝒯, connect((I[i - 1], I[i], I[i + 1])))
          # 2. remove the vertex from 𝒫
          I = I[setdiff(1:n, mod1(i, n))]
          𝒫 = Ring(stdpts[I])
          n = nvertices(𝒫)
          clipped = true
          break
        end
      end

      # enter in "desperate" mode and clip ears at random
      if !clipped
        # attempt to clip a convex vertex
        isconvex(i) = vexity(v, i) == :CONVEX
        j = findfirst(isconvex, 1:n)
        i = isnothing(j) ? rand(method.rng, 1:n) : j
        # 1. push a new triangle to 𝒯
        push!(𝒯, connect((I[i - 1], I[i], I[i + 1])))
        # 2. remove the vertex from 𝒫
        I = I[setdiff(1:n, mod1(i, n))]
        𝒫 = Ring(stdpts[I])
        n = nvertices(𝒫)
        clipped = true
      end
    end
  end

  # remaining polygonal area is the last triangle
  push!(𝒯, connect((I[1], I[2], I[3])))

  SimpleMesh(points, 𝒯)
end

# return index of all ears of 𝒫 assuming that 𝒫 is
# has counter-clockwise orientation
earsccw(𝒫) = filter(i -> isearccw(𝒫, i), 1:nvertices(𝒫))

# tells whether or not vertex i is an ear of 𝒫
# assuming that 𝒫 has counter-clockwise orientation
function isearccw(𝒫::Ring, i)
  v = vertices(𝒫)

  # CE1.1: classify angle as convex vs. reflex
  isconvex = vexity(v, i) == :CONVEX

  # CE1.2: check if segment vᵢ-₁ -- vᵢ+₁ intersects 𝒫
  λ(I) = !(type(I) == CornerTouching || type(I) == NotIntersecting)
  sᵢ = Segment(v[i - 1], v[i + 1])
  hasintersect = false
  for j in 1:nvertices(𝒫)
    sⱼ = Segment(v[j], v[j + 1])
    if intersection(λ, sᵢ, sⱼ)
      hasintersect = true
      break
    end
  end

  # CE1.3: check if vᵢ-1 ∈ C(vᵢ, vᵢ+1, vᵢ+2) and vᵢ+1 ∈ C(vᵢ-2, vᵢ-1, vᵢ)
  incones = incone(v, i - 1, i + 1) && incone(v, i + 1, i - 1)

  isconvex && !hasintersect && incones
end

# helper function to compute the vexity of vertex i
function vexity(v, i)
  α = ∠(v[i - 1], v[i], v[i + 1]) # oriented angle
  θ = α > 0 ? oftype(α, 2π) - α : -α # inner angle
  θ < π ? :CONVEX : :REFLEX
end

# helper function to check if vertex j is inside cone i
function incone(v, j, i)
  s1 = sideof(v[j], Line(v[i], v[i - 1]))
  s2 = sideof(v[j], Line(v[i], v[i + 1]))
  if vexity(v, i) == :CONVEX
    s1 != LEFT && s2 != RIGHT
  else
    s1 != LEFT || s2 != RIGHT
  end
end
