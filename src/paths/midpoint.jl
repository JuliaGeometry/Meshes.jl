# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct MidPointPath <: Path end

sumofdists(point, points) = sum(evaluate(Euclidean(), point, p) for p in points)

function traverse(object, ::MidPointPath)
  path = [1]
  nelms = nelements(object)
  points = [centroid(object, i) for i in 1:nelms]
  for _ in 2:nelms
    pathpoints = points[path]
    otherinds = setdiff(1:nelms, path)
    next = argmax(i -> sumofdists(points[i], pathpoints), otherinds)
    push!(path, next)
  end
  path
end
