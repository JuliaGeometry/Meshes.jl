# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, r::Ray{3,T}, p::Polyhedron) where {T}
    mesh = boundary(p)
    intersections = intersection.((r,), elements(mesh))
    filter!(I -> type(I)!=NoIntersection, intersections)
    isempty(intersections) && (return @IT NoIntersection nothing f)
    length(intersections) == 1 && (return only(intersections))
    return intersections
end
