struct Pyramid{T} <: GeometryPrimitive{3,T}
    middle::Point{3,T}
    length::T
    width::T
end

function coordinates(p::Pyramid{T}) where {T}
    leftup = Point{3,T}(-p.width, p.width, 0) / 2
    leftdown = Point(-p.width, -p.width, 0) / 2
    tip = Point{3,T}(p.middle + Point{3,T}(0, 0, p.length))
    lu = Point{3,T}(p.middle + leftup)
    ld = Point{3,T}(p.middle + leftdown)
    ru = Point{3,T}(p.middle - leftdown)
    rd = Point{3,T}(p.middle - leftup)
    return Point{3,T}[tip, rd, ru, tip, ru, lu, tip, lu, ld, tip, ld, rd, rd, ru, lu, lu,
                      ld, rd]
end

function faces(::Pyramid)
    return (TriangleFace(triangle) for triangle in TupleView{3}(1:18))
end
