# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct InDomain{CRS} <: CoordinateTransform end

InDomain(CRS) = InDomain{CRS}()

InDomain(code::Type{<:EPSG}) = InDomain(CoordRefSystems.get(code))

InDomain(code::Type{<:ESRI}) = InDomain(CoordRefSystems.get(code))

parameters(::InDomain{CRS}) where {CRS} = (; CRS)

function preprocess(::InDomain{CRS}, d::Domain) where {CRS}
  findall(d) do g
    all(pointify(g)) do p
      indomain(CRS, coords(p))
    end
  end
end

function preprocess(::InDomain{CRS}, d::Mesh) where {CRS}
  findall(d) do g
    all(eachvertex(g)) do p
      indomain(CRS, coords(p))
    end
  end
end

function preprocess(::InDomain{CRS}, d::GeometrySet{<:Any,<:Any,<:Union{Polytope,MultiPolytope}}) where {CRS}
  findall(d) do g
    all(eachvertex(g)) do p
      indomain(CRS, coords(p))
    end
  end
end

preprocess(::InDomain{CRS}, d::PointSet) where {CRS} = findall(p -> indomain(CRS, coords(p)), d)

apply(t::InDomain, d::Domain) = view(d, preprocess(t, d)), nothing
