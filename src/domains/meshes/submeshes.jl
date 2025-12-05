# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SubMesh{M,CRS}

A submesh of geometries in a given manifold `M` with point
coordinates specified in a coordinate reference system `CRS`.
"""
const SubMesh{M<:Manifold,C<:CRS} = SubDomain{M,C,<:Mesh{M,C}}

"""
    SubGrid{M,CRS}

A subgrid of geometries in a given manifold `M` with point
coordinates specified in a coordinate reference system `CRS`.
"""
const SubGrid{M<:Manifold,C<:CRS} = SubDomain{M,C,<:Grid{M,C}}
