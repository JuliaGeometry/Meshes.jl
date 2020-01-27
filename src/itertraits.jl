"""
    vertices(mesh)

An iterator for the vertices of the `mesh`.
"""
vertices(mesh::M) where M =
  [v for elm in elements(mesh) for v in vertices(mesh, elm)]

"""
    faces(mesh)

An iterator for the faces of the `mesh`.
"""
faces(mesh::M) where M =
  [f for elm in elements(mesh) for f in faces(mesh, elm)]
