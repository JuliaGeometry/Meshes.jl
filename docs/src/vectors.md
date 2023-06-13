# Vectors

```@example vectors
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
```

```@example vectors
using Meshes, MeshViz # hide
import WGLMakie as Mke # hide
```

```@docs
Vec
```

```@example vectors
rand(Vec3, 100)
```