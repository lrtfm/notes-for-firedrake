# Utils for mesh
## Informations for the mesh

### Numbers of entities of mesh
```
mesh = Mesh('mymesh.msh')
mesh.num_cells()
mesh.num_vertices()
```

### Size of the mesh

```
mesh = RectangleMesh(10, 10, 1, 1)
h1 = CellSize(mesh)
h2 = CellDiameter(mesh)
```
