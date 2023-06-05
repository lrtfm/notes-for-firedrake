from firedrake.mesh import MeshGeometry
from firedrake.petsc import PETSc
import numpy as np


def remove_pyop2_label(plex: PETSc.DMPlex):
    # Delete lables before create mesh with the plex, 
    # otherwise error occurs when run in parallel
    plex.removeLabel("pyop2_core")
    plex.removeLabel("pyop2_owned")
    plex.removeLabel("pyop2_ghost")
    return plex


def get_plex_with_update_coordinates(mesh: MeshGeometry):
    """
    Update the coordinates of the plex in mesh, and then return a clone without pyop2 label
    """
    mesh.topology.init()
    dm = mesh.topology_dm.clone()
    csec = dm.getCoordinateSection()
    coords_vec = dm.getCoordinatesLocal()

    s, e = dm.getDepthStratum(0)
    sec = mesh._vertex_numbering

    data = mesh.coordinates.dat.data_ro_with_halos
    dest = np.zeros_like(data)
    n = mesh.geometric_dimension()
    m = csec.getFieldComponents(0)
    assert m == n
    for i in range(s, e):
        # dof = sec.getDof(i)
        offset = sec.getOffset(i)
        # cdof = csec.getDof(i)
        coffset = csec.getOffset(i)

        dest[coffset//m] = data[offset, :]

    coords_vec.array_w[:] = dest.flatten()
    # dm.setCoordinatesLocal(coords_vec)
    remove_pyop2_label(dm)

    return dm


# With pr https://github.com/firedrakeproject/firedrake/pull/29330
# we can update the coordinates using section
# # Ref: https://github.com/firedrakeproject/firedrake/pull/2796/files
def get_plex_with_update_coordinates_new(mesh: MeshGeometry):
    tdim = mesh.topological_dimension()
    gdim = mesh.geometric_dimension()
    entity_dofs = np.zeros(tdim + 1, dtype=np.int32)
    entity_dofs[0] = gdim 
    coord_section = mesh.create_section(entity_dofs)
    plex = mesh.topology_dm.clone()
    coord_dm = plex.getCoordinateDM()
    coord_dm.setSection(coord_section)

    coords_local = coord_dm.createLocalVec()
    coords_local.array[:] = np.reshape(
        mesh.coordinates.dat.data_ro_with_halos, coords_local.array.shape
    )
    plex.setCoordinatesLocal(coords_local)
    remove_pyop2_label(plex)
    
    return plex


def test_get_plex_with_update_coordinates():
    import firedrake as fd
    import matplotlib.pyplot as plt

    def save_mesh(mesh, name):
        V = fd.FunctionSpace(mesh, 'CG', 1)
        f = fd.Function(V, name='f')
        fd.File(name).write(f)

    mesh_init = fd.RectangleMesh(5, 5, 1, 1)

    # move mesh
    mesh_init.coordinates.dat.data[:] += 1
    save_mesh(mesh_init, 'pvd/mesh_init.pvd')

    # recreate mesh from the plex
    plex = get_plex_with_update_coordinates(mesh_init)
    mesh = fd.Mesh(plex, distribution_parameters={"partition": False})
    save_mesh(mesh, 'pvd/mesh_with_update_plex.pvd')

    fig, ax = plt.subplots(1, 2, figsize=[9, 4], subplot_kw={})
    tp0 = fd.triplot(mesh_init, axes=ax[0])
    tp1 = fd.triplot(mesh, axes=ax[1])
    t0 = ax[0].set_title('Original mesh')
    t1 = ax[1].set_title('Mesh from Plex')


if __name__ == '__main__':
    test_get_plex_with_update_coordinates()