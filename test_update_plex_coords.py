from firedrake import *
import numpy as np


def remove_pyop2_label(plex):
    plex.removeLabel("pyop2_core")
    plex.removeLabel("pyop2_owned")
    plex.removeLabel("pyop2_ghost")
    return plex


def get_plex_with_update_coordinates(mesh):
    """
    Update the coordinates of the plex in mesh, and then return a clone without pyop2 label
    """
    mesh.topology.init()
    dm = mesh.topology_dm
    # cdm = dm.getCoordinateDM()
    # dim = dm.getCoordinateDim()
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
        dof = sec.getDof(i)
        offset = sec.getOffset(i)
        cdof = csec.getDof(i)
        coffset = csec.getOffset(i)

        dest[coffset//m] = data[offset, :]

    coords_vec.array_w[:] = dest.flatten()
    # dm.setCoordinatesLocal(coords_vec)
    dm = dm.clone()
    remove_pyop2_label(dm)

    return dm

def save_mesh(mesh, name):
    V = FunctionSpace(mesh, 'CG', 1)
    f = Function(V, name='f')
    File(name).write(f)

mesh_init = RectangleMesh(5, 5, 1, 1)
save_mesh(mesh_init, 'pvd/mesh_init.pvd')

plex2 = get_plex_with_update_coordinates(mesh_init)
plex2.view()
mesh2 = Mesh(plex2, distribution_parameters={"partition": False})
save_mesh(mesh2, 'pvd/mesh_with_update_plex.pvd')
plex2.view()
