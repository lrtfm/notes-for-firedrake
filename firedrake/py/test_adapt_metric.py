from firedrake import *
from firedrake.petsc import PETSc, OptionsManager
from firedrake.cython.dmcommon import to_petsc_local_numbering
import numpy as np
import matplotlib.pylab as plt


def to_petsc_local_numbering_for_local_vec(vec, V):
    section = V.dm.getLocalSection()
    out = vec.duplicate()
    varray = vec.array_r
    oarray = out.array
    dim = V.value_size
    idx = 0
    for p in range(*section.getChart()):
        dof = section.getDof(p)
        off = section.getOffset(p)
        # PETSc.Sys.syncPrint(f"dof = {dof}, offset = {off}")
        if dof > 0:
            off = section.getOffset(p)
            off *= dim
            for d in range(dof):
                for k in range(dim):
                    oarray[idx] = varray[off + dim * d + k]
                    idx += 1
    # PETSc.Sys.syncFlush()

    return out


def create_metric_from_indicator(indicator):
    """
    Create an metric for adpative mesh from and indicator

    indicator: Function defined on cells, < 1 coarsen, > 1 refine
    """

    mesh = indicator.ufl_domain()
    V = FunctionSpace(mesh, 'DG', 0)
    U = FunctionSpace(mesh, 'CG', 1)
    W = TensorFunctionSpace(mesh, 'CG', 1)

    metric = Function(W, name='metric')

    tv = Function(U, name='total_volume')
    cv = Function(V, name='cell_volume').interpolate(CellVolume(mesh))

    par_loop(('{[i] : 0 <= i < A.dofs}',  'A[i] = B[0]'),
            dx, {'A' : (tv, INC), 'B' : (cv, READ)})

    mesh = indicator.ufl_domain()

    degree = mesh.coordinates.function_space().ufl_element().degree()
    if degree != 1:
        raise Exception("We only support P1 mesh")

    coords = mesh.coordinates.dat.data_with_halos.real

    data = metric.dat.data_with_halos
    total_volume = tv.dat.data_with_halos.real
    cell_volume = cv.dat.data_with_halos.real
    Vc = mesh.coordinates.function_space()
    ind = indicator.dat.data_with_halos

    dim = mesh.geometric_dimension()

    if dim == 2:
        edge2vec = lambda e: [ e[0]**2, 2*e[0]*e[1], e[1]**2 ]
        vec2tensor = lambda g: [[g[0], g[1]], [g[1], g[2]]]
    elif dim == 3:
        edge2vec = lambda e:  [ e[0]**2, 2*e[0]*e[1], 2*e[0]*e[2], e[1]**2, 2*e[1]*e[2], e[2]**2 ]
        vec2tensor = lambda g: [[g[0], g[1], g[2]], [g[1], g[3], g[4]], [g[2], g[4], g[5]]]
    else:
        raise Exception

    pair = []
    for i in range(dim+1):
        for j in range(i):
            pair.append([i, j])

    for k, nodes in enumerate(Vc.cell_node_list):
        vertex = coords[nodes, :]
        edge = []
        for i, j in pair:
            edge.append(vertex[i, :] - vertex[j, :])

        mat = []
        for e in edge:
            mat.append(edge2vec(e))

        mat = np.array(mat)
        b = np.ones(len(edge))
        g = np.linalg.solve(mat, b)

        c = cell_volume[k] * ind[k]
        for n in nodes:
            data[n, :] += c/total_volume[n]*np.array(vec2tensor(g))

    return metric


def adapt(indicator):
    mesh = indicator.ufl_domain()
    metric = create_metric_from_indicator(indicator)
    v = PETSc.Vec().createWithArray(metric.dat._data,
                                    size=(metric.dat._data.size, None),
                                    bsize=metric.dat.cdim, comm=metric.comm)
    reordered = to_petsc_local_numbering_for_local_vec(v, metric.function_space())
    v.destroy()
    plex_new = mesh.topology_dm.adaptMetric(reordered, "Face Sets", "Cell Sets")
    mesh_new = Mesh(plex_new)

    return mesh_new


def test_adapt(dim=2, factor=2):
    if dim == 2:
        mesh = UnitSquareMesh(5, 5)
    elif dim == 3:
        mesh = UnitCubeMesh(5, 5, 5)
    else:
        raise Exception("")

    V = FunctionSpace(mesh, 'DG', 0)
    indicator = Function(V, name='indicator')
    indicator.dat.data[:] = factor
    mesh_adapt = adapt(indicator)

    mesh.topology_dm.viewFromOptions('-dm_view')
    mesh_adapt.topology_dm.viewFromOptions('-dm_view_new')
    return mesh, mesh_adapt


def get_avaiable_adaptors():
    from firedrake.petsc import get_external_packages
    eps = get_external_packages()
    adaptors = ["pragmatic", "mmg", "parmmg"]
    available_adaptors = []
    for apt in adaptors:
        if apt in eps:
            available_adaptors.append(apt)
    if len(available_adaptors) == 0:
        available_adaptors = None

    return available_adaptors


def test_adapt_with_option(dim=2, factor=0.5, adaptor=None):
    # adaptors: pragmatic, mmg, parmmg
    #   -dm_adaptor pragmatic
    #   -dm_adaptor mmg       # 2d or 3d (seq)
    #   -dm_adaptor parmmg    # 3d (parallel)

    adaptors = get_avaiable_adaptors()
    if adaptors is None:
        PETSc.Sys.Print("No adaptor found. You need to install pragmatic, mmg, or parmmg with petsc.")
        return None

    if adaptor is not None and adaptor not in adaptors:
        PETSc.Sys.Print(f"Adaptor {adaptor} not installed with petsc. Installed adaptors: {adaptors}.")
        return None

    if adaptor is None:
        adaptor = adaptors[0]

    parameters = {
        "dm_adaptor": adaptor,
        # "dm_plex_metric_target_complexity": 400,
        # "dm_view": None,
        # "dm_view_new": None,
    }

    om = OptionsManager(parameters, options_prefix="")
    with om.inserted_options():
        mesh, mesh_new = test_adapt(dim=dim, factor=factor)

    return mesh, mesh_new