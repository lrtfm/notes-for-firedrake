import firedrake as fd
import firedrake.pyplot as fdplt
from firedrake.petsc import PETSc
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.collections import LineCollection, PolyCollection
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
import numpy as np
import json
import logging
import re

__all__ = ('triplot', 'plot_mesh_with_label', 'get_mesh_size', 'printf', 'sync_printf', 'NumpyEncoder', 'plot_errors')


_EMPTY_SUBDOMAIN_WARNING = re.compile(
    r"Subdomain \(\d+(?:,\s*\d+)*,?\) is empty\. "
    r"This is likely an error\. Did you choose the right label\?"
)


class _EmptySubdomainWarningFilter(logging.Filter):
    def filter(self, record):
        return _EMPTY_SUBDOMAIN_WARNING.fullmatch(record.getMessage()) is None


def triplot(mesh, axes=None, interior_kw={}, boundary_kw={}):
    """Call Firedrake's triplot without its empty-facet probe warnings.

    Since Firedrake 2026.10 the plotting functions refuse to run when the
    mesh is distributed (``SerialExecutionOnlyError``). In that case each
    process draws its own local part of the mesh, halo cells included, which
    is what the Firedrake function used to do in parallel.
    """
    logger = logging.getLogger('firedrake')
    log_filter = _EmptySubdomainWarningFilter()
    handlers = tuple(logger.handlers)
    for handler in handlers:
        handler.addFilter(log_filter)
    try:
        if mesh.comm.size > 1:
            return _triplot_local(mesh, axes, interior_kw, boundary_kw)
        return fdplt.triplot(mesh, axes=axes, interior_kw=interior_kw, boundary_kw=boundary_kw)
    finally:
        for handler in handlers:
            handler.removeFilter(log_filter)


def _triplot_local(mesh, axes=None, interior_kw={}, boundary_kw={}):
    """Plot the local part of a distributed mesh (owned and halo cells).

    Mirrors ``firedrake.pyplot.triplot`` as of Firedrake 2026.4: cells are
    drawn from the coordinate data with halos and the marked exterior facets
    are colored per marker. The ``subset`` lookup of each marker is
    collective, so all processes must call this together.
    """
    gdim = mesh.geometric_dimension
    tdim = mesh.topological_dimension
    if gdim == 2 and tdim == 2:
        BoundaryCollection, InteriorCollection = LineCollection, PolyCollection
    elif gdim == 3 and tdim == 2:
        BoundaryCollection, InteriorCollection = Line3DCollection, Poly3DCollection
    elif gdim == 3 and tdim == 3:
        BoundaryCollection, InteriorCollection = Poly3DCollection, Poly3DCollection
    else:
        raise NotImplementedError(f"Cannot plot a mesh with gdim={gdim}, tdim={tdim}")
    quad = mesh.ufl_cell().cellname == "quadrilateral"

    if mesh.extruded:
        raise NotImplementedError("Visualizing extruded meshes not implemented yet!")

    if axes is None:
        figure = plt.figure()
        axes = figure.add_subplot(111, projection='3d' if gdim == 3 else None)

    coordinates = mesh.coordinates
    element = coordinates.function_space().ufl_element()
    if element.degree() != 1:
        V = fd.VectorFunctionSpace(mesh, element.family(), 1)
        coordinates = fd.assemble(fd.interpolate(coordinates, V))

    coords = np.real(coordinates.dat.data_ro_with_halos)
    result = []
    interior_kw = dict(interior_kw)
    if tdim <= 2:
        cell_node_map = coordinates.cell_node_map().values_with_halo
        idx = (tuple(range(tdim + 1)) if not quad else (0, 1, 3, 2)) + (0,)
        vertices = coords[cell_node_map[:, idx]]
        interior_kw.setdefault("edgecolors", "k")
        interior_kw.setdefault("linewidths", 1.0)
        if gdim == 2:
            interior_kw.setdefault("facecolors", "none")
        interior_collection = InteriorCollection(vertices, **interior_kw)
        axes.add_collection(interior_collection)
        result.append(interior_collection)

    topology = coordinates.function_space().finat_element.cell.get_topology()

    def facet_data(typ):
        if typ == "interior":
            facets = mesh.interior_facets
            node_map = coordinates.interior_facet_node_map()
            node_map = node_map.values_with_halo[:, :node_map.arity//2]
            local_facet_ids = facets.local_facet_dat.data_ro_with_halos[:, :1].reshape(-1)
        else:
            facets = mesh.exterior_facets
            local_facet_ids = facets.local_facet_dat.data_ro_with_halos
            node_map = coordinates.exterior_facet_node_map().values_with_halo
        mask = np.zeros(node_map.shape, dtype=bool)
        for facet_index, local_facet_index in enumerate(local_facet_ids):
            mask[facet_index, topology[tdim - 1][local_facet_index]] = True
        faces = node_map[mask].reshape(-1, tdim)
        return facets, faces

    markers = mesh.exterior_facets.unique_markers
    color_key = "colors" if tdim <= 2 else "facecolors"
    boundary_kw = dict(boundary_kw)
    boundary_colors = boundary_kw.pop(color_key, None)
    if boundary_colors is None:
        cmap = mpl.colormaps["Dark2"]
        num_markers = len(markers)
        colors = cmap([k / num_markers for k in range(num_markers)])
    else:
        colors = mpl.colors.to_rgba_array(boundary_colors)
    if tdim == 3:
        boundary_kw.setdefault("edgecolors", "k")
        boundary_kw.setdefault("linewidths", 1.0)
    for marker, color in zip(markers, colors):
        vertices = []
        for typ in ["interior", "exterior"]:
            facets, faces = facet_data(typ)
            face_indices = facets.subset(int(marker)).indices
            vertices.append(coords[faces[face_indices, :]])
        vertices = np.concatenate(vertices)
        marker_collection = BoundaryCollection(vertices, **{color_key: color, "label": marker}, **boundary_kw)
        axes.add_collection(marker_collection)
        result.append(marker_collection)

    for collection in result:
        if isinstance(collection, Poly3DCollection):
            collection._facecolors2d = PolyCollection.get_facecolor(collection)
            collection._edgecolors2d = PolyCollection.get_edgecolor(collection)

    axes.autoscale_view()
    if isinstance(axes, Axes3D):
        for setter, idx in zip(["set_xlim", "set_ylim", "set_zlim"], range(coords.shape[1])):
            amin, amax = coords[:, idx].min(), coords[:, idx].max()
            extra = (amax - amin) / 20 or 0.5
            getattr(axes, setter)(amin - extra, amax + extra)
    return result


def plot_mesh_with_label(mesh, axes=None):
    if axes is None:
        fig, axes = plt.subplots(figsize=[4, 4])
    triplot(mesh, axes=axes, boundary_kw={'lw': 3})
    axes.set_aspect(aspect='equal')
    # ax.set_axis_off()
    axes.legend(loc='upper left', bbox_to_anchor=(1, 1))
    axes.set_xlabel('x')
    axes.set_ylabel('y')


def get_mesh_size(mesh):
    degree = mesh.coordinates.function_space().ufl_element().degree()
    if degree > 1:
        V_c = fd.VectorFunctionSpace(mesh, 'CG', 1)
        coords = fd.Function(V_c).interpolate(mesh.coordinates)
        mesh = fd.Mesh(coords)

    V = fd.FunctionSpace(mesh, 'DG', 0)
    h = fd.Function(V)
    h.interpolate(fd.CellDiameter(mesh))

    with h.dat.vec_ro as vec:
        _, h_max = vec.max()

    return h_max


# print for parallel
def printf(*args, **kwargs):
    PETSc.Sys.Print(*args, **kwargs)


def sync_printf(*args, **kwargs):
    # comm = kwargs.get('comm', PETSc.COMM_WORLD)
    # size, rank = comm.size, comm.rank
    # prefix = f'[{rank}/{size}]'
    # PETSc.Sys.syncPrint(prefix, *args, **kwargs)
    PETSc.Sys.syncPrint(*args, **kwargs)
    PETSc.Sys.syncFlush()


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def plot_errors(hs, errors, expect_order=2, filename=None, axes=None):

    hs = np.array(hs)
    errors = np.array(errors)
    mpl.rcParams['font.size'] = 14
    # 
    # default colors:
    # mpl.rcParams['axes.prop_cycle']
    # colors = ('tab:blue', 'tab:orange', 'tab:green', 'tab:red', 
    #           'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 
    #           'tab:olive', 'tab:cyan')

    if axes is None:
        fig, ax = plt.subplots(figsize=[5, 4])
    else:
        ax = axes
    ax.loglog(hs, errors, marker='v', ls='-', color='tab:blue', label=f'$L^2$')

    # plot reference line
    epo = expect_order
    ax.loglog(hs, 2*errors[0]*hs**epo/hs[0]**epo, ls='-.', color='tab:blue', label=f'$O(h^{epo})$')

    def my_log_formatter(s, pos):
        return '$1/' + str(int(1/s)) + '$'

    ax.xaxis.set_major_locator(ticker.LogLocator(2))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.xaxis.set_major_formatter(my_log_formatter)
    ax.xaxis.set_minor_formatter(my_log_formatter)

    ax.set_xlabel('$h$')
    ax.set_ylabel('$L^2$-Error')
    ax.legend(loc='upper left', fontsize=10)
    # ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=10)
    if filename:
        plt.savefig(filename, bbox_inches='tight')
