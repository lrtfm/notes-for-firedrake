import firedrake as fd
from firedrake.petsc import PETSc
import numpy as np
import json

__all__ = ('get_mesh_size', 'printf', 'sync_printf', 'NumpyEncoder', 'plot_errors')


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


def plot_errors(hs, errors, expect_order=2, filename=None):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import ticker

    hs = np.array(hs)
    errors = np.array(errors)
    mpl.rcParams['font.size'] = 14
    # 
    # default colors:
    # mpl.rcParams['axes.prop_cycle']
    # colors = ('tab:blue', 'tab:orange', 'tab:green', 'tab:red', 
    #           'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 
    #           'tab:olive', 'tab:cyan')

    fig, ax = plt.subplots(figsize=[5, 4])
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
