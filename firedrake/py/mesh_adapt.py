from sqlite3 import paramstyle
from firedrake import *
from firedrake.petsc import PETSc, flatten_parameters
import numpy as np

def plot_with_cell_number(mesh):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    triplot(mesh, axes=ax)

    plex = mesh.topology_dm

    cs, ce = plex.getHeightStratum(0)
    ps, pe = plex.getHeightStratum(2)
    csec = plex.getCoordinateSection()
    coords = plex.getCoordinatesLocal()
    for i in range(cs, ce):
        cl, _ = plex.getTransitiveClosure(i)
        _index = np.where((cl >= ps) & (cl < pe))
        _ps = cl[_index]
        center = np.array([0.0, 0.0])
        for _p in _ps:
            offset = csec.getOffset(_p)
            ndof = csec.getDof(_p)
            if ndof > 0 and offset > 0:
                center += coords[offset:offset+2]
        center /= 3

        ax.text(center[0], center[1], str(i), ha='center', va='center')

class YAOptionManager(object):
    commandline_options = frozenset(PETSc.Options().getAll())
    
    def __init__(self, parameters, options_prefix=None, cmd_high_priority=True):
        """Yet Another OptionsManager
        
        Args:
            parameters: a dict of parameters
            options_prefix: options prefix, default: None, i.e ''
            cmd_high_priority: commandline args with highest priority if true, else with lowest priority.
        """

        if options_prefix == None:
            options_prefix = '' 

        if len(options_prefix) and not options_prefix.endswith("_"):
            options_prefix += "_"
        
        self.options_prefix = options_prefix

        if parameters == None:
            parameters = {}
        else:
            parameters = flatten_parameters(parameters)
        if cmd_high_priority:
            parameters = {k: v for k, v in parameters.items() 
                          if self.options_prefix + k not in self.commandline_options}

        self.parameters = parameters 

    def __enter__(self):
        self.opts = PETSc.Options()
        self.to_restore = {}
        self.to_delete = {}
        for k, v in self.parameters.items():
            key = self.options_prefix + k
            if key in self.opts: 
                if v != self.opts[key]:
                    self.to_restore[key] = self.opts[key]
            else:
                self.to_delete[key] = v
            self.opts[key] = v

    def __exit__(self ,type, value, traceback): 
        for k in self.to_delete.keys():
            self.opts.delValue(k)
        for k, v in self.to_restore.items():
            self.opts[k] = v


def adapt1():
    mesh = RectangleMesh(5, 5, 1, 1)
    plex = mesh.topology_dm

    plex.createLabel('adapt')
    adaptLabel = plex.getLabel('adapt')

    cs, ce = plex.getHeightStratum(0)
    for i in range(cs, ce):
        plex.setLabelValue('adapt', i, 0)
        plex.getTransitiveClosure(i)
    plex.setLabelValue('adapt', 25, 1)

    with YAOptionManager({'dm_plex_transform': {'view': None, 'type': 'refine_sbr' }}):
        plex_sbr = plex.adaptLabel('adapt')

    with YAOptionManager({'dm_plex_transform': {'view': None, 'type': 'refine_alfeld'}}):
        plex_alfeld = plex.adaptLabel('adapt')

    # plex.view()
    # plex_sbr.view()
    # plex_alfeld.view()

    m_sbr = Mesh(plex_sbr)
    m_alfeld = Mesh(plex_alfeld)

    # plot_with_cell_number(mesh)
    # plot_with_cell_number(m_sbr)
    # plot_with_cell_number(m_alfeld)

def adapt2():
    parameters = {
        'dm_adaptor': 'mmg',
        'dm_plex_metric': {
            'h_min': 0.1,
            'h_max': 0.2,
        }
    }

    with YAOptionManager(parameters):
        mesh = RectangleMesh(5, 5, 1, 1, name='test')
        plex = mesh.topology_dm
        plex.metricSetFromOptions()

        V = TensorFunctionSpace(mesh, 'CG', 1)
        V1 = FunctionSpace(mesh, 'CG', 1)
        metric = Function(V, name='metric')
        ometric = Function(metric, name='ometric')
        determinant = Function(V1, name='det')

        x, y = SpatialCoordinate(mesh)

        metric.interpolate(100*as_matrix([[1, 0], [0, 1]]))

        with metric.dat.vec_ro as vec, ometric.dat.vec as ovec, determinant.dat.vec as det:
            # vec.view()
            print(vec.max(), vec.min())
            plex.metricEnforceSPD(vec, ovec, det)
            print(ovec.max(), ovec.min())
            plex_metric = plex.adaptMetric(ovec)
        m_metric = Mesh(plex_metric)

if __name__ == '__main__':
    import signal
    import sys
    from time import sleep

    old_handler = None
    def handler(sig_num, frame):
        # global gframe
        print('Sig received with number %d'%sig_num)
        sys.exit(1)

    PETSc.Sys.popErrorHandler()
    # old_handler = signal.signal(signal.SIGSEGV, handler)
    adapt2()

