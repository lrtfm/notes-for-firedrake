from firedrake import *
from firedrake.petsc import PETSc
from pyop2.datatypes import IntType, RealType, ScalarType, \
                            as_cstr, as_ctypes, as_numpy_dtype

import numpy as np
import matplotlib.pyplot as plt

try:
    from petsctools.options import OptionsManager
except ImportError:
    from firedrake.petsc import OptionsManager

def solve_possion(mesh, u_handle, f_handle):
    x = SpatialCoordinate(mesh)
    u_e = u_handle(x)
    f = f_handle(x)
    
    V = FunctionSpace(mesh, 'CG', 1)
    u, v = TrialFunction(V), TestFunction(V)

    L = inner(f, v)*dx
    a = inner(grad(u), grad(v))*dx
    sol = Function(V, name='u_h')
    
    bc = DirichletBC(V, u_e, 'on_boundary')

    solve(a == L, sol, bcs=bc)
    
    err = errornorm(u_e, sol, norm_type='H1')/norm(u_e, norm_type='H1')
    
    return sol, err


def estimate(mesh, sol, u_handle, f_handle, alpha, beta):
    x = SpatialCoordinate(mesh)
    u_e = u_handle(x)
    f = f_handle(x)

    V_eta_K = FunctionSpace(mesh, 'DG', 0)
    V_eta_e = FunctionSpace(mesh, 'HDivT', 0)

    phi_K = TestFunction(V_eta_K)
    phi_e = TestFunction(V_eta_e)

    phi = div(grad(sol)) + f
    g = jump(grad(sol), FacetNormal(mesh))

    ksi_K = assemble(inner(phi**2, phi_K)*dx)
    ksi_e = assemble(inner(g**2, avg(phi_e))*dS)
    ksi_outer = assemble(inner((sol-u_e)**2, phi_e)*ds)

    h_e = assemble(conj(phi_e)*ds)
    h_K = Function(V_eta_K).interpolate(CellDiameter(mesh))
    
    eta_K = assemble_eta_K_op2(ksi_K, ksi_e, ksi_outer, h_K, h_e, alpha=alpha, beta=beta)
    # eta_K2 = assemble_eta_K_py(ksi_K, ksi_e, ksi_outer, h_K, h_e, alpha=alpha, beta=beta)
    # assert np.allclose(eta_K.dat.data_ro_with_halos, eta_K2.dat.data_ro_with_halos)
    
    return eta_K


def assemble_eta_K_py(ksi_K, ksi_e, ksi_outer, h_K, h_e, alpha, beta):
    V_eta_K = ksi_K.function_space()
    V_eta_e = ksi_e.function_space()
    
    cell_node_list_K = V_eta_K.cell_node_list
    cell_node_list_e = V_eta_e.cell_node_list    
    
    ne_per_cell = V_eta_e.cell_node_list.shape[1]

    s1 = np.zeros_like(ksi_K.dat.data_ro_with_halos)
    for i in range(0, ne_per_cell):
        s1 += ksi_e.dat.data_ro_with_halos[cell_node_list_e[:, i]]
    
    s2 = np.zeros_like(ksi_K.dat.data_ro_with_halos)
    for i in range(0, ne_per_cell):
        s2 += h_e.dat.data_ro_with_halos[cell_node_list_e[:, i]] * ksi_outer.dat.data_ro_with_halos[cell_node_list_e[:, i]]

    eta_K = Function(V_eta_K)
    eta_K.dat.data_with_halos[:] = np.sqrt(
        alpha * h_K.dat.data_ro_with_halos**2 * ksi_K.dat.data_ro_with_halos \
        + beta * h_K.dat.data_ro_with_halos * s1
        # + beta * (h_K.dat.data_ro_with_halos * s1 + s2)
    )
    
    return eta_K
    

def assemble_eta_K_op2(ksi_K, ksi_e, ksi_outer, h_K, h_e, alpha, beta):
    V_eta_K = ksi_K.function_space()
    V_eta_e = ksi_e.function_space()

    kernel_str = '''
void assemble_eta_K({type} eta_K[1], {type} ksi_K[1], {type} ksi_e[{dim}], 
              {type} ksi_outer[{dim}], {type} h_K[1], {type} h_e[{dim}]) 
{{
    {type} s = 0;
    {type} s1 = 0;
    for (int i = 0; i < {dim}; i++) s += ksi_e[i];
    for (int i = 0; i < {dim}; i++) s1 += h_e[i]*ksi_outer[i];
    eta_K[0] = sqrt({alpha}*h_K[0]*h_K[0]*ksi_K[0] + {beta}*h_K[0]*s);
}}
'''.format(type=as_cstr(ScalarType), dim=V_eta_e.cell_node_list.shape[1],
           # double or complex       number of edges per element
           alpha=alpha, beta=beta)
    
    kernel = op2.Kernel(kernel_str, 'assemble_eta_K')
    
    eta_K = Function(ksi_K)
    iterset = eta_K.cell_node_map().iterset
    with PETSc.Log.Event("assemble_eta_K"):
        op2.par_loop(kernel, iterset, \
                     eta_K.dat(op2.WRITE, eta_K.cell_node_map()), \
                     ksi_K.dat(op2.READ, ksi_K.cell_node_map()), \
                     ksi_e.dat(op2.READ, ksi_e.cell_node_map()), \
                     ksi_outer.dat(op2.READ, ksi_outer.cell_node_map()), \
                     h_K.dat(op2.READ, h_K.cell_node_map()),
                     h_e.dat(op2.READ, h_e.cell_node_map())
                    )

    return eta_K


def mark_cells(mesh, eta_K, theta):
    plex = mesh.topology_dm
    cell_numbering = mesh._cell_numbering

    if plex.hasLabel('adapt'):
        plex.removeLabel('adapt')

    with eta_K.dat.vec_ro as vec:
        eta = vec.norm()
        eta_max = vec.max()[1]

    cell_node_list_K = eta_K.function_space().cell_node_list
    tol = theta*eta_max
    eta_K_data = eta_K.dat.data_ro_with_halos
    with PETSc.Log.Event("ADD_ADAPT_LABEL"):
        plex.createLabel('adapt')
        cs, ce = plex.getHeightStratum(0)
        for i in range(cs, ce):
            c = cell_numbering.getOffset(i)
            dof  = cell_node_list_K[c][0]
            if eta_K_data[dof] > tol:
                plex.setLabelValue('adapt', i, 1)

    return plex


def remove_pyop2_label(plex: PETSc.DMPlex):
    # Delete lables before create mesh with the plex, 
    # otherwise error occurs when run in parallel
    plex.removeLabel("pyop2_core")
    plex.removeLabel("pyop2_owned")
    plex.removeLabel("pyop2_ghost")
    return plex


def adapt_possion_Lshape():
    om = OptionsManager({"dm_plex_transform_type": "refine_sbr"}, options_prefix='')

    def u_exact(x):
        mesh = x.ufl_domain()
        U = FunctionSpace(mesh, 'CG', 1)
        u = Function(U)
        coords = mesh.coordinates
        x1, x2 = np.real(coords.dat.data_ro[:, 0]), np.real(coords.dat.data_ro[:, 1])
        r = np.sqrt(x1**2 + x2**2)
        theta = np.arctan2(x2, x1)
        u.dat.data[:] = r**(2/3)*np.sin(2*theta/3)
        return u
    
    def f_handle(x):
        return Constant(0)
    
    mesh = Mesh('gmsh/Lshape.msh')
    result = []
    parameters = {}
    parameters["partition"] = False

    for i in range(10):
        if i != 0:
            eta_K = estimate(mesh, sol, u_exact, f_handle, alpha=0.15, beta=0.15)
            plex = mark_cells(mesh, eta_K, theta=0.2)

            with PETSc.Log.Event("ADAPT"):
                with om.inserted_options():
                    new_plex = plex.adaptLabel('adapt')
                # Remove labels to avoid errors
                new_plex.removeLabel('adapt')
                remove_pyop2_label(new_plex)

            new_plex.viewFromOptions('-dm_view')

            # mesh = Mesh(new_plex, distribution_parameters=parameters)
            mesh = Mesh(new_plex)

        sol, err = solve_possion(mesh, u_exact, f_handle)
        ndofs = sol.function_space().dim()
        result.append((ndofs, np.real(err)))
        
    return result

def plot_adapt_result(result):
    data = np.array(result)
    if COMM_WORLD.rank == 0:
        fig, ax = plt.subplots(figsize=[6, 4])
        # plt.rcParams.update({'font.size': 12})
        ax.loglog(data[:, 0], data[:, 1], '*-', label=r'$||u_h - u||_1/||u||_1$')

        ax.loglog(data[:, 0], data[:, 0]**(-1/2)*data[-1, 1]/data[-1, 0]**(-1/2), 'r--', label='$O(N^{-1/2})$')
        ax.set_xlabel('N')
        ax.set_ylabel('H1 Rel. Err.')
        ax.legend()

        return fig

    return None


if __name__ == '__main__':
    ret = adapt_possion_Lshape()
    plot_adapt_result(ret)
