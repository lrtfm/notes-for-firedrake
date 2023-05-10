from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
import json

from intro_utils import NumpyEncoder, plot_errors


def solve_possion_equation(mesh, degree):
    V = FunctionSpace(mesh, 'CG', degree=degree)
    u, v = TrialFunction(V), TestFunction(V)
    a = dot(grad(u), grad(v))*dx

    x, y = SpatialCoordinate(mesh)
    f = sin(pi*x)*sin(pi*y)
    L = f*v*dx
    
    bc = DirichletBC(V, g=0, sub_domain='on_boundary')
    sol = Function(V, name='u_h')

    solve(a == L, sol, bcs=bc, options_prefix='possion')
    
    return sol

        
def test_possion(N, refine, degree):
    base = RectangleMesh(N, N, 1, 1)
    meshes = MeshHierarchy(base, refinement_levels=refine)

    sols = []
    for msh in meshes:
        sol = solve_possion_equation(msh, degree)
        sols.append(sol)

    errors = []
    for i, sol in enumerate(sols[:-1]):
        V_ref = sols[i+1].function_space()
        sol_int = project(sol, V_ref)
        errors.append(errornorm(sol_int, sols[i+1]))

    hs = np.array([1/N/2**i for i in range(len(sols)-1)])
    errors = np.array(errors)
    orders = np.log(errors[1:]/errors[:-1])/np.log(hs[1:]/hs[:-1])
    data = {'hs': hs, 'errors': errors, 'orders': orders}
    
    return data


if __name__ == '__main__':
    opts = PETSc.Options()
    N = opts.getInt('N', default=8)
    refine = opts.getInt('refine', default=3)
    degree = opts.getInt('degree', default=1)
    
    data = test_possion(N, refine, degree)
    
    PETSc.Sys.Print('data = ', json.dumps(data, indent=2, cls=NumpyEncoder))

    if COMM_WORLD.rank == 0:
        np.save('errors.npy', data)
        plot_errors(data['hs'], data['errors'], degree+1, filename='figures/errors.pdf')
