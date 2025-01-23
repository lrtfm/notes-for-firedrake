from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
import json

from intro_utils import get_mesh_size, NumpyEncoder, plot_errors


_default_exact = lambda x: sin(pi*(x[0]**2 + x[1]**2))


def solve_possion_equation(mesh, degree, exact_handle=None):
    V = FunctionSpace(mesh, 'CG', degree=degree)
    u, v = TrialFunction(V), TestFunction(V)
    a = dot(grad(u), grad(v))*dx

    if exact_handle is None:
        exact_handle = _default_exact

    x = SpatialCoordinate(mesh)
    u_exact = exact_handle(x)
    f = - div(grad(u_exact))
    L = f*v*dx

    bc = DirichletBC(V, g=0, sub_domain='on_boundary')
    sol = Function(V, name='u_h')

    solve(a == L, sol, bcs=bc, options_prefix='possion')

    return sol


def points2bdy(points):
    r = np.linalg.norm(points, axis=1).reshape([-1, 1])
    return points/r


def make_high_order_mesh_map_bdy(m, p):
    coords = m.coordinates
    V_p = VectorFunctionSpace(m, 'CG', p)
    coords_p = Function(V_p, name=f'coords_p{i}').interpolate(coords)

    bc = DirichletBC(V_p, 0, 'on_boundary')
    points = coords_p.dat.data_ro_with_halos[bc.nodes]
    coords_p.dat.data_with_halos[bc.nodes] = points2bdy(points)

    return Mesh(coords_p)


def make_high_order_mesh_simple(m, p):
    if p == 1:
        return m

    coords_1 = m.coordinates
    coords_i = coords_1
    for i in range(2, p+1):
        coords_im1 = coords_i
        V_i = VectorFunctionSpace(m, 'CG', i)
        bc = DirichletBC(V_i, 0, 'on_boundary')
        coords_i = Function(V_i, name=f'coords_p{i}').interpolate(coords_im1)
        coords_i.dat.data_with_halos[bc.nodes] = points2bdy(coords_i.dat.data_ro_with_halos[bc.nodes])

    return Mesh(coords_i)


def test_possion_circle(refine, degree, exact_handle=None, iso=False, only_move_bdy=False):
    meshes = []
    for i in range(refine + 1):
        meshes.append(UnitDiskMesh(3+i))

    if exact_handle is None:
        exact_handle = _default_exact

    if iso and degree > 1:
        new_meshes = []
        for mesh in meshes:
            if only_move_bdy:
                new_mesh = make_high_order_mesh_map_bdy(mesh, p=degree)
            else:
                new_mesh = make_high_order_mesh_simple(mesh, p=degree)
            new_meshes.append(new_mesh)
        meshes = new_meshes

    sols = []
    for msh in meshes:
        sol = solve_possion_equation(msh, degree, exact_handle=exact_handle)
        sols.append(sol)

    errors_H1 = []
    errors_L2 = []
    hs = []
    for i, sol in enumerate(sols[:-1]):
        mesh = sol.function_space().mesh()
        x = SpatialCoordinate(mesh)
        u_exact = exact_handle(x)

        hs.append(get_mesh_size(mesh))
        errors_H1.append(errornorm(u_exact, sol, norm_type='H1')/norm(u_exact, norm_type='H1'))
        errors_L2.append(errornorm(u_exact, sol, norm_type='L2')/norm(u_exact, norm_type='L2'))

    hs = np.array(hs)
    errors_H1 = np.array(errors_H1)
    orders_H1 = np.log(errors_H1[1:]/errors_H1[:-1])/np.log(hs[1:]/hs[:-1])
    errors_L2 = np.array(errors_L2)
    orders_L2 = np.log(errors_L2[1:]/errors_L2[:-1])/np.log(hs[1:]/hs[:-1])
    data = {'hs': hs, 'errors_L2': errors_L2, 'orders_L2': orders_L2, 'errors_H1': errors_H1, 'orders_H1': orders_H1}

    return data


if __name__ == '__main__':
    import os, sys
    # print(os.getpid(), sys.argv)
    # for runing within jupyter notebook by magic command %run
    argstr = ''
    for arg in sys.argv[1:]:
        argstr += " "
        if " " in arg:
            if '"' not in arg:
                argstr += '"' + arg + '"'
            elif "'" not in arg:
                argstr += "'" + arg + "'"
            else:
                raise Exception("Can not process the args")
        else:
            argstr += arg

    opts = PETSc.Options().create()
    opts.insertString(argstr)
    refine = opts.getInt('refine', default=3)
    max_degree = opts.getInt('max_degree', default=3)
    exact_str = opts.getString('exact', default="1 - (x[0]**2 + x[1]**2)**4")
    # opts.view()
    # iso = opts.getBool('iso', default=False)
    # only_move_bdy = opts.getBool('only_move_bdy', default=False)

    # exact_handle = lambda x: sin(pi*(x[0]**2 + x[1]**2))
    # exact_handle = lambda x: 1 + cos(pi*(x[0]**2 + x[1]**2))
    exact_handle = lambda x: eval(exact_str)
    PETSc.Sys.Print("Exact solution: ", exact_str, '\n')

    for degree in range(1, max_degree+1):
        iso_list = [False, True] if degree > 1 else [False]
        for iso in iso_list:
            move_bdy_list = [True, False] if (iso and degree > 2) else [False]
            for only_move_bdy in move_bdy_list:
                data = test_possion_circle(refine, degree,
                                           exact_handle=exact_handle,
                                           iso=iso, only_move_bdy=only_move_bdy)

                PETSc.Sys.Print(f'p = {degree}; Use iso: {iso}; Only move bdy: {only_move_bdy}.')
                PETSc.Sys.Print(f'    Rel. H1 errors: {data["errors_H1"]}')
                PETSc.Sys.Print(f'            orders: {data["orders_H1"]}')
                PETSc.Sys.Print(f'    Rel. L2 errors: {data["errors_L2"]}')
                PETSc.Sys.Print(f'            orders: {data["orders_L2"]}\n')
                # PETSc.Sys.Print('data = ', json.dumps(data, indent=2, cls=NumpyEncoder))

                if COMM_WORLD.rank == 0:
                    string = f'degree{degree}-iso{iso}-onlyMoveBdy{only_move_bdy}'
                    np.save(f'data/errors-circle-{string}.npy', data)
                    # plot_errors(data['hs'], data['errors'], degree+1, filename=f'figures/errors-circle-{string}.pdf')
