from firedrake import *
from firedrake.petsc import PETSc
from mpi4py import MPI
import os 
import sys
import gmsh
import numpy as np
from scipy.interpolate import BarycentricInterpolator
import matplotlib.pyplot as plt


__all__ = (
    'make_circle_in_rect',  
    'make_coords_order_plus_one',   
    'make_high_order_coords_for_circle_in_rect', 
    'make_high_order_coords_for_circle_in_rect_simple', 
)

def make_rect(h, filename, p=1, gui=False):
    gmsh.initialize()

    gmsh.model.add('rect')

    L = 1
    rect = gmsh.model.occ.add_rectangle(0, 0, 0, L, L)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)

    __old_verbosity = gmsh.option.getNumber("General.Verbosity")
    gmsh.option.setNumber("General.Verbosity", 1)
    # TODO: catch the exception of gmsh
    gmsh.model.mesh.generate()
    gmsh.option.setNumber("General.Verbosity", __old_verbosity)

    filename = filename or 'test.msh'

    gmsh.write(filename) 
    if gui:
        gmsh.fltk.run()
    gmsh.finalize()


def make_circle_in_rect(h, filename, p=1, gui=False):
    gmsh.initialize()

    gmsh.model.add('circle_in_rect')

    L = 1
    r = 1/4
    rect = gmsh.model.occ.add_rectangle(0, 0, 0, L, L)
    circle = gmsh.model.occ.add_disk(L/2, L/2, 0, r, r)
    omega = gmsh.model.occ.cut([[2, rect]], [[2, circle]], removeObject=True, removeTool=False)

    gmsh.model.occ.synchronize()

    # interface = gmsh.model.get_boundary([[2, circle]])
    
    outer = []
    interface = []
    bdys = gmsh.model.get_boundary(omega[0])
    for dim, tag in bdys:
        tag = abs(tag)
        t = gmsh.model.get_type(dim, tag)
        # gmsh.model.get_value # coordinates
        if t == 'Line':
            outer.append(tag)
        else: # 'Ellipse'q
            interface.append(tag)

    gmsh.model.add_physical_group(dim=1, tags=outer, tag=1)
    gmsh.model.set_physical_name(dim=1, tag=1, name="bdy")
    gmsh.model.add_physical_group(dim=1, tags=interface, tag=2)
    gmsh.model.set_physical_name(dim=1, tag=2, name="interface")

    gmsh.model.add_physical_group(dim=2, tags=[omega[0][0][1]], tag=1)
    gmsh.model.set_physical_name(dim=2, tag=1, name="omega_plus")

    gmsh.model.add_physical_group(dim=2, tags=[circle], tag=2)
    gmsh.model.set_physical_name(dim=2, tag=2, name="omega_minus")

    if p > 1:
        factor = 1/p
        d0 = 1/8
        c0 = d0**(factor - 1)
        h_min = h**(1/factor)
        def meshSizeCallback(dim, tag, x, y, z, lc):
            r1 = x**2 + y**2
            r2 = (x-1)**2 + y**2
            r3 = (x-1)**2 + (y-1)**2
            r4 = x**2 + (y-1)**2

            if r1 < d0**2:
                if r1 > h_min**2:
                    return c0*r1**((1-factor)/2)*h
                else:
                    return h_min

            if r2 < d0**2:
                if r2 > h_min**2:
                    return c0*r2**((1-factor)/2)*h
                else:
                    return h_min

            if r3 < d0**2:
                if r3 > h_min**2:
                    return c0*r3**((1-factor)/2)*h
                else:
                    return h_min

            if r4 < d0**2:
                if r4 > h_min**2:
                    return c0*r4**((1-factor)/2)*h
                else:
                    return h_min

            return h

        gmsh.model.mesh.setSizeCallback(meshSizeCallback)
    else:
        gmsh.option.setNumber("Mesh.MeshSizeMin", h)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)

    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    __old_verbosity = gmsh.option.getNumber("General.Verbosity")
    gmsh.option.setNumber("General.Verbosity", 1)
    # TODO: catch the exception of gmsh
    gmsh.model.mesh.generate()
    gmsh.option.setNumber("General.Verbosity", __old_verbosity)

    filename = filename or 'test.msh'

    gmsh.write(filename) 
    if gui:
        gmsh.fltk.run()
    gmsh.finalize()


def coords_to_interface(v):
    c = np.array([[0.5, 0.5]])
    r = 1/4
    points = v - c
    scale = r/np.linalg.norm(points, axis=1).reshape([-1, 1])
    return c + scale*points


def get_interface_cells(m, interface=2):
    V_flag = FunctionSpace(m, 'CG', 1)
    V_flag_dg0 = FunctionSpace(m, 'DG', 0)
    flag = Function(V_flag, name='flag')
    flag_dg0 = Function(V_flag_dg0, name='interface_cell_flag') # which has two point on the interface
    bc_flag = DirichletBC(V_flag, 1, interface)
    bc_flag.apply(flag)
    flag_dg0.interpolate(flag)
    interface_cells, = np.where(flag_dg0.dat.data_with_halos > 1/2)
    non_interface_nodes_local = np.zeros_like(interface_cells)
    data = flag.dat.data_ro_with_halos[:]
    for i, cell in enumerate(interface_cells):
        nodes = V_flag.cell_node_list[cell]
        index, = np.where(data[nodes] < 1/2)
        non_interface_nodes_local[i] = index
    return interface_cells, non_interface_nodes_local


def get_bary_locs(V_coords):
    import numpy
    element = V_coords.finat_element.fiat_equivalent
    locations = []
    # The order here is the numbering in entity_dofs
    for dual in element.dual_basis():
        loc, = dual.get_point_dict().keys()
        locations.append(loc)

    locations = numpy.asarray(locations)

    # Now we can convert to barycentric coordinates.
    cell = element.ref_el
    vertices = numpy.asarray(cell.get_vertices())
    bary2cart = numpy.vstack([vertices.T, numpy.ones(len(vertices))])
    cart2bary = numpy.linalg.inv(bary2cart)

    bary_locs = numpy.hstack([locations, numpy.ones((len(locations), 1))]).dot(cart2bary.T)
    
    return bary_locs


def make_high_order_coords_for_circle_in_rect_simple(m, p):
    if p == 1:
        return m.coordinates
    coords_2 = make_high_order_coords_for_circle_in_rect(m, 2)
    if p == 2:
        return coords_2
    
    coords_i = coords_2
    for i in range(3, p+1):
        coords_i_1 = coords_i
        V_i = VectorFunctionSpace(m, 'CG', i)
        bc = DirichletBC(V_i, 0, 2)
        coords_i = Function(V_i, name=f'coords_p{i}').interpolate(coords_i_1)
        coords_i_l = Function(V_i, name=f'coords_p{i}_l').interpolate(m.coordinates)
        coords_i.dat.data_with_halos[bc.nodes] = coords_to_interface(coords_i_l.dat.data_ro_with_halos[bc.nodes])

    return coords_i


def make_high_order_coords_for_circle_in_rect(m, p):
    """
    Ref: M. Lenoir, Optimal Isoparametric Finite Elements and Error Estimates for 
        Domains Involving Curved Boundaries, SINUM Vol. 23, No. 3 (1986), pp. 562-580
    """
    if p == 1:
        return m.coordinates
    if p == 2:
        V_c = VectorFunctionSpace(m, 'CG', 2)
        bc = DirichletBC(V_c, 0, 2)
        coords = Function(V_c, name='coords_h')
        coords.interpolate(m.coordinates)
        coords.dat.data_with_halos[bc.nodes] = coords_to_interface(coords.dat.data_ro_with_halos[bc.nodes])
        return coords
    
    coords_2 = make_high_order_coords_for_circle_in_rect(m, 2)
    coords_i = coords_2
    for i in range(3, p+1):
        PETSc.Sys.Print(f'i = {i}')
        coords_i_1 = coords_i
        coords_i = make_coords_order_plus_one(coords_i_1)

    return coords_i


def make_coords_order_plus_one(coords_p_1, p=None, show=False):
    
    m = coords_p_1.ufl_domain()
    p = p or coords_p_1.ufl_element().degree()+1
    
    interface_cells, non_interface_nodes_local = get_interface_cells(m, interface=2)

    V_coords = VectorFunctionSpace(m, 'CG', p)
    entity_dofs = V_coords.finat_element.entity_dofs()
    entity_permutations = V_coords.finat_element.entity_permutations
    entities = entity_dofs[2][0]

    bary_locs = get_bary_locs(V_coords)
    bary = bary_locs[entities]

    coords_p = Function(V_coords).interpolate(coords_p_1)

    indexmap = [[1, 2], [0, 2], [0, 1]]
    new_bary = np.zeros_like(bary)
    max_shift = 0
    if show:
        labels = ['interface nodes', 'inner nodes', 'inter p', 'inter p-1']
    coords_p_data = coords_p.dat.data_ro_with_halos[:].copy()
    for i, cell in enumerate(interface_cells):
        nodes = V_coords.cell_node_list[cell]
        non_interface_node = non_interface_nodes_local[i]
        interface_nodes = indexmap[non_interface_node]
        new_bary[:] = bary[:]
        s = np.sum(bary[:, interface_nodes], axis=1).reshape(-1, 1)
        new_bary[:, non_interface_node] = 0
        new_bary[:, interface_nodes] = new_bary[:, interface_nodes]/s
        coords_local = coords_p_data[nodes[:3]]

        new_coords = new_bary@coords_local
        x = coords_local[interface_nodes, :]

        def get_new_coords_interface(p):
            coords_1d = np.arange(p+1)/p
            bary_1d = np.vstack([1-coords_1d, coords_1d]).T
            v = bary_1d@x
            v_i = coords_to_interface(v)

            inp_x_m = BarycentricInterpolator(coords_1d, v_i[:, 0])
            inp_y_m = BarycentricInterpolator(coords_1d, v_i[:, 1])

            bary_seg = new_bary[:, interface_nodes]

            new_x = inp_x_m(bary_seg[:, 1])
            new_y = inp_y_m(bary_seg[:, 1])
            new_coords_interface = np.vstack([new_x, new_y]).T

            return new_coords_interface, v_i

        new_coords_interface_p, coords_nodes_interface = get_new_coords_interface(p)
        new_coords_interface_p_1, _ = get_new_coords_interface(p-1)
        shift = s**p*(new_coords_interface_p - new_coords_interface_p_1)
        
        nis = entity_dofs[1][non_interface_node]
        coords_p_data[nodes[nis]] = coords_nodes_interface[1:-1, :] # dofs on boundary
        coords_p_data[nodes[entities]] = coords_p_data[nodes[entities]] + shift
        
        smax = np.max(np.linalg.norm(shift, axis=1))
        if smax > max_shift:
            max_shift = smax
        if show:
            center = np.sum(coords_local, axis=0)/3
            xc, yc = center
            def scale(x, y):
                return 0.7*(x-xc)+xc, 0.7*(y-yc)+yc
            plt.text(center[0], center[1], f'c')
            plt.fill(coords_local[:, 0], coords_local[:, 1], ls='-.', fill=False, alpha=0.5)

            for i, (x, y) in enumerate(coords_nodes_interface[1:-1]):
                plt.plot(x, y, 'd', fillstyle='none', label=labels[0])
                labels[0] = None
                plt.text(*scale(x, y), f'{i+1}', ha='center', ma='center')
            for i, (x, y) in enumerate(coords_p_data[nodes[:3]]):
                # plt.plot(x, y, 'x')
                plt.text(*scale(x, y), f'n{i}', ha='center', ma='center')
            for i, (x, y) in enumerate(coords_p_data[nodes[entities]]):
                plt.plot(x, y, 'x', label=labels[1])
                labels[1] = None
                # plt.text(x, y, f'{i}')
            for i, (x, y) in enumerate(new_coords_interface_p):
                plt.plot(x, y, '+', label=labels[2])
                labels[2] = None
                # plt.text(x, y, f'new{i}', rotation=60)
            for i, (x, y) in enumerate(new_coords_interface_p_1):
                plt.plot(x, y, 'o', fillstyle='none', label=labels[3])
                labels[3] = None
                # plt.text(x, y, f'm1-{i}', rotation=35)
    coords_p.dat.data_with_halos[:] = coords_p_data[:]
    if show:
        plt.legend()
        plt.title(f'p={p}')
    max_shift = m.comm.reduce(max_shift, MPI.MAX, root=0)
    PETSc.Sys.Print(f'Max shift: {max_shift}')
    return coords_p


def test_cmp(p=4, h=0.1, filename=None):
    filename = filename or f'circle_in_rect_h{h}.msh'
    make_circle_in_rect(h, filename, gui=False)
    m = Mesh('circle_in_rect_h0.2.msh')
    coords_p  = make_high_order_coords_for_circle_in_rect(m, p)
    coords_p_s = make_high_order_coords_for_circle_in_rect_simple(m, p)
    err = np.max(np.linalg.norm(coords_p.dat.data_ro_with_halos - coords_p_s.dat.data_ro_with_halos, axis=1))
    PETSc.Sys.Print(f'Difference of two method: {err}')
    return err
    

def test(h=0.1, filename=None):
    filename = filename or f'circle_in_rect_h{h}.msh'
    make_circle_in_rect(h, filename, gui=False)

    m = Mesh(filename)
    x, y = SpatialCoordinate(m)
    v0_expr = conditional((x-1/2)**2 + (y-1/2)**2 > 1/16, 0, 1)
    V_dg0 = FunctionSpace(m, 'DG', 0)
    v0 = Function(V_dg0, name='init')

    coords_2 = make_high_order_coords_for_circle_in_rect(m, 2)

    plt.figure(figsize=[8, 8])
    coords_3 = make_coords_order_plus_one(coords_2, show=True)
    plt.axis('equal')
    plt.savefig(f'{filename}_p3.pdf', bbox_inches='tight')

    plt.figure(figsize=[8, 8])
    coords_4 = make_coords_order_plus_one(coords_3, show=True)
    plt.axis('equal')
    plt.savefig(f'{filename}_p4.pdf', bbox_inches='tight')

    plt.figure(figsize=[8, 8])
    coords_5 = make_coords_order_plus_one(coords_4, show=True)
    plt.axis('equal')
    plt.savefig(f'{filename}_p5.pdf', bbox_inches='tight')
    

def test_cmp(p=4, h=0.1, filename=None):
    filename = filename or f'circle_in_rect_h{h}.msh'
    make_circle_in_rect(h, filename, gui=False)
    m = Mesh('circle_in_rect_h0.2.msh')
    coords_p  = make_high_order_coords_for_circle_in_rect(m, p)
    coords_p_s = make_high_order_coords_for_circle_in_rect_simple(m, p)
    err = np.max(np.linalg.norm(coords_p.dat.data_ro_with_halos - coords_p_s.dat.data_ro_with_halos, axis=1))
    PETSc.Sys.Print(f'Difference of two method: {err}')


if __name__ == '__main__':
    test(h=0.1)
