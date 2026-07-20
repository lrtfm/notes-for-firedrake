"""Temporal and spatial convergence tests for the Cahn--Hilliard example.

Run this file from the repository root with

    python firedrake/py/cahn_hilliard_convergence.py

The notebook may instead pass its own ``solve_cahn_hilliard`` function to
``run_convergence_tests`` so that the displayed solver remains the single
implementation used by the chapter examples.
"""

from pathlib import Path

from firedrake import (
    Constant,
    Function,
    FunctionSpace,
    NonlinearVariationalProblem,
    NonlinearVariationalSolver,
    PeriodicRectangleMesh,
    SpatialCoordinate,
    TestFunctions,
    assemble,
    cos,
    div,
    dx,
    exp,
    grad,
    inner,
    pi,
    split,
    sqrt,
)
from firedrake.petsc import PETSc
import numpy as np


def solve_cahn_hilliard(
    mesh, num_steps, initial, degree=1, T=0.02,
    epsilon_value=0.05, source=None, u_exact=None, time=None,
):
    """Minimal solver used when this file is run as a standalone program."""
    V = FunctionSpace(mesh, "CG", degree)
    W = V * V
    dt_value = T / num_steps
    dt = Constant(dt_value)
    t = Constant(0.0) if time is None else time
    epsilon = Constant(epsilon_value)

    state = Function(W)
    u, mu = split(state)
    v, q = TestFunctions(W)
    state_old = Function(W)
    u_old, _ = state_old.subfunctions
    u_old.interpolate(initial)

    residual = (
        (u - u_old) / dt * v * dx
        + inner(grad(mu), grad(v)) * dx
        + mu * q * dx
        - epsilon**2 * inner(grad(u), grad(q)) * dx
        - (u**3 - u_old) * q * dx
    )
    if source is not None:
        residual -= source * v * dx(degree=8)

    problem = NonlinearVariationalProblem(residual, state)
    solver = NonlinearVariationalSolver(
        problem,
        solver_parameters={"snes_error_if_not_converged": True},
    )

    for step in range(num_steps):
        t.assign((step + 1) * dt_value)
        state.assign(state_old)
        solver.solve()
        state_old.assign(state)

    return {"mesh": mesh, "u": u_old, "exact": u_exact}


def solve_manufactured(solver, N, num_steps, degree, T):
    """Solve the manufactured-solution problem on an ``N`` by ``N`` mesh."""
    mesh = PeriodicRectangleMesh(N, N, 1.0, 1.0)
    t = Constant(0.0)
    epsilon_value = 0.05
    x, y = SpatialCoordinate(mesh)
    u_exact = exp(-t) * cos(2 * pi * x) * cos(2 * pi * y)
    mu_exact = -epsilon_value**2 * div(grad(u_exact)) + u_exact**3 - u_exact
    source = -u_exact - div(grad(mu_exact))
    return solver(
        mesh,
        num_steps,
        u_exact,
        degree=degree,
        T=T,
        epsilon_value=epsilon_value,
        source=source,
        u_exact=u_exact,
        time=t,
    )


def l2_error(result):
    """Return the final-time L2 error from a solver result dictionary."""
    error = result["u"] - result["exact"]
    return sqrt(assemble(error**2 * dx(domain=result["mesh"], degree=5)))


def convergence_rates(errors, scales):
    """Compute consecutive rates for errors measured at decreasing scales."""
    rates = [np.nan]
    for previous, current, coarse, fine in zip(
        errors[:-1], errors[1:], scales[:-1], scales[1:]
    ):
        rates.append(np.log(previous / current) / np.log(coarse / fine))
    return np.asarray(rates)


def _time_convergence(solver):
    T = 0.0025
    steps = np.asarray([4, 8, 16, 32])
    scales = T / steps
    errors = np.asarray([
        l2_error(solve_manufactured(solver, 16, n, degree=4, T=T))
        for n in steps
    ])
    rates = convergence_rates(errors, scales)
    return {
        "title": "Temporal convergence",
        "headers": ("steps", "dt", "L2 error", "rate"),
        "rows": [
            (n, f"{dt:.3e}", f"{error:.3e}",
             "--" if np.isnan(rate) else f"{rate:.2f}")
            for n, dt, error, rate
            in zip(steps, scales, errors, rates)
        ],
        "scales": scales,
        "errors": errors,
        "rates": rates,
        "order": 1,
        "xlabel": r"\tau",
    }


def _space_convergence(solver):
    T = 0.001
    N_values = np.asarray([16, 24, 32, 48])
    scales = 1.0 / N_values
    steps = np.ceil(10 * T / scales**2).astype(int)
    errors = np.asarray([
        l2_error(solve_manufactured(solver, N, n, degree=1, T=T))
        for N, n in zip(N_values, steps)
    ])
    rates = convergence_rates(errors, scales)
    return {
        "title": "Spatial convergence",
        "headers": ("N", "steps", "h", "L2 error", "rate"),
        "rows": [
            (N, n, f"{h:.3e}", f"{error:.3e}",
             "--" if np.isnan(rate) else f"{rate:.2f}")
            for N, n, h, error, rate
            in zip(N_values, steps, scales, errors, rates)
        ],
        "scales": scales,
        "errors": errors,
        "rates": rates,
        "order": 2,
        "xlabel": "h",
    }


def run_convergence_tests(solver=solve_cahn_hilliard):
    """Run both tests, optionally using a solver supplied by the notebook."""
    return {
        "time": _time_convergence(solver),
        "space": _space_convergence(solver),
    }


def print_convergence_result(result):
    """Print one convergence table with aligned headings and columns."""
    rows = [[str(value) for value in row] for row in result["rows"]]
    widths = [
        max(len(header), *(len(row[j]) for row in rows))
        for j, header in enumerate(result["headers"])
    ]
    lines = [result["title"]]
    lines.append("  ".join(
        header.rjust(width)
        for header, width in zip(result["headers"], widths)
    ))
    lines.extend(
        "  ".join(value.rjust(width) for value, width in zip(row, widths))
        for row in rows
    )
    PETSc.Sys.Print("\n".join(lines))


def plot_convergence_results(results):
    """Plot the time and space errors together with reference slopes."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.2), constrained_layout=True)
    for ax, result in zip(axes, (results["time"], results["space"])):
        scales = result["scales"]
        errors = result["errors"]
        order = result["order"]
        xlabel = result["xlabel"]
        reference = errors[-1] * (scales / scales[-1])**order
        ax.loglog(scales, errors, "o-", label="error")
        ax.loglog(scales, reference, "--", label=fr"$O({xlabel}^{order})$")
        ax.set_xlabel(fr"${xlabel}$")
        ax.set_ylabel("$L^2$ error")
        ax.invert_xaxis()
        indices = [0, len(scales) // 2, len(scales) - 1]
        ticks = scales[indices]
        ax.set_xticks(ticks, [f"{value:.1e}" for value in ticks])
        ax.tick_params(axis="x", which="minor", labelbottom=False)
        ax.grid(which="both", alpha=0.3)
        ax.legend()
    return fig


def main():
    import matplotlib

    matplotlib.use("Agg")
    results = run_convergence_tests()
    print_convergence_result(results["time"])
    print_convergence_result(results["space"])
    output = Path(PETSc.Options().getString(
        "output", "cahn_hilliard_convergence.pdf"
    ))
    fig = plot_convergence_results(results)
    fig.savefig(output)
    PETSc.Sys.Print(f"Saved convergence plot to {output}")


if __name__ == "__main__":
    main()
