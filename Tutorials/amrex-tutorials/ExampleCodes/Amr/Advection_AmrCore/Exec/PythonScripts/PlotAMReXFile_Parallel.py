import yt
import numpy as np
import matplotlib.pyplot as plt
import argparse
from matplotlib.patches import Rectangle
from mpi4py import MPI
import os
import glob
import time  # <-- added

def slice_and_filled_contour(plotfile, field="phi", axis="z", coord=0.0,
                             outfile_plot="contour.png",
                             nx=512, ny=512):
    """
    Extract a 2D slice from a multi-level AMReX plotfile and plot filled contours.
    Mesh boundaries are overlaid, aligned with the slice axis.
    Uses physical coordinates for FRB and mesh.
    """

    # Load dataset
    ds = yt.load(plotfile)

    # Slice
    sl = ds.slice(axis, coord)

    # Determine plane axes
    plane_axes_dict = {"x": [1, 2], "y": [0, 2], "z": [0, 1]}
    plane_axes = plane_axes_dict[axis]

    # Physical widths along plane axes
    widths = [ds.domain_right_edge[i] - ds.domain_left_edge[i] for i in plane_axes]

    # FRB in physical coordinates
    frb = sl.to_frb(widths, (nx, ny))
    arr = np.array(frb[field])

    # Coordinates along plane axes
    X = np.linspace(ds.domain_left_edge[plane_axes[0]],
                    ds.domain_right_edge[plane_axes[0]], nx)
    Y = np.linspace(ds.domain_left_edge[plane_axes[1]],
                    ds.domain_right_edge[plane_axes[1]], ny)
    X, Y = np.meshgrid(X, Y)

    # Labels
    labels_dict = {
        "x": ("y (code units)", "z (code units)"),
        "y": ("x (code units)", "z (code units)"),
        "z": ("x (code units)", "y (code units)")
    }
    xlabel, ylabel = labels_dict[axis]

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    levels = 20
    cs = ax.contourf(X, Y, arr, levels=levels, cmap="rainbow", extend="both")
    cbar = plt.colorbar(cs, ax=ax, orientation="vertical", shrink=0.95)
    cbar.set_label(field)

    # Overlay mesh boundaries
    for g in ds.index.grids:
        left = g.LeftEdge[plane_axes[0]]
        bottom = g.LeftEdge[plane_axes[1]]
        width = g.RightEdge[plane_axes[0]] - g.LeftEdge[plane_axes[0]]
        height = g.RightEdge[plane_axes[1]] - g.LeftEdge[plane_axes[1]]
        rect = Rectangle((left, bottom), width, height,
                         fill=False, edgecolor="black", linewidth=0.3)
        ax.add_patch(rect)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"Filled Contour of {field} at {axis}={coord}")
    ax.set_aspect("equal")

    plt.savefig(outfile_plot, dpi=150, bbox_inches="tight")
    plt.close()


# ========================
# MPI-enabled main
# ========================
if __name__ == "__main__":
    start_time = time.time()  # total start time

    parser = argparse.ArgumentParser(description="MPI: Plot 2D slices of AMReX data with filled contours and mesh")
    parser.add_argument("folder", help="Path to folder containing AMReX plotfiles (e.g. plt*)")
    parser.add_argument("--var", default="phi", help="Variable name to plot")
    parser.add_argument("--axis", type=str, default="z", choices=["x","y","z"],
                        help="Axis along which to slice (x, y, or z)")
    parser.add_argument("--location", type=float, default=0.0,
                        help="Coordinate along the chosen axis for the slice")
    parser.add_argument("--nx", type=int, default=512, help="Number of points along first slice axis")
    parser.add_argument("--ny", type=int, default=512, help="Number of points along second slice axis")
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Gather all plotfiles
    all_files = sorted(glob.glob(os.path.join(args.folder, "plt*")))

    if rank == 0:
        print(f"Found {len(all_files)} plotfiles")

    # Distribute files across ranks
    files_per_rank = [f for i, f in enumerate(all_files) if i % size == rank]

    # Create Images folder if not exists
    outdir = os.path.join(args.folder, "Images")
    if rank == 0 and not os.path.exists(outdir):
        os.makedirs(outdir)
    comm.Barrier()  # wait until folder exists

    # Process files assigned to this rank
    for plotfile in files_per_rank:
        step = os.path.basename(plotfile).replace("plt", "")
        outfile_plot = os.path.join(outdir, f"{args.var}_slice_{step}.png")

        t0 = time.time()
        print(f"Rank {rank} starting {os.path.basename(plotfile)}")
        slice_and_filled_contour(
            plotfile=plotfile,
            field=args.var,
            axis=args.axis,
            coord=args.location,
            outfile_plot=outfile_plot,
            nx=args.nx,
            ny=args.ny,
        )
        t1 = time.time()
        print(f"Rank {rank} finished {os.path.basename(plotfile)} in {t1-t0:.2f} seconds")

    comm.Barrier()
    if rank == 0:
        end_time = time.time()
        print(f"All ranks finished. Total elapsed time: {end_time - start_time:.2f} seconds")

