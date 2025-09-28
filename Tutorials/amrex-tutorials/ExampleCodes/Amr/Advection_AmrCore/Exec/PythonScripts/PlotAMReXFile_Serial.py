import yt
import numpy as np
import matplotlib.pyplot as plt
import argparse
from matplotlib.patches import Rectangle

def slice_and_filled_contour(plotfile, field="phi", axis="z", coord=0.0,
                             outfile_data="slice.npy", outfile_plot="contour.png",
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

    # Save raw data
    np.save(outfile_data, arr)
    print(f"Slice data saved to {outfile_data}")

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
    print(f"Contour plot saved to {outfile_plot}")

    return arr


# ========================
# Command-line interface
# ========================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot 2D slice of AMReX data with filled contours and mesh")
    parser.add_argument("plotfile", help="Path to AMReX plotfile")
    parser.add_argument("--var", default="phi", help="Variable name to plot")
    parser.add_argument("--axis", type=str, default="z", choices=["x","y","z"],
                        help="Axis along which to slice (x, y, or z)")
    parser.add_argument("--location", type=float, default=0.0,
                        help="Coordinate along the chosen axis for the slice")
    parser.add_argument("--nx", type=int, default=512, help="Number of points along first slice axis")
    parser.add_argument("--ny", type=int, default=512, help="Number of points along second slice axis")
    args = parser.parse_args()

    outfile_data = f"{args.var}_slice.npy"
    outfile_plot = f"{args.var}_slice.png"

    slice_and_filled_contour(
        plotfile=args.plotfile,
        field=args.var,
        axis=args.axis,
        coord=args.location,
        outfile_data=outfile_data,
        outfile_plot=outfile_plot,
        nx=args.nx,
        ny=args.ny,
    )

