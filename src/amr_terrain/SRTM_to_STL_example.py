from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import utm
from matplotlib.patches import Rectangle
from scipy.interpolate import RectBivariateSpline

from amr_terrain.terrain import SRTM, Terrain


def SRTM_Converter(
    outputDir,
    refLat,
    refLon,
    refHeight,
    left,
    right,
    bottom,
    top,
    slope_west,
    slope_east,
    slope_south,
    slope_north,
    flat_west,
    flat_east,
    flat_south,
    flat_north,
    use_tiff,
    write_stl,
    longmin,
    longmax,
    latmin,
    latmax,
    product="SRTM1",
    do_plot=True,
    ds=None,
):
    # # Terrain-Resolved Domain Setup
    # This notebook will generate a surface geometry file to be used by the microscale solver (e.g., SOWFA) to conform the solver domain to the terrain (e.g., with `moveDynamicMesh`).
    #
    # Notebook written by Eliot Quon, modified by Regis Thedin\
    # {eliot.quon,regis.thedin}@nrel.gov

    # ## 1. User input
    # Output directory (absolute dir)
    outdir = outputDir

    # The following cell should be modified by the user, following the examples given. The cell contains information about the actual location to be saved.
    #
    # The `refloc` variable is the location coresponding to (0,0) in SOWFA.
    #
    # Set the fringe on each side of the domain. This fringe region will used to blend the high-resolution SRTM domain data with either i) low-reosolution WRF (mesoscale) digital elevation model (DEM), or ii) flat.
    #
    # If `getWRFdata` above is `True`, then blending to mesoscale will occur; otherwise, the domain will be blended to flat. If blending to flat, the user can specify an extra fringe region of completely flat domain (`fringe_flat`). Additionally, if blending to flat, the terrain surface data can be shifted vertically such that the flat region is at $z=0$ by setting `shiftFlatToZero` to `True`.
    #
    # With respect to the bounding box, it is nice to have the boundaries exactly where the mesh would go because of the blending. For instance, a 5x5 km domain needs to match all the levels of cells: 20, 40, 80, 160 m. Essentially, 5000/160 needs to be an integer number, the same way 5000/80 needs to be as well. However, we only really need to match the coarsest resolution because they are multiple. Finally, an extra fringe of width `ds` (the resolution set above) is added, half on each side. The desired bounding box should go into the `xmin`,`xmax`,`ymin`,`ymax` variables, ignoring the `ds` addition. This extra fringe is to ensure that the STL is slighly larger than the bounding box that will be set on the microscale solver (needed in OpenFOAM to avoid numerical issues).

    refloc = (refLat, refLon, refHeight)
    xmin, xmax = -left, right
    ymin, ymax = -bottom, top
    fringe_flat_w = flat_west
    fringe_flat_s = flat_south
    fringe_flat_n = flat_north
    fringe_flat_e = flat_east
    shiftFlatToZero = True
    fringe_w = slope_west
    fringe_s = slope_south
    fringe_n = slope_north
    fringe_e = slope_east
    tiffile = use_tiff

    srtm_bounds = (
        refloc[1] + longmin,
        refloc[0] + latmin,
        refloc[1] + longmax,
        refloc[0] + latmax,
    )

    if tiffile == " ":
        available_products = list(SRTM.data_products.keys())
        assert product in available_products, available_products
        ds = SRTM.data_products[product]
    else:
        assert ds is not None, ds
        assert ds > 0, ds

    case = f"wfip_xm{abs(int(xmin))}to{int(xmax)}_ym{abs(int(ymin))}to{int(ymax)}_blendFlat3N3S3E3W_ff{fringe_flat_w}"

    # this will be downloaded:
    # srtm_output=f'{outdir}/{case}.tif' # need absolute path for GDAL
    # Get the absolute path needed for CGAL
    # Force to native path (Windows) using str().
    srtm_output = str(Path(f"{outdir}/{case}.tif").resolve())  # need absolute path for GDAL
    if tiffile == " ":
        pass
    else:
        srtm_output = tiffile

    if tiffile == " ":
        srtm = SRTM(srtm_bounds, fpath=srtm_output, product=product)
        srtm.download()
        print(f"output tiff: {tiffile}", flush=True)
    else:
        srtm = Terrain(srtm_bounds, fpath=srtm_output)

    x1 = np.arange(xmin, xmax, ds)
    y1 = np.arange(ymin, ymax, ds)
    xsurf, ysurf = np.meshgrid(x1, y1, indexing="ij")
    x, y, z = srtm.to_terrain(dx=ds)
    xref, yref, _, _ = utm.from_latlon(*refloc[:2], force_zone_number=srtm.zone_number)
    vmin, vmax = 1500, 2500
    if refloc[0] < 0:
        yref = yref - 10000000

    if np.amin(z) < 0:
        z[z < 0] = 0

    if write_stl:
        fig, ax = plt.subplots(figsize=(12, 8))
        cm = ax.pcolormesh(x - xref, y - yref, z, cmap="terrain")  # ,vmin=vmin,vmax=vmax)
        cb = fig.colorbar(cm, ax=ax)
        cb.set_label("elevation [m]", fontsize="x-large")
        ax.tick_params(labelsize="large")
        ax.set_xlabel("easting [m]")
        ax.set_ylabel("northing [m]")
        ax.set_title(f"{product} DEM projection")
        ax.axis("scaled")

        # bounding box for microscale region
        les = Rectangle(
            (xmin, ymin),
            xmax - xmin,
            ymax - ymin,
            edgecolor="r",
            lw=3,
            facecolor="0.5",
            alpha=0.5,
        )
        ax.add_patch(les)

    # ### 3.1 Downscale to output grid

    interpfun = RectBivariateSpline(x[:, 0] - xref, y[0, :] - yref, z)

    # resampled SRTM data stored in 'zsrtm'
    zsrtm = interpfun(x1, y1, grid=True)

    if write_stl and do_plot:
        fig, ax = plt.subplots(figsize=(12, 8))
        cm = ax.pcolormesh(xsurf, ysurf, zsrtm, cmap="terrain")  # ,vmin=vmin,vmax=vmax)
        cb = fig.colorbar(cm, ax=ax)
        cb.set_label("elevation [m]", fontsize="x-large")
        ax.tick_params(labelsize="large")
        ax.set_xlabel("easting [m]")
        ax.set_ylabel("northing [m]")
        ax.set_title(f"{product} terrain height")
        ax.axis("scaled")

        fig.savefig(f"{outdir}/elevation_{product.lower()}_{case}.png", dpi=150, bbox_inches="tight")

    # ## 4. Get the low-resolution terrain from the mesoscale
    # This part is only relevant if the user chose to blen the high-resolution SRTM terrain data with WRF

    # ## 5. Blend surface definitions

    # check distance from west boundary
    blend_w = np.ones(xsurf.shape)
    if fringe_w > 0:
        blend_w = np.minimum(np.maximum((xsurf - xmin - fringe_flat_w - fringe_w) / fringe_w, 0), 1)

    # check distance from east boundary
    blend_e = np.ones(xsurf.shape)
    if fringe_e > 0:
        blend_e = np.minimum(np.maximum((xmax - xsurf - fringe_flat_e - fringe_w) / fringe_e, 0), 1)

    # check distance from south boundary
    blend_s = np.ones(xsurf.shape)
    if fringe_s > 0:
        blend_s = np.minimum(np.maximum((ysurf - ymin - fringe_flat_s - fringe_s) / fringe_s, 0), 1)

    # check distance from north boundary
    blend_n = np.ones(xsurf.shape)
    if fringe_n > 0:
        blend_n = np.minimum(np.maximum((ymax - ysurf - fringe_flat_n - fringe_n) / fringe_n, 0), 1)

    # combine blending functions
    blend = blend_w * blend_e * blend_s * blend_n

    if write_stl and do_plot:
        fig, ax = plt.subplots(figsize=(12, 8))
        cm = ax.pcolormesh(xsurf, ysurf, blend, cmap="magma")
        cb = fig.colorbar(cm, ax=ax)
        ax.tick_params(labelsize="large")
        ax.set_xlabel("easting [m]")
        ax.set_ylabel("northing [m]")
        ax.set_title("blending function")
        ax.axis("scaled")

    # create flat surface to be blended
    # SRTM data is unlikely to be around the z=0 mark, so get the average
    z0 = np.amin(zsrtm)  # 0 #np.mean(zsrtm)
    zflat = np.full(zsrtm.shape, z0)

    # surface to blend
    zlowres = zflat

    # now, blend the high/low resolution elevations
    zblend = blend * zsrtm + (1 - blend) * zlowres

    if write_stl and do_plot:
        fig, ax = plt.subplots(figsize=(12, 8))
        cm = ax.pcolormesh(xsurf, ysurf, zblend, cmap="terrain")  # ,vmin=vmin,vmax=vmax)
        cb = fig.colorbar(cm, ax=ax)
        cb.set_label("elevation [m]", fontsize="x-large")
        ax.tick_params(labelsize="large")
        ax.set_xlabel("easting [m]")
        ax.set_ylabel("northing [m]")
        ax.set_title("blended terrain height")
        ax.axis("scaled")
        fig.savefig(f"{outdir}/elevation_blended_{case}.png", dpi=150, bbox_inches="tight")

    # # 6. Shift terrain
    # Shifts the terrain data so that the flat borders are at $z=0$

    if shiftFlatToZero:
        zTerrainRef = zblend[0, 0]
        zblend = zblend - zblend[0, 0]
        case = case + "_flatz0"

    if shiftFlatToZero:
        if write_stl and do_plot:
            fig, ax = plt.subplots(figsize=(12, 8))
            cm = ax.pcolormesh(xsurf, ysurf, zblend, cmap="terrain")  # ,vmin=vmin,vmax=vmax)
            cb = fig.colorbar(cm, ax=ax)
            cb.set_label("elevation [m]", fontsize="x-large")
            ax.tick_params(labelsize="large")
            ax.set_xlabel("easting [m]")
            ax.set_ylabel("northing [m]")
            ax.set_title("shifted terrain height")
            ax.axis("scaled")
            fig.savefig(f"{outdir}/elevation_blended_{case}.png", dpi=150, bbox_inches="tight")

    # ## 6. Write out terrain surface STL
    x1 = xsurf.flatten(order="F")
    y1 = ysurf.flatten(order="F")
    z1 = zblend.flatten(order="F")
    data = np.column_stack([x1, y1, z1])
    mesh = pv.PolyData(data)
    vtkout = f"{outdir}/terrain.vtk"
    mesh.save(vtkout)
    return xref, yref, zTerrainRef, srtm, srtm.zone_number, srtm_output
