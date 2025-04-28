from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import rasterio.features
import rasterio.transform
from matplotlib.patches import Rectangle
from scipy.interpolate import NearestNDInterpolator, RectBivariateSpline
from shapely.geometry import box as sbox

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
    dst_crs=None,
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
    shiftFlatToZero = True
    tiffile = use_tiff

    latlon_bounds = (
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

    case = f"wfip_xm{abs(int(xmin))}to{int(xmax)}_ym{abs(int(ymin))}to{int(ymax)}_blendFlat3N3S3E3W_ff{flat_west}"

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
        srtm = SRTM(latlon_bounds, fpath=srtm_output, product=product)
        srtm.download()
        print(f"output tiff: {tiffile}", flush=True)
    else:
        srtm = Terrain(latlon_bounds, fpath=srtm_output, dst_crs=dst_crs)

    x1 = np.arange(xmin, xmax, ds)
    y1 = np.arange(ymin, ymax, ds)

    x, y, z = srtm.to_terrain(dx=ds)

    xref, yref = srtm.to_xy(*refloc[:2])
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
    zsrtm = interpfun(x1, y1, grid=True)
    del interpfun, z

    if write_stl and do_plot:
        fig, ax = plt.subplots(figsize=(12, 8))
        xsurf, ysurf = np.meshgrid(x1, y1, indexing="ij")
        cm = ax.pcolormesh(xsurf, ysurf, zsrtm, cmap="terrain")  # ,vmin=vmin,vmax=vmax)
        cb = fig.colorbar(cm, ax=ax)
        cb.set_label("elevation [m]", fontsize="x-large")
        ax.tick_params(labelsize="large")
        ax.set_xlabel("easting [m]")
        ax.set_ylabel("northing [m]")
        ax.set_title(f"{product} terrain height")
        ax.axis("scaled")
        fig.savefig(f"{outdir}/elevation_{product.lower()}_{case}.png", dpi=150, bbox_inches="tight")

    xsurf, ysurf = np.meshgrid(x1, y1, indexing="xy")
    zsrtm = np.flipud(zsrtm.T)
    transform = rasterio.transform.from_bounds(
        x1[0] - 0.5 * ds, y1[0] - 0.5 * ds, x1[-1] + 0.5 * ds, y1[-1] + 0.5 * ds, x1.size, y1.size
    )
    domain_area = sbox(
        xmin + (flat_west + slope_west),
        ymin + (flat_south + slope_south),
        xmax - (flat_east + slope_east),
        ymax - (flat_north + slope_north),
    )
    im = (
        rasterio.features.rasterize(
            [(domain_area, 1)], out_shape=zsrtm.shape, fill=0, transform=transform, all_touched=True
        )
        == 1
    )
    # print(100.0 * np.sum(im) / im.size, flush=True)

    interp = NearestNDInterpolator(
        list(zip(xsurf[im == 1].flatten(), ysurf[im == 1].flatten())),
        zsrtm[im == 1].flatten(),
    )
    zsrtm = interp(np.column_stack([xsurf.flatten(), ysurf.flatten()])).reshape(zsrtm.shape)

    if write_stl and do_plot:
        fig, ax = plt.subplots(figsize=(12, 8))
        cm = ax.pcolormesh(xsurf, ysurf, np.flipud(zsrtm), cmap="terrain")  # ,vmin=vmin,vmax=vmax)
        # cm = ax.imshow(zsrtm, cmap="terrain", extent=(xmin, xmax, ymin, ymax))  # ,vmin=vmin,vmax=vmax)
        cb = fig.colorbar(cm, ax=ax)
        cb.set_label("elevation [m]", fontsize="x-large")
        ax.tick_params(labelsize="large")
        ax.set_xlabel("easting [m]")
        ax.set_ylabel("northing [m]")
        ax.set_title(f"{product} terrain height (blanked)")
        ax.axis("scaled")
        fig.savefig(f"{outdir}/elevation_{product.lower()}_blanked_{case}.png", dpi=150, bbox_inches="tight")

    # ## 5. Blend surface definitions

    # check distance from west boundary
    blend_w = np.ones(xsurf.shape)
    if slope_west > 0:
        blend_w = np.maximum(
            np.minimum((xsurf - (xmin + flat_west)) / slope_west, 1),
            0,
        )

    # check distance from east boundary
    blend_e = np.ones(xsurf.shape)
    if slope_east > 0:
        blend_e = np.maximum(
            np.minimum((xmax - flat_east - xsurf) / slope_east, 1),
            0,
        )

    # check distance from south boundary
    blend_s = np.ones(xsurf.shape)
    if slope_south > 0:
        blend_s = np.flipud(
            np.maximum(
                np.minimum(
                    (ysurf - (ymin + flat_south)) / slope_south,
                    1,
                ),
                0,
            )
        )

    # check distance from north boundary
    blend_n = np.ones(xsurf.shape)
    if slope_north > 0:
        blend_n = np.flipud(
            np.maximum(
                np.minimum(
                    ((ymax - flat_north) - ysurf) / slope_north,
                    1,
                ),
                0,
            )
        )

    # combine blending functions
    blend = blend_w * blend_e * blend_s * blend_n
    if write_stl and do_plot:
        fig, ax = plt.subplots(figsize=(12, 8))
        cm = ax.pcolormesh(xsurf, ysurf, np.flipud(blend), cmap="magma")
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
        cm = ax.pcolormesh(xsurf, ysurf, np.flipud(zblend), cmap="terrain")  # ,vmin=vmin,vmax=vmax)
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
            cm = ax.pcolormesh(xsurf, ysurf, np.flipud(zblend), cmap="terrain")  # ,vmin=vmin,vmax=vmax)
            cb = fig.colorbar(cm, ax=ax)
            cb.set_label("elevation [m]", fontsize="x-large")
            ax.tick_params(labelsize="large")
            ax.set_xlabel("easting [m]")
            ax.set_ylabel("northing [m]")
            ax.set_title("shifted terrain height")
            ax.axis("scaled")
            fig.savefig(f"{outdir}/elevation_blended_{case}.png", dpi=150, bbox_inches="tight")

    # ## 6. Write out terrain surface STL
    xsurf, ysurf = np.meshgrid(x1, y1, indexing="ij")
    x1 = xsurf.flatten(order="F")
    y1 = ysurf.flatten(order="F")
    z1 = np.flipud(zblend).T.flatten(order="F")
    data = np.column_stack([x1, y1, z1])
    mesh = pv.PolyData(data)
    vtkout = f"{outdir}/terrain.vtk"
    mesh.save(vtkout)
    return xref, yref, zTerrainRef, srtm, srtm.zone_number, srtm_output
