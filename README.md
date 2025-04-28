# amrTerrain
This python utility provides a workflow to create input files for immersed forcing terrain simulation in AMR-Wind. 

## Code Execution 

The code is executed by running python src/backend/backendinterface.py pathtoyamlfile

The code can be run on most versions of python which has the required python libraries available. Currently, it has been tested on python3.9 and python3.12 though conda-force. 

## YAML file 

Most of the variables used to generate the AMR-Wind input file has a default value and only minimal arguments are required if the user does not want to change the default input. A minimal yaml input file is given by 

caseParent: "/Users/hgopalan/Documents/P101_AMR-Wind/Data/tempGUI/tutorials/"
caseFolder: "1_precursor_rans"
centerLat: 36.38
centerLon: -97.40
farmRadius: 25000 
metMastLatLon: True
metMastNames: ["mast1"]
metMastLat: [36.38]
metMastLon: [-97.40]
metMastWind: [7.071,7.071]
metMastHeight: [100]

The default inputs run a no precursor simulation in RANS model for a neutral atmospheric boundary layer. 

## Input Data 

There are two optional input data from the user: (i) Terrain data and (ii) Roughness data. 

### Terrain Data

The default method for including the terrain data taps the python-elevation package. This package downloads the SRTM data at 30 or 90 m resolution. The elevation method works when the domain radius is less then 90 kms. For larger domains, the user needs to specify a till file as follows: 

useTiff: "........pathtofillfile.tif"

### Roughness Data 

Currently there is only one option for supporting heterogenous roughness data. This involves downloading the netCDF file from: https://cds.climate.copernicus.eu/datasets/satellite-land-cover?tab=overview 

This database provides a 300 m roughness data resolution. If you need a finer roughness data source, drop me an e-mail and specify the dataformat for input source. A method can be provided to include newer dataset. AMR-Wind requires input data for roughness in a specific format for efficient reading and a converter will be written inside to conver to this format. 

## Sample Workflow: 

We will start with the minimal example 

caseParent: "/Users/hgopalan/Documents/P101_AMR-Wind/Data/tempGUI/tutorials/" -- Parent folder 
caseFolder: "1_precursor_rans" -- Case folder 
centerLat: 36.38 -- Center Latitude of the Farm 
centerLon: -97.40 -- Center Longitude of the Farm 
farmRadius: 25000  -- Farm Radius 
metMastLatLon: True -- Input for met mast is in lat lon format 
metMastNames: ["mast1"] -- Name of the names. Multiple values are written as follows: ["mast1","mast2","mast3",..]
metMastLat: [36.38] -- Latitude for the met masts. Multiple values can be provided. 
metMastLon: [-97.40] -- Longitude for the met masts. Multiple values can be provided. 
metMastWind: [7.071,7.071] -- Reference wind (if used). The values are applied at the first mast in metMastNames. 
metMastHeight: [100] -- The height of the met mast above terrain in m. 

## Advanced Workflow 


### farmRadius 

The domain in AMR-Wind is always rectangular. farmRadius sets the domain to be of equal size in both Horizontal direction (X and Y). If the aspect ratio of the farm is not close to 1, the user can specify separate inouts in each direction as follows (in m): 

north: 25000 
south: 25000 
east: 25000 
west: 25000 

The domain size is set to be 25 km in each direction from centerLat and centerLon. 

### Fringe Zone 

AMR-Wind does not use terrain conforming mesh. The mesh is instead smoothed to a flat surface at the horizontal boundaries. The input for the fringe zone requires two inputs: 

northSlope: 3000 
northFlat: 1000

The first value is the distance over which the terrain smoothes to a flat surface from the existing terrain and the second value is a buffer region of the flat terrain. A similar value has to be specified for south, east and west directions. If no inputs are provided by the user, the fringe zone is set to me 10% with 5% for slope and 5% for flat. This number may be adjusted in future based on further testing. 

The vertical height of the domain is decided based on the location of the maximum terrain height in the domain and is set to: maxTerrainHt + 2 * max(maxTerrainHt,2048).

### Mesh 

AMR-Wind is a Cartesian mesh code and no mesh generation is required. The marking of the immersed forcing region to represent the terrain is done during the solver initialization stage. There are two inputs to the meshing: 

cellSize: 96 
verticalAR: 4 

The cellSize variable specifies the cell size in the horizontal direction (X and Y). verticalAR specifies the cell size in the vertical direction as cellSize/verticalAR. AMR - Wind can be run with aspect ratio of 4, 8 or 16. However, the change in the aspect ratio also increases the number of steps in the multi-grid solver at each step. It is recommended to keep the aspect ratio to 4 and include multi-level grid to refine the terrain. 

### Turbulence Model 

turbulenceModel: "LES"/"RANS"

If no turbulence model is specified, an one-equation TKE model RANS model is employed (https://link.springer.com/article/10.1023/A:1011560202388). Changing the option to "LES" enables the use of the zero-equation non-linear sub-grid scale model of Kosovic (https://journals.ametsoc.org/view/journals/mwre/138/11/2010mwr3286.1.xml). There are other turbulence models which are available in AMR - Wind and they are not supported with terrain immersed forcing method. 

### 1-D RANS Solver 

This is specified by adding the following values in the yaml file: 

rans1D: True 

The 1-D solver is optional for the LES model while it is run by default for the RANS model. The 1-D solver step runs the RANS model in single column mode and generates the vertical profile of wind speed, temperature and turbulence for initializating the simulations and using it as a boundary condition. The use of the RANS 1-D solver speeds up the solution significiantly and should be used. The solver requires the following inputs from the user: (i) reference wind speed, (ii) reference temperature, (iii) reference roughness and (iv) stability condition. 

The refinement wind speed is taken from: 

metMastWind: [7.071,7.071] 

while the temperature input (in K), roughness length (in m) and stability conditions are specified as follows:  

refTemperature: 300 
refRoughness: 0.1 
molLength: -1e30 

When these values are not specified, the optional values shown above is used. 

If a sweep angle is specified 

sweepAngle: 30 

the metMastWind speed magnitude is set to be  10 m/s. This will be discussed later. 

It is optionally possible to limit the domain height used in the 1-D solver by including the variable (in m): 

ransDomainTop: 4096 

As the 1-D solver runs fast, it is not recommended to enable this option except for testing the inversion heights and other canonical ABL effects. 

### Terrain Preprocessing 

By default, no user-input is required from the user for setting up the terrain in the solver. However, the following optional inputs can be included: 

terrainSTL: True 

This writes out the terrain as a STL file for visuvalization. 

### Terrain Refinement 

The common method in most commercial codes is to have an uniform mesh in the area of interest and to refine the mesh elsewhere. AMR - Wind does not have a general grid stretching option and uses levels to provide a refined mesh. The meshing starts with the following two inputs: 

cellSize: 96 
verticalAR: 4 

This is referred as level0 mesh. Next mesh is refined at each entries in here: 

metMastNames: ["mast1","...."]

The mesh refinement is applied in two parts: (i) cylindrical mesh around the mast and (ii) Adaptive mesh refinement (AMR) near the terrain. The cylindrical mesh around the mast is controlled with following inputs: 

metMastRadius: 500 
metMastRefinementLevel: 3 

The default refinement around each met mast creates 3 levels of refinement with a horizontal radius of 500 m. The vertical extents of the refinement runs from terrain to 100 m above the met mast and is currently fixed. AMR refinement near the terrain is also set to 3 levels of refinement. However, the user can override the AMR refinement level as follows: 

refineLow: [-3000,-3000,300]
refineHigh: [3000,3000,600]
refineTerrainMaxLevel: 3 

This is useful when there is an area of steep terrain which is not near the regions of met mast. 

In addition to the above two refinement levels, a box refinement can also be provided which requires the following inputs: 

refinementMinX: [-500,,,,,,]
refinementMaxX: [500,,,,,,]
refinementMinY: [-500,,,,,,]
refinementMaxY: [500,,,,,,]
heightAboveTerrain: [200,...]
refinementLevels: [3,...]

This capability works only in the X-Y coordinate system. It has been deprecated in future releases and it is recommended that the user specify dummy met-masts to create refinement zones which are well aligned with AMR refinements. 

### Postprocessing 

There are two different post-processing options: (i) line plots and (ii) terrain-aligned sampling. The line plots are generated at each met mast or virtual met mast locations by adding 

metMastLineSampling: True 

The terrain-aligned sampling is useful for creating speed-up maps and is enabled as follows: 

writeTerrainSampling: true 
verticalLevels: [10,80,100,200]

The distances in verticalLevels are measured from the terrain in m. A postprocessing script is included to convert the amrex format of output to vtk for plotting. 






