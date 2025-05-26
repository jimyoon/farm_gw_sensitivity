# Yoon_et_al_2025_Nature_Water

## The persistence of groundwater depletion in the United States across widely ranging futures
Jim Yoon <sup>1</sup>, Stephen B. Ferencz <sup>1</sup>, Travis Thurber <sup>1</sup> 
1. Earth Systems Science Division, Pacific Northwest National Laboratory, Richland, WA, USA

Abstract: Groundwater resources are essential for crop production, enabling irrigation in many arid regions of the world. With a growing dependence on groundwater, concerns have been raised regarding alarming rates of groundwater depletion at the global scale. Modeling efforts focused on projecting future groundwater depletion have been limited by key factors: exogenous assumptions of human demand for groundwater that are irresponsive to changing groundwater availability and cost, simplistic treatment of the physics of groundwater depletion, narrow application to localized contexts, and limited exploration of future uncertainties. To address these gaps, we introduce a new modeling approach for process-rich, broad-exploration of groundwater depletion futures, initially deployed for the continental United States. Our results indicate persistent groundwater depletion across a range of hydrologic-economic futures, with many of the key agricultural regions of the U.S. encountering substantial depletion even in the absence of adverse hydrologic or economic change.

### Input Data 
1. **Farm agent based model (ABM) inputs**. Source for 
2. **Water Table Depth**. Fan, Y., Li, H., & Miguez-Macho, G. (2013). Global Patterns of Groundwater Table Depth. Science, 339(6122), 940-943. https://doi.org/10.1126/science.1229881
3. **Porosity and Permeabilty**. Gleeson, T., Moosdorf, N., Hartmann, J., & van Beek, L. P. H. (2014). A glimpse beneath earth's surface: GLobal HYdrogeology MaPS (GLHYMPS) of permeability and porosity. Geophysical Research Letters, 41(11), 3891-3898. https://doi.org/10.1002/2014GL059856
4. **Recharge**. Döll, P., & Fiedler, K. (2008). Global-scale modeling of groundwater recharge. Hydrol. Earth Syst. Sci., 12(3), 863-885. https://doi.org/10.5194/hess-12-863-2008
5. **Aquifer Depth** de Graaf, Inge, Laura Condon, and Reed Maxwell. "Hyper‐resolution continental‐scale 3‐D aquifer parameterization for groundwater modeling." Water Resources Research 56.5 (2020): e2019WR026004. 
6. **Depth to Bedrock** Shangguan, Wei, et al. "Mapping the global depth to bedrock for land surface modeling." Journal of Advances in Modeling Earth Systems 9.1 (2017): 65-88.

### Output Data 
3.  Yoon, J., Ferencz, S. (2025). Name for Data Repo (Version v1) [Data set]. MSD-LIVE Data Repository. DOI

## Code Reference 
Code for Executing the processing and analysis steps in "Reproduce my Experiment" provided in the `workflow` folder on this meta-repository. Yoon, J., Ferencz, S., and Thurber, T. (2025). Name for Code Repo (Version X.X). Zenodo. DOI 

## Contributing Modeling Software 
-Standard Python Packages along with the Pyomo package (https://www.pyomo.org/). Pyomo is a Python-based, open-source optimization modeling language.  
-IPOPT (https://coin-or.github.io/Ipopt/). Ipopt (Interior Point Optimizer) is an open source software package for large-scale nonlinear optimization.   
-QGIS (https://qgis.org/) for geospatial processing and visualization.   
-Superwell V1.1 (https://zenodo.org/records/14583794). Modified version of the Supwerwell V1.1 Python code is used to generate the groundwater cost curves.   

## Recreate my experiment 
# creating the environment:
- `module load python/miniconda4.12`
- `source /share/apps/python/miniconda4.12/etc/profile.d/conda.sh`
- `conda create -n farmwell python=3.11`
- `conda activate farmwell`
- `conda install -c conda-forge numpy pandas scipy SALib pyomo matplotlib openpyxl ipopt`
- `conda deactivate`

# running the experiment
- Update the total number of nodes to use in `q.sl`
- Update the expected walltime in `q.sl`
- Update the partition in `q.sl`
- Update the total number of farms to run in `driver.sh`
- `sbatch q.sl`
- results and errors are written to `./output` directory per farm

# groundwater farm ABM workflow 
- input files 

# processing output files 
- use

# reproduce my figures 
- Figure 1
- Figure 2
- Figure 3
- Figure 4 
