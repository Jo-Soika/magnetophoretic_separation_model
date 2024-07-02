# Magnetophoretic Separation Model
`magnetophoretic_separation_model` is a Python package for estimating separation efficiency of 2D and 3D magnetophoretic separation systems.
It uses analytical descriptions for particle motion, velocity profile within the separation channel and magnetic field.

***
## Code
### Installation
After employing the file and create a local environment, run `setup.py`. Make sure that the directory contains `requirements.txt`.

### First run
Run the `main.py` code in `model/` folder. The initial settings represent the separation efficiency of a specific separation channel.
The results shall be displayed, including a 3D visualization of the trajectories, with attracted (blue) and escaping (red) trajectories.


### Change the Design Variables
The design variables of the channel can be changed in the `main.py` code. For that, jump to the `##Design Variables` section:
```
## Design Variables
# Channel dimensions
h_ch = 0.0035  # [m] channel height
b_ch = 0.0035  # [m] channel width
l_ch = 0.015  # [m] channel length

# Flow Conditions
V_dot = 1e-07  # [m^3/s] Volumenstrom

# Magnet Configuration
magnet_shape = 'cylinder'  # 'cuboid', 'cylinder', (not working for parameter studies!)
polarization_val = 1500  # [mT] magnetic polarization
length_magnet = 0.01  # [m] length of magnet (in polarization dir.)
d_mag = b_ch  # [m] size magnet (diameter or side length)
distance_mag_to_fluid = 0.0001  # [m] nearest distance of magnet to fluid within channel
x_mag = l_ch - d_mag / 2  # [m] position of magnet in x-dir

# Magnetic nano particle properties / beads properties
V_mnp = 4.93E-18  # [m^3] volume magnet particle (for calc_mode = 0)
V_mnp_dist = [7.40231E-18, 1.16425E-17, 1.8308E-17, 2.87763E-17, 4.52746E-17, 7.1186E-17]
V_mnp_frac = [0.04244, 0.182, 0.319, 0.2898, 0.1386, 0.02818]
density_mnp = 5240  # [kg/m^3]
magnetic_susceptibility = 2.8

# Non-magnetic particle properties (e.g. bacterium)
V_bacterium = 0  # [m^3] volume bacterium (if 0 -> no bacterium)
bacterium_density = 1100  # [kg/m^3]
```
Find more information about the definitions of the variables, the setup and the computation method in our publication:
`DOI............`

### Change calculation mode
In addition to determining the separation efficiency, the code has other functions. Certain starting positions can be selected for a mode so that the behavior of the particles is displayed in different areas.
The calculation mode can be changed in the `main.py` code:
```
calculation_mode = 0
## Calculation Modes:
# 0: Separation efficiency for one V_mnp
# 1: Total separation efficiency for distributed V_mnp (V_mnp_dist)
# 2: Plot particle trajectories
```
Mode 0 evaluates the separation efficiency of one magnetic nanoparticle volume.

Mode 1 evaluates the separation efficiency of different magnetic nanoparticle volumes. This is usefull to take different agglomerate sizes into account. Change `V_mnp_dist` for the volume of the agglomerates and `V_mnp_frac` to set the fraction of MNP with this specific volume to the total amount of MNP.

Mode 2 plots trajectories of particles released from a grid.

### Study settings
Besides a 3D rectangular channel, the model can simulate 2D systems. Additionally, parameter studies can be done, that, however, need to be defined explicitly.
It can be activated if the trajectories shall be plotted. For parameter studies, the plots shall be deactivated.
For that, turn the following settings to `True` or `False`.
```
## Study settings
activate_parameter_study = False  # Change design variable ranges below
activate_3d_dimensions = True  # 'False' 2D simulations
activate_plots = False  # 'True': displays trajectories and separation lines
activate_parabolic_flow = True  # only for 2D, if 'False': 2D averaged (constant) flow profile
```
### Parameter studies
As mentioned, parameter studies are possible. For the studies, the plots are deactivated by default.
The model considers five different design variables to sweep. The ranges can be modified here: 
```
# Insert here the sweep ranges:
h_channel_sweep = [0.0035]
V_dot_sweep = [1.0e-7]
V_mnp_sweep = [4.93e-18]
polarization_val_sweep = [1000, 1100, 1200, 1300, 1400, 1500]
b_channel_sweep = [0.0035]
```
It sweeps through every possible combination, so that total simulation time can be very large. For large sweeps, the 2D models are preferred as simulation are faster than in 3D.


***
## Abstract of the accompanying publication
Continuous flow magnetophoresis represents a common technique for actively separating particles within a fluid. For separation systems design, accurately predicting particle behaviour helps to characterise system performance, typically measured by the separation efficiency (SE). While Finite Element Method (FEM) simulations offer high accuracy, they demand extensive computational resources. Alternatively, results can be achieved more quickly with simplified numerical models that use analytical descriptions of fluid flow, magnetic fields, and particle movement.

In this research, we model a millifluidic system that separates magnetic particles using magnetophoresis. Therefore, we (1) develop a simple numerical model that can simulate continuous flow separation for rectangular channels in two and three dimensions, (2) introduce a novel and simple approach to calculate the SE, and (3) quantify the effects of model assumptions in flow profile and dimensions on SE. Our method for estimating SE considers particle flux variation across the channel's cross-section due to the flow profile. The results are compared to a FEM model developed in COMSOL.

The obtained three-dimensional simulation model computes results in seconds, around 180 times faster than the FEM approach, while deviating less than 2\% from the FEM model. A comparison of the different two-dimensional and three-dimensional models underscores the significant influence of the flow profile and the SE calculation method on the result. The two dimensional models generally overestimate the SE of up to 15\% due to their lower peak flow velocity. However, using an averaged flow velocity leads to good agreement for high SE due to the overlap of differences in flow profile and SE calculation.


## Citation
If you use the provided code or find it helpful for your own work, please cite:

TODO
```
@article{soika2024magnetophoretic,
  title={Efficient Numerical Modelling of Magnetophoresis in Millifluidic Systems},
  author={Soika, Johannes and Wanninger, Tobias and Muschak, Patrick and Schwaminger, Sebastian and Berensmeier, Sonja and Zimmermann, Markus},
  journal={},
  volume={},
  number={},
  pages={},
  year={2024},
  publisher={}
}
```

