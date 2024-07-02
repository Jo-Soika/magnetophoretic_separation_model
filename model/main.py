## Magnetophoretic Separation Model
#  2023.06. Tobias Wanninger
#  from 2023.07. Johannes Soika
#  LPL - TUM
##
# %% Import
import magpylib as magpy
import magneto_sim_3d
import magneto_sim_2d
import numpy as np
import pandas as pd
import time

start = time.time()
# %% Settings
## Calculation Modes:
# 0: Separation efficiency for one V_mnp
# 1: Total separation efficiency for distributed V_mnp (V_mnp_dist)
# 2: Plot particle trajectories
calculation_mode = 0

## Study settings
activate_parameter_study = False  # Change design variable ranges below
activate_3d_dimensions = True  # 'False' 2D simulations
activate_plots = False  # 'True': displays trajectories and separation lines
activate_parabolic_flow = True  # only for 2D, if 'False': 2D averaged flow profile

## Design Variables
# Channel dimensions
h_ch = 0.0035  # [m] channel height
b_ch = 0.0035  # [m] channel width
l_ch = 0.015  # [m] channel length

# Flow Conditions
V_dot = 1e-07  # [m^3/s] Volumenstrom

# Magnet Configuration
magnet_shape = 'cylinder'  # 'cuboid', 'cylinder', (not working for parameter studies!)
polarization_val = 1040  # [mT] magnetic polarization
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

if activate_parameter_study:
    # Insert here the sweep ranges:
    h_channel_sweep = [0.0035]  # np.arange(0.0025, 0.00451, 0.001)
    b_channel_sweep = [0.0035]
    V_dot_sweep = [1.0e-7, 2.0e-7, 3.0e-7]
    V_mnp_sweep = [4.93e-18]
    polarization_val_sweep = [1000]  # [1000, 1100, 1200, 1300, 1400, 1500]

#%%
if magnet_shape == 'cuboid':  # Not for parameter studies!
    magnet_conf = (
        magpy.magnet.Cuboid(magnetization=(0, 0, polarization_val),
                            dimension=(d_mag * 1e3, b_ch * 1e3, length_magnet * 1e3),
                            position=(x_mag * 1e3, 0, -length_magnet * 1e3 / 2 - distance_mag_to_fluid)),)
else:
    magnet_conf = (magpy.magnet.Cylinder(magnetization=(0, 0, polarization_val),
                                         dimension=(d_mag * 1e3, length_magnet * 1e3),
                                         position=(x_mag * 1e3, 0, - length_magnet * 1e3 / 2 - distance_mag_to_fluid)),)
## Hyperparameters
evaluation_planes_y = 6  # num. Evaluation planes of the 3D models
number_iterations = 25  # num. iter. to find h_sep

# %% Simulations

if not activate_parameter_study:
    if activate_3d_dimensions:
        millifluidic_channel_object = magneto_sim_3d.MillifluidicChannel(
            flow_rate=V_dot,  # [m^3/s]
            channel_height=h_ch,  # [m]
            channel_width=b_ch,  # [m]
            channel_length=l_ch,
            magnetic_susceptibility=magnetic_susceptibility,
            magnetic_volume=V_mnp,  # [m^3]
            magnetic_volume_dist=V_mnp_dist,
            magnetic_volume_frac=V_mnp_frac,
            particle_density=density_mnp,
            bacterium_volume=V_bacterium,
            bacterium_density=bacterium_density,
            magnet=magnet_conf,
            evaluation_planes=evaluation_planes_y,
        )
        millifluidic_channel_object.main_run(calc_mod=calculation_mode, iterations=number_iterations,
                                             plotit=activate_plots,
                                             plot_3d_traj=activate_plots)
    elif not activate_3d_dimensions:
        millifluidic_channel_object = magneto_sim_2d.MillifluidicChannel(
            flow_rate=V_dot,  # [m^3/s]
            channel_height=h_ch,  # [m]
            channel_width=b_ch,  # [m]
            channel_length=l_ch,
            magnetic_susceptibility=magnetic_susceptibility,
            magnetic_volume=V_mnp,  # [m^3]
            magnetic_volume_dist=V_mnp_dist,
            magnetic_volume_frac=V_mnp_frac,
            particle_density=density_mnp,
            bacterium_volume=V_bacterium,
            bacterium_density=bacterium_density,
            magnet=magnet_conf,
            parabolic_flow=activate_parabolic_flow,
        )
        millifluidic_channel_object.main_run(calc_mod=calculation_mode, iterations=number_iterations,
                                             plotit=activate_plots)
elif activate_parameter_study:
    h_ch_used = []
    V_dot_used = []
    V_mnp_used = []
    polarization_val_used = []
    b_ch_used = []

    total_cb_res = []
    dpdx_res = []
    v_max = []
    loop_count = 1
    max_loops = len(h_channel_sweep) * len(V_dot_sweep) * len(V_mnp_sweep) * len(polarization_val_sweep) * len(
        b_channel_sweep)
    h_ref = h_ch

    for b_chan in b_channel_sweep:
        for h_chan in h_channel_sweep:
            for V_dot in V_dot_sweep:
                for V_mnp in V_mnp_sweep:
                    for polarization_val in polarization_val_sweep:
                        message1 = '----- Loop Nr. ' + str(loop_count) + ' of ' + str(max_loops) + ' -----'
                        print(message1)

                        V_dot_used.append(V_dot)
                        h_ch_used.append(h_chan)
                        V_mnp_used.append(V_mnp)
                        polarization_val_used.append(polarization_val)
                        b_ch_used.append(b_chan)
                        if activate_3d_dimensions:
                            millifluidic_channel_object = magneto_sim_3d.MillifluidicChannel(
                                flow_rate=V_dot,  # [m^3/s]
                                channel_height=h_chan,  # [m]
                                channel_width=b_chan,  # [m]
                                channel_length=l_ch,
                                magnetic_susceptibility=magnetic_susceptibility,
                                magnetic_volume=V_mnp,  # [m^3]
                                magnetic_volume_dist=V_mnp_dist,
                                magnetic_volume_frac=V_mnp_frac,
                                particle_density=density_mnp,
                                bacterium_volume=V_bacterium,
                                bacterium_density=bacterium_density,
                                magnet=(
                                    magpy.magnet.Cylinder(magnetization=(0, 0, polarization_val),
                                                          dimension=(d_mag * 1e3, length_magnet * 1e3),
                                                          position=(x_mag * 1e3, 0,
                                                                    - length_magnet * 1e3 / 2 - distance_mag_to_fluid)),),
                                evaluation_planes=evaluation_planes_y,
                            )
                            millifluidic_channel_object.main_run(calc_mod=calculation_mode,
                                                                 iterations=number_iterations,
                                                                 plotit=False, plot_3d_traj=False)
                        elif not activate_3d_dimensions:
                            millifluidic_channel_object = magneto_sim_2d.MillifluidicChannel(
                                flow_rate=V_dot,  # [m^3/s]
                                channel_height=h_chan,  # [m]
                                channel_width=b_chan,  # [m]
                                channel_length=l_ch,
                                magnetic_susceptibility=magnetic_susceptibility,
                                magnetic_volume=V_mnp,  # [m^3]
                                magnetic_volume_dist=V_mnp_dist,
                                magnetic_volume_frac=V_mnp_frac,
                                particle_density=density_mnp,
                                bacterium_volume=V_bacterium,
                                bacterium_density=bacterium_density,
                                magnet=(
                                    magpy.magnet.Cylinder(magnetization=(0, 0, polarization_val),
                                                          dimension=(d_mag * 1e3, length_magnet * 1e3),
                                                          position=(x_mag * 1e3, 0,
                                                                    - length_magnet * 1e3 / 2 - distance_mag_to_fluid)),),
                                parabolic_flow=activate_parabolic_flow,
                            )

                            millifluidic_channel_object.main_run(calc_mod=calculation_mode,
                                                                 iterations=number_iterations,
                                                                 plotit=False)

                        total_cb_res.append(millifluidic_channel_object.total_frac_captured_bacteria)
                        dpdx_res.append(millifluidic_channel_object.dp_dx)
                        v_max.append(millifluidic_channel_object.v_max[0])
                        loop_count += 1

    h_ch_used_ar = np.array(h_ch_used)
    V_dot_used_ar = np.array(V_dot_used)
    V_mnp_used_ar = np.array(V_mnp_used)
    polarization_val_used_ar = np.array(polarization_val_used)
    b_ch_used_ar = np.array(b_ch_used)
    total_cb_res_ar = np.array(total_cb_res)
    dpdx_res_ar = np.array(dpdx_res)
    v_max_ar = np.array(v_max)

    results = [h_ch_used_ar, b_ch_used_ar, V_dot_used_ar, V_mnp_used_ar, polarization_val_used_ar, v_max_ar,
               dpdx_res_ar,
               total_cb_res_ar]
    if activate_3d_dimensions:
        df = pd.DataFrame(
            {'h_ch': h_ch_used_ar, 'b_ch': b_ch_used_ar, 'V_dot': V_dot_used_ar, 'V_mnp': V_mnp_used_ar,
             'polarization': polarization_val_used_ar, 'v_max': v_max_ar, 'dpdx': dpdx_res_ar,
             'SE 3D:': total_cb_res_ar})
        df_name = 'results/parameter_study_3D.xlsx'
        df.to_excel(df_name, sheet_name='results3d', index=False)
        np.save('results/results_sweep_3D', np.array(results))
    elif not activate_3d_dimensions:
        df = pd.DataFrame(
            {'h_ch': h_ch_used_ar, 'b_ch': b_ch_used_ar, 'V_dot': V_dot_used_ar, 'V_mnp': V_mnp_used_ar,
             'polarization': polarization_val_used_ar, 'v_max': v_max_ar, 'dpdx': dpdx_res_ar,
             'SE 2D:': total_cb_res_ar})
        df_name = 'results/parameter_study_2D.xlsx'
        df.to_excel(df_name, sheet_name='results2d', index=False)
        np.save('results/results_sweep_2D', np.array(results))

end = time.time()
print("Simulation Time: " + str(round((end - start), 3)) + " sec")