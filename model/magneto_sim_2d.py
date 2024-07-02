# %%
import magpylib as magpy
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import solve_ivp
import Event
from matplotlib.collections import LineCollection
from math import pi
import pickle

# %% Millifluidic class
class MillifluidicChannel:
    def __init__(
            self,
            flow_rate=None,
            channel_height=None,
            channel_width=None,
            channel_length=None,
            magnetic_susceptibility=None,
            magnetic_volume=None,
            magnetic_volume_dist=None,
            magnetic_volume_frac=None,
            bacterium_volume=4.71e-18,
            bacterium_density=1100,
            magnet=None,
            particle_density=5240,
            fluid_density=997,
            mu_0=4e-7 * pi,
            gravitational_acceleration=np.array([0, 0, -9.81]),
            mu=0.001,
            initial_state=None,
            parabolic_flow=True,
            number_of_steps=12,
    ):
        self.channel_height = channel_height
        self.channel_width = channel_width
        self.channel_length = channel_length
        self.magnetic_susceptibility = magnetic_susceptibility
        self.magnetic_volume = magnetic_volume
        self.magnet = magnet
        self.magnet_collection = magpy.Collection()
        self.particle_density = particle_density
        self.fluid_density = fluid_density
        self.mu_0 = mu_0
        self.gravitational_acceleration = gravitational_acceleration
        self.mu = mu
        self.hydraulic_drag_diameter = None
        self.particle_path = []
        self.fig1 = None
        self.ax1 = None
        self.aggregate_density = (particle_density * magnetic_volume + fluid_density * (
                1e-6 / 2) ** 2 * math.pi * 6e-6) / (magnetic_volume + (1e-6 / 2) ** 2 * math.pi * 6e-6)
        self.bacterium_volume = bacterium_volume
        self.bacterium_density = bacterium_density
        self.h_max = None
        self.F_mag = []
        self.F_g = []
        self.speed = []
        self.h_grad_h = []
        self.sol = None
        self.initial_state = initial_state
        self.flow_rate = flow_rate
        self.magnetic_volume_dist = magnetic_volume_dist
        self.magnetic_volume_frac = magnetic_volume_frac
        self.dp_dx = None
        self.frac_captured_bacteria = None
        self.total_frac_captured_bacteria = None
        self.max_step = None
        self.v_max = None
        self.parabolic_flow = parabolic_flow
        self.n_step = number_of_steps

    def main_run(self, calc_mod=None, iterations=25, plotit=False, msg=True):
        self.magnet_collection = magpy.Collection()
        for n in range(len(self.magnet)):
            self.magnet_collection.add(self.magnet[n])

        self.v_max = np.round(self.flow_model([0, 0, self.channel_height / 2]), 4)
        print("max flow velocity: " + str(self.v_max) + " m/s")

        if calc_mod == 0:
            self.hydraulic_drag_diameter = (6 / pi * (self.magnetic_volume + 0 * self.bacterium_volume)) ** (1 / 3)
            if plotit is True:
                self.fig1, self.ax1 = plt.subplots()
                self.plot_magnetic_field()
            self.h_max = self.determine_h_max(iterations, plotit)
            if msg:
                print('Determined h_max: {}'.format(round(self.h_max, 4)))
            if plotit:
                plt.tight_layout()  #
                plt.show()
            self.frac_captured_bacteria = self.captured_bacteria()
            self.total_frac_captured_bacteria = self.frac_captured_bacteria
            if msg:
                print('Total fraction of captured bacteria: ' + str(round(self.frac_captured_bacteria * 100, 4)) + '%')

        elif calc_mod == 1:
            # calculates the (total) bacterium loss and separation height; V_mag is different for every loop;
            self.frac_captured_bacteria = []
            if plotit is True:
                # Create a figure before entering the loop
                self.fig1, ax = plt.subplots(figsize=(10, len(self.magnetic_volume_dist) * 5))
            for ind_vmag, V_magnetic in enumerate(self.magnetic_volume_dist):
                self.particle_path = []
                self.hydraulic_drag_diameter = (6 / pi * (V_magnetic + self.bacterium_volume)) ** (1 / 3)
                self.h_max = []
                self.magnetic_volume = V_magnetic
                if plotit is True:
                    self.ax1 = plt.subplot(len(self.magnetic_volume_dist), 1, ind_vmag + 1)
                    self.plot_magnetic_field()
                self.h_max = self.determine_h_max(iterations, plotit)
                self.frac_captured_bacteria.append(self.captured_bacteria())
                print('SE = ' + str(round(self.captured_bacteria() * 100, 4)) + '%')
            if plotit is True:
                plt.tight_layout()  # This ensures the plots do not overlap
                plt.show()  # Display the plots
            self.total_frac_captured_bacteria = sum(np.multiply(self.frac_captured_bacteria, self.magnetic_volume_frac))
            if msg:
                print('The total Separation Efficiency is: SE = ' + str(
                    round(self.total_frac_captured_bacteria * 100, 4)) + '%')
        elif calc_mod == 2:
            z_n = np.linspace(self.channel_height * 0.02, self.channel_height * 0.98, 20)
            initial_particle_position = []
            for j, j_n in enumerate(z_n):
                initial_particle_position.append((0.0, 0.0, j_n))
            initial_particle_position = np.array(initial_particle_position)
            # initial_particle_position = ((0.0, 0.0, 0.003), (0.0, 0.0, 0.00175), (0.0, 0.0, 0.0005),)

            self.hydraulic_drag_diameter = (6 / pi * (self.magnetic_volume + self.bacterium_volume)) ** (1 / 3)
            if plotit is True:
                self.fig1, self.ax1 = plt.subplots()
                self.plot_magnetic_field()
                self.ax1.plot([-0.2, self.channel_length * 1e3 * 1.2], [0, 0], 'yellow')
                self.ax1.plot([-0.2, self.channel_length * 1e3 * 1.2],
                              [self.channel_height * 1000, self.channel_height * 1000], 'black')

            for init_particle_pos in initial_particle_position:
                initial_particle_pos = np.asarray(init_particle_pos)
                initial_state = np.hstack((initial_particle_pos, [0, 0, 0]))
                self.solver(initial_state=initial_state)

            if plotit is True:
                for traj_solutions in self.particle_path:
                    if traj_solutions != 'error':
                        if traj_solutions.t_events[0].size > 0:
                            color_path = 'deepskyblue'
                        else:
                            color_path = 'tomato'
                        particle_path = traj_solutions.y
                        self.ax1.plot(particle_path[0, :] * 1e3, particle_path[2, :] * 1e3, color=color_path,
                                      linewidth=2.0)

    def calculate_h_field(self, particle_location):
        # Evaluate H field on particle location
        h_field = magpy.getH(self.magnet_collection, particle_location)
        # h_field is given in kA/m, thus we transform it to the SI unit A/m
        h_field = h_field * 1000
        return h_field

    def plot_magnetic_field(self):

        # create grid
        xmin = -0.2
        xmax = self.channel_length * 1e3 * 1.2
        ymin = -2.0
        ymax = self.channel_height * 1e3 + 1.3

        ts = np.linspace(xmin, xmax, 20)
        zs = np.linspace(ymin, ymax, 20)
        self.ax1.set_xlim([xmin, xmax])
        self.ax1.set_ylim([ymin, ymax])
        grid = np.array([[(x, 0, z) for x in ts] for z in zs])

        # compute and plot field of coil2
        b_field = magpy.getH(self.magnet_collection, grid)
        b_amplitude = np.linalg.norm(b_field, axis=2)
        b_amplitude /= np.amax(b_amplitude)

        cp = self.ax1.contourf(
            grid[:, :, 0], grid[:, :, 2], b_amplitude,
            levels=100,
            cmap='coolwarm',
        )
        self.ax1.streamplot(
            grid[:, :, 0], grid[:, :, 2], b_field[:, :, 0], b_field[:, :, 2],
            density=1,
            color='black',
        )

        # figure styling
        self.ax1.set(
            # title='Magnetic field of cylinder',
            xlabel='x-position [mm]',
            ylabel='z-position [mm]',
            aspect=1,
        )

        # plt.colorbar(cp, ax=self.ax1, label='% of max field strength')

    def plot_particle_path(self):
        for sol in self.particle_path:
            particle_path = sol.y * 1000
            particle_path = particle_path[[0, 2], :]
            # self.ax1.plot(particle_path[0, :] * 1000, particle_path[1, :] * 1000, 'red')
            self.ax1.plot([-0.2, self.channel_length * 1e3 * 1.2], [0, 0], 'yellow')
            self.ax1.plot([-0.2, self.channel_length * 1e3 * 1.2],
                          [self.channel_height * 1000, self.channel_height * 1000], 'black')

            particle_path = np.transpose(particle_path)
            xy = particle_path.reshape(-1, 1, 2)
            segments = np.hstack([xy[:-1], xy[1:]])
            if sol.t_events[0].size > 0:
                color_path = 'deepskyblue'
            else:
                color_path = 'tomato'
            coll = LineCollection(segments, colors=color_path)  # , array=speed, cmap=plt.cm.autumn)

            self.ax1.add_collection(coll)

    def movement_model(self, time, solution):
        x = solution[0:3]
        velocity_flow = self.flow_model(x)
        dx_dt = (self.get_gravitational_force() + self.get_magnetic_force(x)) / (
                3 * self.mu * math.pi * self.hydraulic_drag_diameter) + velocity_flow
        v = dx_dt
        solution = np.hstack((dx_dt, v))
        return solution

    def get_magnetic_force(self, particle_position):
        h_field = self.calculate_h_field(
            (particle_position[0] * 1000, particle_position[1] * 1000, particle_position[2] * 1000))

        eps = 0.0000001
        pos1 = (particle_position[0] * 1000 + eps * h_field[0], particle_position[1] * 1000 + eps * h_field[1],
                particle_position[2] * 1000 + eps * h_field[2])
        pos2 = (particle_position[0] * 1000 - eps * h_field[0], particle_position[1] * 1000 - eps * h_field[1],
                particle_position[2] * 1000 - eps * h_field[2])

        dot_prod = (self.calculate_h_field(pos1) - self.calculate_h_field(pos2)) / (2 * eps) * 1e3
        # 1e3 for transforming mm^-1 in m^-1
        K = 3 * self.magnetic_susceptibility / (3 + self.magnetic_susceptibility)
        maximum_magnetization = 86 * self.particle_density * 3 / 2  # [Am^2/kg*kg/m^3], mass magnetization saturation was set to 86 [Am^2/kg], according to iron oxide
        if np.linalg.norm(h_field) < K * maximum_magnetization:
            relative_perm = K
        else:
            relative_perm = maximum_magnetization / np.linalg.norm(h_field)  # [-]
        magnetic_force = self.mu_0 * self.magnetic_volume * relative_perm * dot_prod
        return magnetic_force

    def get_gravitational_force(self):
        gravitational_force = self.gravitational_acceleration * self.magnetic_volume * (
                self.particle_density - self.fluid_density)

        return gravitational_force

    def flow_model(self, particle_position):
        if self.parabolic_flow:
            velocity_flow = np.array([6 * self.flow_rate / (self.channel_width * self.channel_height ** 3) * (
                    self.channel_height * particle_position[2] - particle_position[2] ** 2), 0, 0])
        else:
            velocity_flow = np.array([self.flow_rate / (self.channel_width * self.channel_height), 0, 0])
        return velocity_flow

    def solver(self, initial_state):
        Event.hit_bottom.terminal = True
        max_distance = self.channel_length * 1.1
        travel_event = lambda t, x: Event.travel_too_far(t, x, max_distance)
        travel_event.terminal = True
        hit_upper_event = lambda t, x: Event.hit_upper_wall(t, x, self.channel_height)
        hit_upper_event.terminal = True
        Event.travel_too_far.terminal = True

        self.max_step = self.channel_length / self.v_max[0] / self.n_step

        try:
            sol = solve_ivp(fun=self.movement_model, t_span=[0, 1000], y0=initial_state, max_step=self.max_step,
                            events=(Event.hit_bottom, travel_event, hit_upper_event))
        except:
            sol = 'error'
            print('solver error')
            pass
        self.particle_path.append(sol)
        return sol

    def determine_h_max(self, iterations, plotit):
        x0 = 0.0
        y0 = 0.0
        hit_bottom = False
        termination = False
        step_decay = 0.99
        step_size = self.channel_height / 15
        calculation_position = self.channel_height * 0.99
        termination_criterion = 1
        i = 0
        j = 0
        initial_state = [x0, y0, calculation_position, 0, 0, 0]

        sol = self.solver(initial_state=initial_state)
        if sol == 'error':
            return calculation_position
        else:
            if plotit is True:
                self.plot_particle_path()
            for n in range(3):
                if sol.t_events[n].size > 0:
                    termination_criterion = n + 1
            if termination_criterion == 1:
                calculation_position = self.channel_height
                termination = True
            step_size = step_size * step_decay

        while i < iterations and j < 12:
            if termination is True or sol == 'error':
                break
            if hit_bottom == True:
                step_decay = 0.6
            if sol.message == 'The solver successfully reached the end of the integration interval.':
                calculation_position = calculation_position - step_size
            elif termination_criterion == 1:
                calculation_position = calculation_position + step_size
                if calculation_position > self.channel_height:
                    calculation_position = self.channel_height
                hit_bottom = True
            elif termination_criterion == 2:
                calculation_position = calculation_position - step_size
            elif termination_criterion == 3:
                calculation_position = calculation_position - step_size
            else:
                raise 'Value Error: No termination criterion reached'

            if sol == 'error':
                return calculation_position
            else:
                initial_state = [x0, y0, calculation_position, 0, 0, 0]
                # print(calculation_position)
                sol = self.solver(initial_state=initial_state)
                if sol == 'error':
                    return calculation_position
                if plotit is True:
                    self.plot_particle_path()
                for n in range(3):
                    if sol.t_events[n].size > 0:
                        termination_criterion = n + 1
                step_size = step_size * step_decay

                i += 1
                if hit_bottom:
                    j += 1
        return calculation_position

    def captured_bacteria(self):
        # Separation efficiency
        if self.parabolic_flow:
            separation_efficiency = 6 / self.channel_height ** 3 * (
                        self.channel_height * self.h_max ** 2 / 2 - self.h_max ** 3 / 3)
        else:
            separation_efficiency = self.h_max / self.channel_height
        return separation_efficiency

    def force_monitor(self):

        for i in range(len(self.particle_path)):
            particle_path = self.particle_path[i].y[0:3, :]

            for j in range(len(particle_path[1, :])):
                self.F_g.append(self.get_gravitational_force())
                self.F_mag.append(self.get_magnetic_force(particle_path[:, j]))

    def force_plot(self):
        activate_force_z = False
        if activate_force_z:
            fig, ax = plt.subplots()
            plt.title('Forces in z')
            ax.set_xlabel('z [m]')
            ax.set_ylabel('F [N]')
            plt.yscale('log')
            F_mag_ar = np.asarray(self.F_mag)
            F_g = min(self.F_g[0])
            prev_count = 0
            ax.hlines(y=-F_g, xmin=0.0, xmax=0.003, linewidth=2, colors='r', label='gravity force')
            for i in range(len(self.particle_path)):
                particle_path = np.asarray(self.particle_path[i].y[0:3, :])
                path_len = len(particle_path[1, :])
                F_mag = F_mag_ar[prev_count:prev_count + path_len, 2]
                label_f_mag = 'Force path' + str(i)
                ax.plot(particle_path[2, :], -F_mag, label=label_f_mag)
                prev_count = prev_count + path_len
            plt.legend()

        activate_grav_force = False
        if activate_grav_force:
            fig, ax = plt.subplots()
            plt.title('grav Forces')
            ax.set_xlabel('z [m]')
            ax.set_ylabel('F [N]')
            plt.yscale('log')
            F_mag_ar = np.asarray(self.F_mag)
            prev_count = 0
            for i in range(len(self.particle_path)):
                contained = False
                particle_path = np.asarray(self.particle_path[i].y[0:3, :])
                path_len = len(particle_path[1, :])
                list_paths = [10, 12]
                if i in list_paths:
                    F_mag_x = F_mag_ar[prev_count:prev_count + path_len, 0]
                    F_mag_y = F_mag_ar[prev_count:prev_count + path_len, 1]
                    F_mag_z = F_mag_ar[prev_count:prev_count + path_len, 2]
                    label_f_x = 'Force path x' + str(i)
                    label_f_y = 'Force path y' + str(i)
                    label_f_z = 'Force path z' + str(i)
                    ax.plot(particle_path[2, :], abs(F_mag_x), label=label_f_x)
                    ax.plot(particle_path[2, :], abs(F_mag_y), label=label_f_y)
                    ax.plot(particle_path[2, :], abs(F_mag_z), label=label_f_z)
                prev_count = prev_count + path_len
            plt.legend()
