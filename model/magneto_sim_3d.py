# %% Import
import magpylib as magpy
import numpy as np
import matplotlib.pyplot as plt
import math
from math import pi, sin, cos,tanh, sinh, cosh
from scipy.integrate import solve_ivp
import Event
from matplotlib.collections import LineCollection
from matplotlib import cm

# %% Def
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
            particle_density=5240,
            bacterium_volume=4.71e-18,
            bacterium_density=1100,
            magnet=None,
            mu_0=4e-7 * pi,
            gravitational_acceleration=np.array([0, 0, -9.81]),
            fluid_density=997,
            mu=0.001,  # [Pas] viscosity fluid
            number_of_steps=12,
            evaluation_planes=6,
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
        self.bacterium_volume = bacterium_volume
        self.bacterium_density = bacterium_density
        self.h_max = []
        self.F_mag = []
        self.F_g = []
        self.speed = []
        self.h_grad_h = []
        self.sol = None
        self.flow_rate = flow_rate
        self.poly_h_sep = None
        self.magnetic_volume_dist = magnetic_volume_dist
        self.magnetic_volume_frac = magnetic_volume_frac
        self.dp_dx = None
        self.frac_captured_bacteria = None
        self.total_frac_captured_bacteria = None
        self.max_step = None
        self.v_max = None
        self.n_step = number_of_steps
        self.evaluation_planes = evaluation_planes

    def main_run(self, calc_mod=None, iterations=25, plotit=False, msg=True, plot_3d_traj=True):

        for n in range(len(self.magnet)):
            self.magnet_collection.add(self.magnet[n])

        self.flow_model([0, 0, 0])
        if plotit:
            self.plot_velocity_field()
        print('pressure gradient: ' + str(round(self.dp_dx, 3)) + " Pa/m")
        self.hydraulic_drag_diameter = (6 / pi * (self.magnetic_volume + self.bacterium_volume)) ** (1 / 3)
        self.v_max = np.round(self.flow_model([0, 0, self.channel_height / 2]), 4)
        print("max flow velocity: " + str(self.v_max) + " m/s")

        # Plots the separation height & fraction of captured bacteria in 3D over channel width (y-Dir);
        # plotit-True: Particle traces are displayed; V_mag is one constant value
        if calc_mod == 0:
            y_var = np.linspace(- self.channel_width / 2 * 0.95, 0.0, self.evaluation_planes)
            self.hydraulic_drag_diameter = (6 / pi * (self.magnetic_volume + self.bacterium_volume)) ** (1 / 3)
            if plot_3d_traj:
                fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                ax.set_xlabel('x [m]')
                ax.set_ylabel('y [m]')
                ax.set_zlabel('z [m]')
                ax.set_box_aspect((self.channel_length * 1.1, self.channel_width, self.channel_height))

            for ind_l, y_pos in enumerate(y_var):
                if plotit:
                    self.fig1, self.ax1 = plt.subplots()
                    self.plot_magnetic_field(y_pos=y_pos)
                self.particle_path = []
                h_sep = self.determine_h_sep(iterations, plotit, y_pos)

                if plot_3d_traj:
                    for traj_solutions in self.particle_path:
                        if traj_solutions != 'error':
                            if traj_solutions.t_events[0].size > 0:
                                color_path = 'deepskyblue'
                            else:
                                color_path = 'tomato'
                            particle_path = traj_solutions.y
                            ax.plot(particle_path[0, :], particle_path[1, :], particle_path[2, :], color=color_path)

                if h_sep > self.channel_height:
                    self.h_max.append(self.channel_height)
                elif h_sep < 0:
                    self.h_max.append(0)
                else:
                    self.h_max.append(h_sep)
                if msg:
                    print('Determined h_sep: {}'.format(round(self.h_max[ind_l], 6)))
                if plotit:
                    plt.tight_layout()
                    plt.show()
                    fig_name = 'img/h_sep_particle_trac_at_y=' + str(round(y_pos, 5))
                    plt.savefig(fig_name + '.pdf', format='pdf')
                    plt.savefig(fig_name + '.svg', format='svg')

            # Plotting h_sep over the channel width
            h_max_rev = self.h_max[::-1]
            h_max_plt = np.concatenate((self.h_max, h_max_rev[1:]), axis=None)
            y_var_rev = -y_var[::-1]
            y_var_plt = np.concatenate((y_var, y_var_rev[1:]), axis=None)
            poly_results = np.polyfit(y_var_plt, h_max_plt, 4)
            self.poly_h_sep = np.poly1d(poly_results)
            plot_h_sep = plot_3d_traj
            if plot_h_sep:
                fig, ax = plt.subplots()
                plt.title('Line of Separation')
                ax.set_xlabel('y [m]')
                ax.set_ylabel('h_sep [m]')
                ax.scatter(y_var_plt, h_max_plt, color='k', marker='x')
                y_var_high_res = np.linspace(min(y_var_plt), max(y_var_plt), 100)
                ax.plot(y_var_high_res, self.poly_h_sep(y_var_high_res), color='k')
                ax.set_xlim([-self.channel_width / 2, self.channel_width / 2])
                ax.set_ylim([0, self.channel_height])
                ax.set(aspect=1,)
                fig_name = 'img/h_sep_over_y'
                plt.savefig(fig_name + '.pdf', format='pdf')
                plt.savefig(fig_name + '.svg', format='svg')

            # calculating the fraction of captured bacteria
            self.frac_captured_bacteria = self.captured_bacteria()
            self.total_frac_captured_bacteria = self.frac_captured_bacteria
            if msg:
                print('The total Separation Efficiency is: SE = ' + str(
                    round(self.frac_captured_bacteria * 100, 4)) + '%')

        elif calc_mod == 1:
            # calculates the (total) fraction of captured bacteria and separation height; V_mag is different for every loop;
            self.frac_captured_bacteria = []
            y_var = np.linspace(- self.channel_width / 2 * 0.95, 0.0, self.evaluation_planes)
            for ind_vmag, V_magnetic in enumerate(self.magnetic_volume_dist):
                self.hydraulic_drag_diameter = (6 / pi * (V_magnetic + self.bacterium_volume)) ** (1 / 3)
                self.h_max = []
                self.magnetic_volume = V_magnetic
                for ind_l, y_pos in enumerate(y_var):
                    if plotit is True:
                        self.fig1, self.ax1 = plt.subplots()
                        self.plot_magnetic_field(y_pos=y_pos)
                    self.particle_path = []
                    h_sep = self.determine_h_sep(iterations, plotit, y_pos)
                    if h_sep > self.channel_height:
                        self.h_max.append(self.channel_height)
                    elif h_sep < 0:
                        self.h_max.append(0)
                    else:
                        self.h_max.append(h_sep)
                    if plotit is True:
                        plt.tight_layout()
                        plt.show()
                        fig_name = 'img/h_sep_particle_trac_at_y=' + str(round(y_pos, 5))
                        plt.savefig(fig_name + '.pdf', format='pdf')
                        plt.savefig(fig_name + '.svg', format='svg')
                        # plt.close()

                # Plotting h_sep over the channel width
                h_max_rev = self.h_max[::-1]
                h_max_plt = np.concatenate((self.h_max, h_max_rev[1:]), axis=None)
                y_var_rev = -y_var[::-1]
                y_var_plt = np.concatenate((y_var, y_var_rev[1:]), axis=None)
                poly_results = np.polyfit(y_var_plt, h_max_plt, 4)
                self.poly_h_sep = np.poly1d(poly_results)

                plot_h_sep = False  # 'True' for separation lines
                if plot_h_sep:
                    fig, ax = plt.subplots()
                    plt.title('Line of Separation')
                    ax.set_xlabel('y [m]')
                    ax.set_ylabel('h_sep [m]')
                    y_var_high_res = np.linspace(min(y_var_plt), max(y_var_plt), 100)
                    ax.scatter(y_var_plt, h_max_plt, color='k', marker='x')
                    ax.plot(y_var_high_res, self.poly_h_sep(y_var_high_res), color='k')
                    ax.set_xlim([-self.channel_width / 2, self.channel_width / 2])
                    ax.set_ylim([0, self.channel_height])
                    fig_name = 'img/h_sep_over_y'
                    plt.savefig(fig_name + '.pdf', format='pdf')
                    plt.savefig(fig_name + '.svg', format='svg')

                # calculating the fraction of captured bacteria
                self.frac_captured_bacteria.append(self.captured_bacteria())
                print('SE = ' + str(round(self.captured_bacteria() * 100, 4)) + '%')
            self.total_frac_captured_bacteria = sum(np.multiply(self.frac_captured_bacteria, self.magnetic_volume_frac))
            if msg:
                print('The total Separation Efficiency is: SE = ' + str(
                    round(self.total_frac_captured_bacteria * 100, 4)) + '%')

        elif calc_mod == 2:
            # force and particle tracking depending on release position

            # Release Position Grid
            y_n = np.linspace(-self.channel_width / 2 * 0.8, self.channel_width / 2 * 0.8, 15)
            z_n = np.linspace(self.channel_height * 0.02, self.channel_height * 0.98, 20)
            ipp_grid = []
            for i, i_n in enumerate(y_n):
                for j, j_n in enumerate(z_n):
                    ipp_grid.append((0, round(i_n, 7), round(j_n, 7)))
            initial_particle_position = ipp_grid

            self.hydraulic_drag_diameter = (6 / pi * (self.magnetic_volume + self.bacterium_volume)) ** (1 / 3)

            if plot_3d_traj:
                fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                ax.set_xlabel('x [m]')
                ax.set_ylabel('y [m]')
                ax.set_zlabel('z [m]')
                ax.set_box_aspect((self.channel_length * 1.1, self.channel_width, self.channel_height))
                ax.set_xlim([0, self.channel_length * 1.1])
                ax.set_ylim([0, self.channel_width / 2])
                ax.set_zlim([0, self.channel_height])
                ax.set_xticks([0, 0.015 - 0.0035, 0.015])  # ([0, 0.015-0.0035, 0.015])
                ax.set_yticks([self.channel_width / -2,
                               self.channel_width / 2])  # ([self.channel_width/-2, self.channel_width/2])
                ax.set_zticks([0, self.channel_height])  # ([0, self.channel_height])

                for init_particle_pos in initial_particle_position:
                    initial_particle_pos = np.asarray(init_particle_pos)
                    initial_state = np.hstack((initial_particle_pos, [0, 0, 0]))
                    self.solver(initial_state=initial_state)

                if plot_3d_traj:
                    for traj_solutions in self.particle_path:
                        if traj_solutions != 'error':
                            if traj_solutions.t_events[0].size > 0:
                                color_path = 'deepskyblue'
                            else:
                                color_path = 'tomato'
                            particle_path = traj_solutions.y
                            ax.plot(particle_path[0, :], particle_path[1, :], particle_path[2, :], color=color_path)

            if plot_3d_traj:
                fig_name = 'traj/Praticle_trajectories_3D'
                plt.savefig(fig_name + '.svg', format='svg')

    def flow_model(self, particle_position):
        n = 3
        K_pres = 1
        for num1 in np.arange(1, n + 1, 2):
            K_pres -= 1 / num1 ** 5 * 192 / pi ** 5 * self.channel_height / self.channel_width * tanh(num1 * pi * self.channel_width / (2 * self.channel_height))
        G_dpdx = self.flow_rate * 12 * self.mu / (K_pres * self.channel_height ** 3 * self.channel_width)
        self.dp_dx = -G_dpdx
        K_vel = 0
        for num2 in np.arange(1, n + 1, 2):
            K_vel += 1 / (num2) ** 3 * (1 - cosh(num2 * pi * particle_position[1] / self.channel_height) / cosh(num2 * pi * self.channel_width / (2 * self.channel_height))) * sin(
                num2 * pi * particle_position[2] / self.channel_height)
        u_yz = 4 * self.channel_height ** 2 * G_dpdx / (pi ** 3 * self.mu) * K_vel
        if u_yz < 0:
            velocity_flow = np.array([0, 0, 0])
        else:
            velocity_flow = np.array([u_yz, 0, 0])
        return velocity_flow

    def volume_flux(self):
        n=3
        coeff_flux = 0
        for num2 in np.arange(1, n + 1, 1):
            beta2 = (2 * num2 - 1) * math.pi / self.channel_height
            coeff_flux += 1 / (2 * num2 - 1) ** 5 * (math.cosh(beta2 * self.channel_width) - 1) / math.sinh(
                beta2 * self.channel_width)
        total_volume_flux = self.dp_dx/self.mu*(-self.channel_height ** 3 * self.channel_width / 12 + 16 * self.channel_height ** 4 / pi ** 5 * coeff_flux)
        return total_volume_flux

    def plot_velocity_field(self):
        step_size = 0.00005
        Y_pos = np.arange(-self.channel_width / 2, self.channel_width / 2, step_size)
        Z_pos = np.arange(0, self.channel_height, step_size)
        Y_pos, Z_pos = np.meshgrid(Y_pos, Z_pos)
        vel_fluid = np.zeros([len(Z_pos[:, 0]), len(Y_pos[0, :])])
        for i_y in np.arange(0, len(Y_pos[0, :]), 1):
            for i_z in np.arange(0, len(Z_pos[:, 0]), 1):
                particle_pos = np.array([0, Y_pos[i_z, i_y], Z_pos[i_z, i_y]])
                vel_fluid_attr = self.flow_model(particle_pos)
                u_vel = vel_fluid_attr[0]
                if u_vel < 0:
                    vel_fluid[i_z, i_y] = 0
                else:
                    vel_fluid[i_z, i_y] = u_vel
        # Plot the surface
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.plot_surface(Y_pos, Z_pos, vel_fluid, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        # Customize the axis.
        ax.set_xlabel('y [m]')
        ax.set_ylabel('z [m]')
        ax.set_zlabel('u [m/s]')

        # plt.savefig('velocity_field_channelflow.svg', format='svg')
        # plt.close()
        # plt.show()

    def calculate_h_field(self, particle_location):
        # Evaluate H field on particle location
        h_field = magpy.getH(self.magnet_collection, particle_location)
        # h_field is given in kA/m, thus we transform it to the SI unit A/m
        h_field = h_field * 1000
        return h_field

    def plot_magnetic_field(self, y_pos):
        # create grid
        xmin = -0.2
        xmax = self.channel_length * 1e3 * 1.2
        ymin = -2.0
        ymax = self.channel_height * 1e3 + 1.3

        ts = np.linspace(xmin, xmax, 20)
        zs = np.linspace(ymin, ymax, 20)
        self.ax1.set_xlim([xmin, xmax])
        self.ax1.set_ylim([ymin, ymax])
        grid = np.array([[(x, y_pos, z) for x in ts] for z in zs])

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
            density=2,
            color='black',
        )

        # figure styling
        self.ax1.set(
            title='Magnetic field of cylinder at y=' + str(round(y_pos, 5)),
            xlabel='x-position [mm]',
            ylabel='z-position [mm]',
            aspect=1,
        )

        plt.colorbar(cp, ax=self.ax1, label='% of max field strength')

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
            coll = LineCollection(segments, colors=color_path)

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
        maximum_magnetization = 86 * self.particle_density  # [Am^2/kg*kg/m^3], mass magnetization saturation was set to 86 [Am^2/kg], according to iron oxide

        if np.linalg.norm(h_field) < K*maximum_magnetization:
            relative_perm = K
        else:
            relative_perm = maximum_magnetization / np.linalg.norm(h_field)  # [-]
        magnetic_force = self.mu_0 * self.magnetic_volume * relative_perm * dot_prod
        return magnetic_force

    def get_gravitational_force(self):
        gravitational_force = self.gravitational_acceleration * self.magnetic_volume * (
                self.particle_density - self.fluid_density)  # grav force and buoyancy
        return gravitational_force

    def solver(self, initial_state):
        Event.hit_bottom.terminal = True
        max_distance = self.channel_length * 1.1
        travel_event = lambda t, x: Event.travel_too_far(t, x, max_distance)
        travel_event.terminal = True
        hit_upper_event = lambda t, x: Event.hit_upper_wall(t, x, self.channel_height)
        hit_upper_event.terminal = True
        Event.travel_too_far.terminal = True

        # For higher velocities, the solver needs smaller step sizes. As very small step sizes lead to high simulation durations.
        # Thus, the step size can be selected according to the velocity relative to the channel length.
        # This calculation can be moved to the main()-function.
        max_velocity = self.flow_model([0, initial_state[1], self.channel_height / 2])[0]
        self.max_step = self.channel_length / max_velocity / self.n_step
        try:
            sol = solve_ivp(fun=self.movement_model, t_span=[0, 1000], y0=initial_state, max_step=self.max_step,
                            events=(Event.hit_bottom, travel_event, hit_upper_event))
        except:
            sol = 'error'
            print('solver error')
            pass

        self.particle_path.append(sol)
        return sol

    def determine_h_sep(self, iterations, plotit, y_pos):
        x0 = 0
        y0 = y_pos
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
            # for i in range(iterations):
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
        # Separation Efficiency calculation
        y_var = np.linspace(0, self.channel_width, 40)

        # Total Volume Flux
        coeff_vol = 0
        for num in np.arange(1, 5, 1):
            beta = (2 * num - 1) * pi / self.channel_height
            coeff_vol += 1 / (2 * num - 1) ** 5 * (cosh(beta * self.channel_width) - 1) / sinh(
                beta * self.channel_width)
        V_dot_s = -self.channel_height ** 3 * self.channel_width / 12 + 16 * self.channel_height ** 4 / pi ** 5 * coeff_vol

        # Volume Flux throught separation height part
        V_help = []
        for y_val in y_var:
            V_help.append(-1 / 4 * self.channel_height * self.poly_h_sep(
                y_val - self.channel_width / 2) ** 2 + 1 / 6 * self.poly_h_sep(
                y_val - self.channel_width / 2) ** 3 + 4 * self.channel_height ** 2 / pi ** 3 * self.sum_coeff(y_val,
                                                                                                               3))
        V_dot_h_sep = np.trapz(V_help, y_var)
        # fraction of captured bacteria
        c_b = V_dot_h_sep / V_dot_s
        return c_b

    def sum_coeff(self, y_val, n):
        sum_K = 0
        for num in np.arange(1, n + 1, 1):
            beta = (2 * num - 1) * pi / self.channel_height
            sum_K += 1 / beta / (2 * num - 1) ** 3 * (
                    sinh(beta * y_val) + sinh(beta * (self.channel_width - y_val))) / sinh(
                beta * self.channel_width) * (1 - cos(beta * self.poly_h_sep(y_val - self.channel_width / 2)))
        return sum_K