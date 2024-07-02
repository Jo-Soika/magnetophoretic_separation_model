import tkinter as tk
from tkinter import ttk, messagebox

# Import your simulation function from the simulation file
def run_simulation():
    try:
        params = {
            "h_ch": float(h_ch_entry.get()),
            "b_ch": float(b_ch_entry.get()),
            "l_ch": float(l_ch_entry.get()),
            "V_dot": float(V_dot_entry.get()),
            "magnet_shape": magnet_shape_var.get(),
            "polarization_val": int(polarization_val_entry.get()),
            "length_magnet": float(length_magnet_entry.get()),
            "d_mag": float(d_mag_entry.get()),
            "distance_mag_to_fluid": float(distance_mag_to_fluid_entry.get()),
            "activate_3d_dimensions": activate_3d_dimensions_var.get(),
            "activate_plots": activate_plots_var.get(),
            "activate_parabolic_flow": activate_parabolic_flow_var.get(),
            "calculation_mode": int(calculation_mode_var.get())
        }
        # Placeholder for the actual simulation function call
        # result = run_millifluidic_simulation(params)
        # messagebox.showinfo("Simulation Result", result)
        print(params)  # Just to show the parameters collected, remove this in production
    except ValueError as e:
        messagebox.showerror("Input Error", "Please check your inputs. All fields must be filled with the correct type.")
def setup_gui():
    root = tk.Tk()
    root.title("Magnetophoretic Separation Model GUI")

    # Default settings
    defaults = {
        "h_ch": 0.0035,
        "b_ch": 0.0035,
        "l_ch": 0.015,
        "V_dot": 1e-07,
        "magnet_shape": "cylinder",
        "polarization_val": 1500,
        "length_magnet": 0.01,
        "d_mag": 0.0035,
        "distance_mag_to_fluid": 0.0001,
        "activate_3d_dimensions": True,
        "activate_plots": True,
        "activate_parabolic_flow": True,
        "calculation_mode": 0
    }

    # Channel dimensions
    ttk.Label(root, text="Channel Height (m):").grid(row=0, column=0)
    h_ch_entry = ttk.Entry(root)
    h_ch_entry.grid(row=0, column=1)
    h_ch_entry.insert(0, defaults["h_ch"])

    ttk.Label(root, text="Channel Width (m):").grid(row=1, column=0)
    b_ch_entry = ttk.Entry(root)
    b_ch_entry.grid(row=1, column=1)
    b_ch_entry.insert(0, defaults["b_ch"])

    ttk.Label(root, text="Channel Length (m):").grid(row=2, column=0)
    l_ch_entry = ttk.Entry(root)
    l_ch_entry.grid(row=2, column=1)
    l_ch_entry.insert(0, defaults["l_ch"])

    ttk.Label(root, text="Flow Rate (m^3/s):").grid(row=3, column=0)
    V_dot_entry = ttk.Entry(root)
    V_dot_entry.grid(row=3, column=1)
    V_dot_entry.insert(0, defaults["V_dot"])

    # Magnet configuration
    magnet_shape_var = tk.StringVar(value=defaults["magnet_shape"])
    ttk.Label(root, text="Magnet Shape:").grid(row=4, column=0)
    ttk.OptionMenu(root, magnet_shape_var, "cylinder", "cylinder", "cuboid").grid(row=4, column=1)

    ttk.Label(root, text="Polarization (mT):").grid(row=5, column=0)
    polarization_val_entry = ttk.Entry(root)
    polarization_val_entry.grid(row=5, column=1)
    polarization_val_entry.insert(0, defaults["polarization_val"])

    ttk.Label(root, text="Magnet Length (m):").grid(row=6, column=0)
    length_magnet_entry = ttk.Entry(root)
    length_magnet_entry.grid(row=6, column=1)
    length_magnet_entry.insert(0, defaults["length_magnet"])

    ttk.Label(root, text="Magnet Diameter/Side Length (m):").grid(row=7, column=0)
    d_mag_entry = ttk.Entry(root)
    d_mag_entry.grid(row=7, column=1)
    d_mag_entry.insert(0, defaults["d_mag"])

    ttk.Label(root, text="Distance Magnet to Fluid (m):").grid(row=8, column=0)
    distance_mag_to_fluid_entry = ttk.Entry(root)
    distance_mag_to_fluid_entry.grid(row=8, column=1)
    distance_mag_to_fluid_entry.insert(0, defaults["distance_mag_to_fluid"])

    # Additional settings
    activate_3d_dimensions_var = tk.BooleanVar(value=defaults["activate_3d_dimensions"])
    ttk.Checkbutton(root, text="Activate 3D Dimensions", variable=activate_3d_dimensions_var).grid(row=9, column=0)

    activate_plots_var = tk.BooleanVar(value=defaults["activate_plots"])
    ttk.Checkbutton(root, text="Activate Plots", variable=activate_plots_var).grid(row=10, column=0)

    activate_parabolic_flow_var = tk.BooleanVar(value=defaults["activate_parabolic_flow"])
    ttk.Checkbutton(root, text="Activate Parabolic Flow (2D only)", variable=activate_parabolic_flow_var).grid(row=11,
                                                                                                               column=0)

    calculation_mode_var = tk.StringVar(value=str(defaults["calculation_mode"]))
    ttk.Label(root, text="Calculation Mode:").grid(row=12, column=0)
    ttk.OptionMenu(root, calculation_mode_var, "0", "0", "1", "2").grid(row=12, column=1)

    run_button = ttk.Button(root, text="Run Simulation", command=run_simulation)
    run_button.grid(row=13, column=0, columnspan=2)

    return root

def main():
    root = setup_gui()
    root.mainloop()

if __name__ == "__main__":
    main()



