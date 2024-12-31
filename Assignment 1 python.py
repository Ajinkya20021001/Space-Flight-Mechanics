import numpy as np
import tkinter as tk
from tkinter import ttk

def calculate_position_velocity(r0, v0, t):
    """
    Calculate position and velocity using f and g functions
    """
    # Constants
    mu = 398600.4418  # Earth's gravitational parameter (km³/s²)
    
    # Convert time to seconds and inputs to numpy arrays
    t = t * 3600  # convert hours to seconds
    r0 = np.array(r0)
    v0 = np.array(v0)
    
    # Calculate initial magnitudes
    r0_mag = np.linalg.norm(r0)
    v0_mag = np.linalg.norm(v0)
    
    # Calculate initial radial velocity
    v0_rad = np.dot(r0, v0) / r0_mag
    
    # Calculate semi-major axis
    alpha = 2/r0_mag - v0_mag**2/mu
    a = 1/alpha  # semi-major axis
    
    # Calculate initial eccentric anomaly
    chi = np.sqrt(mu) * abs(alpha) * t
    
    # Universal variable iteration
    z = alpha * chi**2
    
    # Calculate f and g functions
    f = 1 - (chi**2/r0_mag) * (1 - z/2)
    g = t - (chi**3/(6*np.sqrt(mu))) * (1 - z/3)
    
    # Calculate position
    r = f * r0 + g * v0
    
    # Calculate derivatives of f and g
    fdot = (np.sqrt(mu)/(r0_mag * np.linalg.norm(r))) * chi * (z - 1)
    gdot = 1 - (chi**2/np.linalg.norm(r)) * (1 - z/2)
    
    # Calculate velocity
    v = fdot * r0 + gdot * v0
    
    return r, v

class OrbitCalculatorUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Orbit Calculator")
        
        # Create and configure main frame
        main_frame = ttk.Frame(root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Position inputs
        ttk.Label(main_frame, text="Initial Position (km)").grid(row=0, column=0, columnspan=3)
        ttk.Label(main_frame, text="x:").grid(row=1, column=0)
        ttk.Label(main_frame, text="y:").grid(row=1, column=1)
        ttk.Label(main_frame, text="z:").grid(row=1, column=2)
        
        self.rx = ttk.Entry(main_frame, width=15)
        self.ry = ttk.Entry(main_frame, width=15)
        self.rz = ttk.Entry(main_frame, width=15)
        self.rx.grid(row=2, column=0, padx=5)
        self.ry.grid(row=2, column=1, padx=5)
        self.rz.grid(row=2, column=2, padx=5)
        
        # Velocity inputs
        ttk.Label(main_frame, text="Initial Velocity (km/s)").grid(row=3, column=0, columnspan=3, pady=(10,0))
        ttk.Label(main_frame, text="vx:").grid(row=4, column=0)
        ttk.Label(main_frame, text="vy:").grid(row=4, column=1)
        ttk.Label(main_frame, text="vz:").grid(row=4, column=2)
        
        self.vx = ttk.Entry(main_frame, width=15)
        self.vy = ttk.Entry(main_frame, width=15)
        self.vz = ttk.Entry(main_frame, width=15)
        self.vx.grid(row=5, column=0, padx=5)
        self.vy.grid(row=5, column=1, padx=5)
        self.vz.grid(row=5, column=2, padx=5)
        
        # Time input
        ttk.Label(main_frame, text="Time (hours):").grid(row=6, column=0, pady=(10,0))
        self.time = ttk.Entry(main_frame, width=15)
        self.time.grid(row=6, column=1, pady=(10,0))
        
        # Calculate button
        ttk.Button(main_frame, text="Calculate", command=self.calculate).grid(row=7, column=0, columnspan=3, pady=10)
        
        # Results display
        self.result_text = tk.Text(main_frame, height=4, width=50)
        self.result_text.grid(row=8, column=0, columnspan=3)
        
        # Set default values
        self.rx.insert(0, "20000")
        self.ry.insert(0, "-105000")
        self.rz.insert(0, "-19000")
        self.vx.insert(0, "0.9")
        self.vy.insert(0, "-3.4")
        self.vz.insert(0, "-1.5")
        self.time.insert(0, "2")
    
    def calculate(self):
        try:
            r0 = [float(self.rx.get()), float(self.ry.get()), float(self.rz.get())]
            v0 = [float(self.vx.get()), float(self.vy.get()), float(self.vz.get())]
            t = float(self.time.get())
            
            r_final, v_final = calculate_position_velocity(r0, v0, t)
            
            result = f"Final position (km):\n"
            result += f"r = {r_final[0]:.2f}i + {r_final[1]:.2f}j + {r_final[2]:.2f}k\n"
            result += f"Final velocity (km/s):\n"
            result += f"v = {v_final[0]:.4f}i + {v_final[1]:.4f}j + {v_final[2]:.4f}k"
            
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, result)
            
        except ValueError as e:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, "Please enter valid numbers")
        except Exception as e:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, f"Error in calculation: {str(e)}")

def main():
    root = tk.Tk()
    app = OrbitCalculatorUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()