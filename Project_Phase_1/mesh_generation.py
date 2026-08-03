import gmsh
import numpy as np
import os

class BubbleArrayMesh:
    def __init__(self, L_electrode, H_electrolyte, theta_deg, alpha_wall, 
                 mesh_size=0.05, dimension=2):
        """
        Parameters:
        -----------
        L_electrode : float
            Electrode length (2D) or width (3D)
        H_electrolyte : float
            Electrolyte gap height
        theta_deg : float
            Contact angle in degrees [30, 120]
        alpha_wall : float
            Coverage ratio [0.1, 0.85]
        mesh_size : float
            Target mesh element size
        dimension : int
            2 or 3
        """
        self.L = L_electrode
        self.H = H_electrolyte
        self.theta = np.radians(theta_deg)
        self.alpha = alpha_wall
        self.mesh_size = mesh_size
        self.dim = dimension
        
        gmsh.initialize()
        gmsh.model.add("bubble_array")
        
    def compute_bubble_geometry(self):
        """Calculate bubble dimensions from coverage ratio"""
        if self.dim == 2:
            # Number of bubbles (approximate)
            N_bubbles = max(1, int(self.alpha * self.L / (2 * np.sin(self.theta))))
            # Adjust to fit exactly
            w_base = self.alpha * self.L / N_bubbles
            R_b = w_base / (2 * np.sin(self.theta))
            return N_bubbles, R_b, w_base
        else:
            # 3D: square array
            N_per_side = max(1, int(np.sqrt(self.alpha * self.L**2 / 
                                            (np.pi * np.sin(self.theta)**2))))
            N_bubbles = N_per_side**2
            A_base = self.alpha * self.L**2 / N_bubbles
            R_b = np.sqrt(A_base / (np.pi * np.sin(self.theta)**2))
            return N_bubbles, R_b, N_per_side
    
    def generate_2d_bubble(self, x_center, R_b):
        """Create a single 2D bubble (circular arc)"""
        # Bubble center (above electrode)
        y_center = R_b * np.cos(self.theta)
        
        # Create circle
        circle = gmsh.model.occ.addCircle(x_center, y_center, 0, R_b)
        
        # Create points at contact angles
        theta_start = np.pi/2 + self.theta
        theta_end = np.pi/2 - self.theta
        
        # Trim circle to create spherical cap
        # Use boolean operations to cut below electrode
        rectangle = gmsh.model.occ.addRectangle(x_center - 2*R_b, -R_b, 0, 
                                                4*R_b, y_center + R_b)
        
        bubble = gmsh.model.occ.cut([(2, rectangle)], [(2, circle)])
        return bubble[0]
    
    def generate_2d_mesh(self, output_file="mesh_2d.msh"):
        """Generate complete 2D mesh with bubble array"""
        N_bubbles, R_b, w_base = self.compute_bubble_geometry()
        
        # Create electrolyte domain
        electrolyte = gmsh.model.occ.addRectangle(0, 0, 0, self.L, self.H)
        
        # Create bubbles
        bubbles = []
        spacing = self.L / N_bubbles
        
        for i in range(N_bubbles):
            x_center = spacing * (i + 0.5)
            bubble_shapes = self.generate_2d_bubble(x_center, R_b)
            bubbles.extend(bubble_shapes)
        
        # Subtract bubbles from electrolyte
        if bubbles:
            fluid_domain = gmsh.model.occ.cut([(2, electrolyte)], bubbles)
        else:
            fluid_domain = [(2, electrolyte)]
        
        gmsh.model.occ.synchronize()
        
        # Physical groups for boundary conditions
        gmsh.model.addPhysicalGroup(1, [1], 1)  # Electrode (bottom)
        gmsh.model.setPhysicalName(1, 1, "electrode")
        
        gmsh.model.addPhysicalGroup(1, [2, 3, 4], 2)  # Other boundaries
        gmsh.model.setPhysicalName(1, 2, "bulk_boundary")
        
        # Bubble surfaces (zero-flux)
        bubble_surfaces = []
        for dim, tag in bubbles:
            if dim == 1:  # Lines
                bubble_surfaces.append(tag)
        
        if bubble_surfaces:
            gmsh.model.addPhysicalGroup(1, bubble_surfaces, 3)
            gmsh.model.setPhysicalName(1, 3, "bubble_interface")
        
        # Fluid domain
        gmsh.model.addPhysicalGroup(2, [fluid_domain[0][1]], 4)
        gmsh.model.setPhysicalName(2, 4, "electrolyte")
        
        # Mesh size field for refinement near bubbles
        field = gmsh.model.mesh.field
        field.add("Distance", 1)
        field.setNumbers(1, "EdgesList", bubble_surfaces)
        
        field.add("Threshold", 2)
        field.setNumber(2, "InField", 1)
        field.setNumber(2, "SizeMin", self.mesh_size * 0.3)
        field.setNumber(2, "SizeMax", self.mesh_size)
        field.setNumber(2, "DistMin", 0.01)
        field.setNumber(2, "DistMax", 0.1)
        
        field.add("Min", 3)
        field.setNumbers(3, "FieldsList", [2])
        field.setNumber(3, "NumThreads", 4)
        
        gmsh.model.mesh.field.setAsBackgroundMesh(3)
        
        # Generate mesh
        gmsh.model.mesh.generate(2)
        gmsh.write(output_file)
        gmsh.finalize()
        
        return {
            'N_bubbles': N_bubbles,
            'R_b': R_b,
            'w_base': w_base,
            'mesh_file': output_file
        }
    
    def generate_3d_mesh(self, output_file="mesh_3d.msh"):
        """Generate complete 3D mesh with bubble array"""
        N_bubbles, R_b, N_per_side = self.compute_bubble_geometry()
        
        # Create electrolyte domain (box)
        electrolyte = gmsh.model.occ.addBox(0, 0, 0, self.L, self.L, self.H)
        
        # Create bubbles (spherical caps)
        bubbles = []
        spacing = self.L / N_per_side
        
        for i in range(N_per_side):
            for j in range(N_per_side):
                x_center = spacing * (i + 0.5)
                y_center = spacing * (j + 0.5)
                z_center = R_b * np.cos(self.theta)
                
                # Create sphere
                sphere = gmsh.model.occ.addSphere(x_center, y_center, z_center, R_b)
                
                # Cut below electrode (z < 0)
                cut_box = gmsh.model.occ.addBox(x_center - 2*R_b, y_center - 2*R_b, 
                                               -2*R_b, 4*R_b, 4*R_b, 2*R_b)
                
                bubble = gmsh.model.occ.cut([(3, sphere)], [(3, cut_box)])
                bubbles.extend(bubble[0])
        
        # Subtract bubbles from electrolyte
        if bubbles:
            fluid_domain = gmsh.model.occ.cut([(3, electrolyte)], bubbles)
        else:
            fluid_domain = [(3, electrolyte)]
        
        gmsh.model.occ.synchronize()
        
        # Physical groups
        # Bottom surface (electrode)
        bottom_surfaces = []
        for dim, tag in gmsh.model.getEntities(2):
            bbox = gmsh.model.getBoundingBox(dim, tag)
            if abs(bbox[2]) < 1e-6:  # z_min = 0
                bottom_surfaces.append(tag)
        
        gmsh.model.addPhysicalGroup(2, bottom_surfaces, 1)
        gmsh.model.setPhysicalName(2, 1, "electrode")
        
        # Bubble surfaces
        bubble_surfaces = []
        for dim, tag in bubbles:
            if dim == 2:
                bubble_surfaces.append(tag)
        
        if bubble_surfaces:
            gmsh.model.addPhysicalGroup(2, bubble_surfaces, 2)
            gmsh.model.setPhysicalName(2, 2, "bubble_interface")
        
        # Fluid domain
        gmsh.model.addPhysicalGroup(3, [fluid_domain[0][1]], 3)
        gmsh.model.setPhysicalName(3, 3, "electrolyte")
        
        # Mesh refinement
        field = gmsh.model.mesh.field
        field.add("Distance", 1)
        field.setNumbers(1, "SurfacesList", bubble_surfaces)
        
        field.add("Threshold", 2)
        field.setNumber(2, "InField", 1)
        field.setNumber(2, "SizeMin", self.mesh_size * 0.3)
        field.setNumber(2, "SizeMax", self.mesh_size)
        field.setNumber(2, "DistMin", 0.01)
        field.setNumber(2, "DistMax", 0.1)
        
        field.add("Min", 3)
        field.setNumbers(3, "FieldsList", [2])
        gmsh.model.mesh.field.setAsBackgroundMesh(3)
        
        # Generate mesh
        gmsh.model.mesh.generate(3)
        gmsh.write(output_file)
        gmsh.finalize()
        
        return {
            'N_bubbles': N_bubbles,
            'R_b': R_b,
            'mesh_file': output_file
        }


# Parameter sweep example
def generate_parameter_sweep():
    """Generate meshes for all parameter combinations"""
    L = 1.0  # electrode length
    H = 0.5  # electrolyte height
    
    theta_values = [30, 60, 90, 120]  # degrees
    alpha_values = [0.1, 0.3, 0.5, 0.7, 0.85]
    
    results = []
    
    os.makedirs("meshes", exist_ok=True)
    
    for theta in theta_values:
        for alpha in alpha_values:
            print(f"Generating mesh: θ={theta}°, α={alpha}")
            
            # 2D mesh
            mesh_gen = BubbleArrayMesh(L, H, theta, alpha, mesh_size=0.02, dimension=2)
            info_2d = mesh_gen.generate_2d_mesh(
                f"meshes/mesh_2D_theta{theta}_alpha{alpha}.msh"
            )
            
            # 3D mesh (optional, more expensive)
            # mesh_gen_3d = BubbleArrayMesh(L, H, theta, alpha, mesh_size=0.05, dimension=3)
            # info_3d = mesh_gen_3d.generate_3d_mesh(
            #     f"meshes/mesh_3D_theta{theta}_alpha{alpha}.msh"
            # )
            
            results.append({
                'theta': theta,
                'alpha': alpha,
                'info_2d': info_2d
            })
    
    return results


if __name__ == "__main__":
    results = generate_parameter_sweep()
    print(f"\nGenerated {len(results)} mesh configurations")
