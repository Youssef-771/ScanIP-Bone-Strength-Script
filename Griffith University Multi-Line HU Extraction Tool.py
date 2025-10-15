from simpleware.scripting import *
import os, math, random
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse

# ---------------- Helpers ----------------
def to_real_point3d(p):
    return RealPoint3D(p.GetX(), p.GetY(), p.GetZ())

def delete_old_measurements(doc, names):
    ms = doc.GetMeasurements()
    sv = StringVector()
    for n in names:
        sv.append(n)
    try:
        ms.DeleteMeasurements(sv)
    except:
        for n in names:
            try:
                ms.DeleteMeasurement(n)
            except:
                pass

def get_counts_from_any(obj):
    try:
        return obj.GetVoxelCountX(), obj.GetVoxelCountY(), obj.GetVoxelCountZ()
    except:
        try:
            doc = obj if hasattr(obj, "GetDimensions") else App.GetInstance().GetActiveDocument()
            if doc:
                dims = doc.GetDimensions()
                return dims.GetVoxelCountX(), dims.GetVoxelCountY(), dims.GetVoxelCountZ()
        except:
            pass
    return None, None, None

def get_spacing_mm(doc):
    dims = doc.GetDimensions()
    return dims.GetSpacingX(), dims.GetSpacingY(), dims.GetSpacingZ()

def get_background_info(bg):
    try:
        return bg.GetPixelType(), bg.GetName()
    except:
        return None, "Unknown Background"

def clamp(v, vmin, vmax):
    return max(vmin, min(v, vmax))

def clear_previous_annotations(doc):
    try:
        annotations = doc.GetAnnotations()
        all_annotations = annotations.GetAll()
        for a in all_annotations:
            try:
                name = a.GetName()
                if name and name.startswith("HU_"):
                    annotations.RemoveAnnotation(a)  
            except:
                continue
    except:
        pass

def _norm(v):
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def _unit(v):
    l = _norm(v)
    if l < 1e-12:
        return (0.0, 0.0, 0.0)
    return (v[0]/l, v[1]/l, v[2]/l)

def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def _plane_axes_from_orientation(orientation):
    try:
        s = str(orientation).lower()
        if "axial" in s or "xy" in s:
            return 0, 1, 2
        if "coronal" in s or "xz" in s:
            return 0, 2, 1
        if "sagittal" in s or "yz" in s:
            return 1, 2, 0
    except:
        pass
    try:
        i = int(orientation)
        if i == 0: return 0, 1, 2
        if i == 1: return 0, 2, 1
        if i == 2: return 1, 2, 0
    except:
        pass
    return 0, 1, 2

def hu_to_density(hu):
    return 0.001029 * hu + 0.114259

def density_to_E(rho):
    return 6850 * (rho ** 1.49)

COLOUR_PALETTE = [
    (255, 0, 0), (0, 255, 0), (0, 0, 255), (255, 128, 0),
    (153, 0, 204), (0, 179, 179), (255, 255, 0), (255, 0, 255),
]

def get_colour_for_index(i):
    if i < len(COLOUR_PALETTE):
        rgb = COLOUR_PALETTE[i]
    else:
        rgb = (random.randint(0,255), random.randint(0,255), random.randint(0,255))
    return Colour(rgb[0], rgb[1], rgb[2])

# --------------- Single Point Circle Command -----------------
class ExtractCircleHUCommand(Command):
    def __init__(self, doc, background, center, circle_id=1, n_radial=8, 
                 radius_mm=10.0, radial_samples=20):
        Command.__init__(self)
        self.doc = doc
        self.background = background
        self.center = center
        self.circle_id = circle_id
        self.n_radial = n_radial
        self.radius_mm = float(radius_mm)
        self.radial_samples = radial_samples
        self.center_data = None
        self.radial_data = {}

    def GetName(self):
        return f"Extract Circle HU+E at point {self.circle_id}"

    def Do(self):
        self._sample_center()
        self._sample_radial_lines()
        self._create_annotations()
        return True

    def _sample_center(self):
        cx, cy, cz = self.center.GetX(), self.center.GetY(), self.center.GetZ()
        hu = self._sample_hu_at_point(cx, cy, cz)
        if not (isinstance(hu, float) and math.isnan(hu)):
            rho = hu_to_density(hu)
            E = density_to_E(rho)
            self.center_data = (cx, cy, cz, hu, rho, E)
        else:
            self.center_data = (cx, cy, cz, float('nan'), float('nan'), float('nan'))

    def _sample_radial_lines(self):
        cx, cy, cz = self.center.GetX(), self.center.GetY(), self.center.GetZ()
        sx_mm, sy_mm, sz_mm = get_spacing_mm(self.doc)
        orientation = self.doc.GetActiveSliceView().GetOrientation()
        axA, axB, axN = _plane_axes_from_orientation(orientation)
        spacing = [sx_mm, sy_mm, sz_mm]
        
        for radial_id in range(self.n_radial):
            angle = 2.0 * math.pi * radial_id / float(self.n_radial)
            cos_a, sin_a = math.cos(angle), math.sin(angle)
            radial_values = []
            
            for r_sample in range(self.radial_samples + 1):
                t = r_sample / float(self.radial_samples)
                distance_mm = t * self.radius_mm
                offset = [0.0, 0.0, 0.0]
                offset[axA] = (distance_mm * cos_a) / spacing[axA]
                offset[axB] = (distance_mm * sin_a) / spacing[axB]
                
                rx = cx + offset[0]
                ry = cy + offset[1]
                rz = cz + offset[2]
                
                hu = self._sample_hu_at_point(rx, ry, rz)
                if not (isinstance(hu, float) and math.isnan(hu)):
                    rho = hu_to_density(hu)
                    E = density_to_E(rho)
                    radial_values.append((rx, ry, rz, hu, rho, E))
                else:
                    radial_values.append((rx, ry, rz, float('nan'), float('nan'), float('nan')))
            
            self.radial_data[radial_id] = radial_values

    def _sample_hu_at_point(self, x, y, z):
        try:
            nx, ny, nz = get_counts_from_any(self.background)
            ix, iy, iz = int(round(x)), int(round(y)), int(round(z))
            if nx is not None: ix = clamp(ix, 0, nx-1)
            if ny is not None: iy = clamp(iy, 0, ny-1)
            if nz is not None: iz = clamp(iz, 0, nz-1)
            pixel_type, _ = get_background_info(self.background)
            if pixel_type and 'Float' in str(pixel_type):
                return self.background.GetVoxelValueReal(ix, iy, iz)
            else:
                try:
                    return self.background.GetVoxelValueReal(ix, iy, iz)
                except:
                    return float(self.background.GetVoxelValue(ix, iy, iz))
        except:
            return float('nan')

    def _create_annotations(self):
        annotations = self.doc.GetAnnotations()
        orientation = self.doc.GetActiveSliceView().GetOrientation()
        
        all_slices = IntVector()
        try:
            total = self.doc.GetSliceCount(orientation)
            for i in range(total):
                all_slices.append(i)
        except:
            all_slices.append(self.doc.GetActiveSliceView().GetActiveSlice())
        
        colour = get_colour_for_index(self.circle_id - 1)
        self._draw_circle(annotations, orientation, all_slices, colour)
        self._draw_radial_lines(annotations, orientation, all_slices, colour)
        self._add_degree_labels(annotations, orientation, all_slices)
        
        try:
            cx, cy, cz = self.center.GetX(), self.center.GetY(), self.center.GetZ()
            _, _, _, hu, _, E = self.center_data
            txt = f"Circle {self.circle_id}\nRadius: {self.radius_mm:.1f}mm\n"
            txt += f"Center HU: {hu:.1f}\nCenter E: {E:.0f} MPa"
            sv = IntVector()
            sv.append(self.doc.GetActiveSliceView().GetActiveSlice())
            annotations.AddTextBox(f"HU_Circle_Text_{self.circle_id}", orientation, sv,
                                 RealPoint3D(cx + 10, cy + 10, cz), RealPoint3D(cx + 100, cy + 60, cz),
                                 RealSize(120, 60), txt, 1.0, True, False)
        except:
            pass

    def _draw_circle(self, annotations, orientation, all_slices, colour):
        cx, cy, cz = self.center.GetX(), self.center.GetY(), self.center.GetZ()
        sx_mm, sy_mm, sz_mm = get_spacing_mm(self.doc)
        axA, axB, axN = _plane_axes_from_orientation(orientation)
        spacing = [sx_mm, sy_mm, sz_mm]
        num_segments = 36
        
        for seg in range(num_segments):
            a1 = 2.0 * math.pi * seg / num_segments
            a2 = 2.0 * math.pi * (seg + 1) / num_segments
            
            if axA == 0 and axB == 1:
                p1 = [cx + (self.radius_mm / spacing[0]) * math.cos(a1), 
                      cy + (self.radius_mm / spacing[1]) * math.sin(a1), cz]
                p2 = [cx + (self.radius_mm / spacing[0]) * math.cos(a2), 
                      cy + (self.radius_mm / spacing[1]) * math.sin(a2), cz]
            elif axA == 0 and axB == 2:
                p1 = [cx + (self.radius_mm / spacing[0]) * math.cos(a1), cy,
                      cz + (self.radius_mm / spacing[2]) * math.sin(a1)]
                p2 = [cx + (self.radius_mm / spacing[0]) * math.cos(a2), cy,
                      cz + (self.radius_mm / spacing[2]) * math.sin(a2)]
            else:
                p1 = [cx, cy + (self.radius_mm / spacing[1]) * math.cos(a1),
                      cz + (self.radius_mm / spacing[2]) * math.sin(a1)]
                p2 = [cx, cy + (self.radius_mm / spacing[1]) * math.cos(a2),
                      cz + (self.radius_mm / spacing[2]) * math.sin(a2)]
            
            seg_name = f"HU_Circle_{self.circle_id}_seg_{seg}"
            line = annotations.AddLine(seg_name, orientation, all_slices,
                                     RealPoint3D(p1[0], p1[1], p1[2]),
                                     RealPoint3D(p2[0], p2[1], p2[2]), True, True)
            line.SetColour(colour)
            line.SetWidth(2.0)

    def _draw_radial_lines(self, annotations, orientation, all_slices, colour):
        cx, cy, cz = self.center.GetX(), self.center.GetY(), self.center.GetZ()
        sx_mm, sy_mm, sz_mm = get_spacing_mm(self.doc)
        axA, axB, axN = _plane_axes_from_orientation(orientation)
        spacing = [sx_mm, sy_mm, sz_mm]
        
        for radial_id in range(self.n_radial):
            angle = 2.0 * math.pi * radial_id / float(self.n_radial)
            
            if axA == 0 and axB == 1:
                end_point = [cx + (self.radius_mm / spacing[0]) * math.cos(angle),
                           cy + (self.radius_mm / spacing[1]) * math.sin(angle), cz]
            elif axA == 0 and axB == 2:
                end_point = [cx + (self.radius_mm / spacing[0]) * math.cos(angle), cy,
                           cz + (self.radius_mm / spacing[2]) * math.sin(angle)]
            else:
                end_point = [cx, cy + (self.radius_mm / spacing[1]) * math.cos(angle),
                           cz + (self.radius_mm / spacing[2]) * math.sin(angle)]
            
            line_name = f"HU_Circle_Radial_{self.circle_id}_{radial_id}"
            line = annotations.AddLine(line_name, orientation, all_slices,
                                     RealPoint3D(cx, cy, cz),
                                     RealPoint3D(end_point[0], end_point[1], end_point[2]),
                                     True, True)
            line.SetColour(colour)
            line.SetWidth(1.5)

    def _add_degree_labels(self, annotations, orientation, all_slices):
        cx, cy, cz = self.center.GetX(), self.center.GetY(), self.center.GetZ()
        sx_mm, sy_mm, sz_mm = get_spacing_mm(self.doc)
        axA, axB, axN = _plane_axes_from_orientation(orientation)
        spacing = [sx_mm, sy_mm, sz_mm]
        # Position labels at the line end
        label_radius_mm = self.radius_mm
        
        for radial_id in range(self.n_radial):
            angle = 2.0 * math.pi * radial_id / float(self.n_radial)
            angle_deg = radial_id * 360.0 / float(self.n_radial)
            
            # Calculate label position at the exact end of the radial line
            if axA == 0 and axB == 1:
                label_point = [cx + (label_radius_mm / spacing[0]) * math.cos(angle),
                             cy + (label_radius_mm / spacing[1]) * math.sin(angle), cz]
            elif axA == 0 and axB == 2:
                label_point = [cx + (label_radius_mm / spacing[0]) * math.cos(angle), cy,
                             cz + (label_radius_mm / spacing[2]) * math.sin(angle)]
            else:
                label_point = [cx, cy + (label_radius_mm / spacing[1]) * math.cos(angle),
                             cz + (label_radius_mm / spacing[2]) * math.sin(angle)]
            
            try:
                # Use the label position as both anchor and corner for minimal box
                text_box = annotations.AddTextBox(
                    f"HU_Circle_Deg_{self.circle_id}_{radial_id}",
                    orientation,
                    all_slices,
                    RealPoint3D(label_point[0], label_point[1], label_point[2]),
                    RealPoint3D(label_point[0] + 20, label_point[1] + 15, label_point[2]),
                    RealSize(40, 20),
                    f"{angle_deg:.0f}°",
                    1.0,
                    True,
                    False
                )
                # Make visible in 3D view
                text_box.SetShowIn3DView(True)
                # Set colors for visibility
                text_box.SetForegroundColour(Colour(255, 255, 0))  # Yellow text
                text_box.SetBackgroundColour(Colour(0, 0, 0))      # Black background
                text_box.SetDrawBackground(True)
                text_box.SetDrawBorder(True)
                text_box.SetBorderColour(Colour(0, 255, 255))      # Cyan border
            except:
                pass

    def CanUndo(self):
        return False

    def OnNativeDelete(self):
        try:
            self.baseline_values = None
            self.radial_data = None
            self.doc = None
            self.background = None
        except:
            pass

# --------------- Single Point Circle User Action -----------------
class ExtractCircleHUUserAction(UserAction):
    def __init__(self):
        UserAction.__init__(self, "Measurements", "Single point circle analysis", "Circle HU Tool")
        tool_icon = os.path.join(App.GetInstance().GetInstallLocation(), "Plugins", "GU Logo.png")
        if os.path.isfile(tool_icon):
            self.SetBitmapFileName(tool_icon)
    
    def OnActivated(self):
        app = App.GetInstance()
        doc = app.GetActiveDocument()
        if doc is None: return
        background = doc.GetActiveBackground()
        if background is None: return
        clear_previous_annotations(doc)

        try:
            n_radial = InputDialog.GetInteger("Number of radial lines:", "8")
        except:
            n_radial = 8
        try:
            radius_mm = float(InputDialog.GetInteger("Circle radius (mm):", "10"))
        except:
            radius_mm = 10.0
        try:
            radial_samples = InputDialog.GetInteger("Samples per radial line:", "20")
        except:
            radial_samples = 20

        center = doc.PickPoint()
        if center is None: return

        cmd = ExtractCircleHUCommand(doc, background, center, circle_id=1, n_radial=n_radial,
                                    radius_mm=radius_mm, radial_samples=radial_samples)
        doc.SubmitCommand(cmd)

        save_dir = InputDialog.ChooseDirectory("Choose folder to save CSV and plots",
                                              os.path.expanduser("~"))
        if not save_dir:
            app.ShowMessage("Export cancelled.", "Cancelled")
            return

        self._export_circle_data(cmd, save_dir, doc, background)
        app.ShowMessage(f"Exported to:\n{save_dir}", "Complete")

    def _export_circle_data(self, cmd, save_dir, doc, background):
        csv_path = os.path.join(save_dir, "circle_center_data.csv")
        with open(csv_path, 'w') as f:
            f.write("X,Y,Z,HU,Density[g/cm3],E[MPa]\n")
            cx, cy, cz, hu, rho, E = cmd.center_data
            f.write(f"{cx:.3f},{cy:.3f},{cz:.3f},{hu:.3f},{rho:.3f},{E:.2f}\n")

        radial_csv = os.path.join(save_dir, "circle_radial_data.csv")
        with open(radial_csv, 'w') as f:
            f.write("RadialID,Sample,Distance[mm],X,Y,Z,HU,Density[g/cm3],E[MPa]\n")
            for radial_id, radial_values in cmd.radial_data.items():
                for sample_idx, (rx, ry, rz, hu, rho, E) in enumerate(radial_values):
                    distance_mm = (sample_idx / float(cmd.radial_samples)) * cmd.radius_mm
                    f.write(f"{radial_id},{sample_idx},{distance_mm:.3f},{rx:.3f},{ry:.3f},{rz:.3f},"
                           f"{hu:.3f},{rho:.3f},{E:.2f}\n")

        self._create_combined_plot(cmd, save_dir, doc, background)
        self._create_radial_profile_plots(cmd, save_dir)

    def _create_combined_plot(self, cmd, save_dir, doc, background):
        fig = plt.figure(figsize=(20, 7))
        orientation = doc.GetActiveSliceView().GetOrientation()
        axA, axB, axN = _plane_axes_from_orientation(orientation)
        cx, cy, cz = cmd.center_data[0], cmd.center_data[1], cmd.center_data[2]
        center_coords = [cx, cy, cz]
        slice_idx = int(round(center_coords[axN]))
        
        orientation_names = {0: "Axial", 1: "Coronal", 2: "Sagittal"}
        try:
            orient_str = str(orientation).lower()
            if "axial" in orient_str or "xy" in orient_str:
                view_name = "Axial"
            elif "coronal" in orient_str or "xz" in orient_str:
                view_name = "Coronal"
            elif "sagittal" in orient_str or "yz" in orient_str:
                view_name = "Sagittal"
            else:
                view_name = orientation_names.get(int(orientation), "Unknown")
        except:
            view_name = "Unknown"
        
        try:
            nx, ny, nz = get_counts_from_any(background)
            sx_mm, sy_mm, sz_mm = get_spacing_mm(doc)
            spacing = [sx_mm, sy_mm, sz_mm]
            
            if axN == 2:
                ct_slice = np.zeros((ny, nx))
                for y in range(ny):
                    for x in range(nx):
                        try:
                            ct_slice[y, x] = background.GetVoxelValueReal(x, y, slice_idx)
                        except:
                            ct_slice[y, x] = 0
                extent_x, extent_y = (0, nx), (0, ny)
                aspect_ratio = (spacing[1] / spacing[0])
            elif axN == 1:
                ct_slice = np.zeros((nz, nx))
                for z in range(nz):
                    for x in range(nx):
                        try:
                            ct_slice[z, x] = background.GetVoxelValueReal(x, slice_idx, z)
                        except:
                            ct_slice[z, x] = 0
                extent_x, extent_y = (0, nx), (0, nz)
                aspect_ratio = (spacing[2] / spacing[0])
            else:
                ct_slice = np.zeros((nz, ny))
                for z in range(nz):
                    for y in range(ny):
                        try:
                            ct_slice[z, y] = background.GetVoxelValueReal(slice_idx, y, z)
                        except:
                            ct_slice[z, y] = 0
                extent_x, extent_y = (0, ny), (0, nz)
                aspect_ratio = (spacing[2] / spacing[1])
        except:
            ct_slice = None
            aspect_ratio = 1.0
        
       # ---------------- CT slice with zoom ----------------
        ax1 = plt.subplot(1, 3, 1)
        if ct_slice is not None:
            ax1.imshow(ct_slice, cmap='gray', origin='lower',
                extent=[extent_x[0], extent_x[1], extent_y[0], extent_y[1]],
                aspect=aspect_ratio)
    
        radius_vox_a = cmd.radius_mm / spacing[axA]
        radius_vox_b = cmd.radius_mm / spacing[axB]
    
        circle = Ellipse((center_coords[axA], center_coords[axB]), 
                    width=2*radius_vox_a, height=2*radius_vox_b,
                    fill=False, color='cyan', linewidth=2)
        ax1.add_patch(circle)
    
        label_frequency = 2 if cmd.n_radial <= 12 else 4
        for radial_id in range(cmd.n_radial):
            angle = 2.0 * math.pi * radial_id / float(cmd.n_radial)
            angle_deg = radial_id * 360.0 / float(cmd.n_radial)
        
            end_x = center_coords[axA] + radius_vox_a * math.cos(angle)
            end_y = center_coords[axB] + radius_vox_b * math.sin(angle)
            ax1.plot([center_coords[axA], end_x], [center_coords[axB], end_y],
                    'c-', linewidth=1, alpha=0.7)
        
            if radial_id % label_frequency == 0:
                label_distance = 2.5
                label_x = center_coords[axA] + radius_vox_a * label_distance * math.cos(angle)
                label_y = center_coords[axB] + radius_vox_b * label_distance * math.sin(angle)
                line_start_x = center_coords[axA] + radius_vox_a * 1.05 * math.cos(angle)
                line_start_y = center_coords[axB] + radius_vox_b * 1.05 * math.sin(angle)
                ax1.plot([line_start_x, label_x], [line_start_y, label_y],
                        'yellow', linewidth=1.5, alpha=0.8)
                ax1.text(label_x, label_y, f'{angle_deg:.0f}°',   # <-- use normal angle here
                        color='yellow', fontsize=10, ha='center', va='center',
                        fontweight='bold',
                        bbox=dict(boxstyle='round,pad=0.5', facecolor='black', alpha=0.8, edgecolor='cyan'))

        # Zooms to circle region
        zoom_factor = 4.0   
        x_min = center_coords[axA] - zoom_factor * radius_vox_a
        x_max = center_coords[axA] + zoom_factor * radius_vox_a
        y_min = center_coords[axB] - zoom_factor * radius_vox_b
        y_max = center_coords[axB] + zoom_factor * radius_vox_b
        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(y_min, y_max)
    
        ax1.set_title(f'CT Slice with Sampling Circle\n{view_name} View - Slice {slice_idx}', 
                    fontsize=11, fontweight='bold')
        ax1.set_xlabel('Voxels')
        ax1.set_ylabel('Voxels')

        
        # ---------------- Polar plot for Elastic Modulus ----------------
        ax2 = plt.subplot(1, 3, 2, projection='polar')
        angles, E_means, E_stds = [], [], []
        for radial_id in range(cmd.n_radial):
            angle = 2.0 * math.pi * radial_id / float(cmd.n_radial)
            angles.append(angle)
            radial_values = cmd.radial_data[radial_id]
            E_vals = [E for _, _, _, _, _, E in radial_values 
                     if not (isinstance(E, float) and math.isnan(E))]
            if E_vals:
                E_means.append(np.mean(E_vals))
                E_stds.append(np.std(E_vals))
            else:
                E_means.append(0)
                E_stds.append(0)
        
        angles.append(angles[0])
        E_means.append(E_means[0])
        E_stds.append(E_stds[0])
        angles = np.array(angles)
        E_means = np.array(E_means)
        E_stds = np.array(E_stds)
        
        ax2.plot(angles, E_means, 'o-', linewidth=2, markersize=6, color='blue', label='Mean E')
        ax2.fill_between(angles, E_means - 0.25*E_stds, E_means + 0.25*E_stds,
                         alpha=0.3, color='blue', label='±0.25 STD')
        ax2.fill(angles, E_means, alpha=0.15, color='blue')
        ax2.set_theta_zero_location("E")
        ax2.set_theta_direction("counterclockwise")
        
        for i in range(cmd.n_radial):
            angle = 2.0 * math.pi * i / float(cmd.n_radial)
            angle_deg = i * 360.0 / float(cmd.n_radial)
            label_radius = max(E_means) * 1.2 if max(E_means) > 0 else 1000
            ax2.text(angle, label_radius, f'{angle_deg:.0f}°',
                   ha='center', va='center', fontsize=9, fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
        
        ax2.set_title(f"Elastic Modulus Distribution\n(Radius = {cmd.radius_mm:.1f} mm)",
                     va='bottom', pad=20, fontsize=11, fontweight='bold')
        ax2.set_ylabel('E [MPa]', labelpad=30, fontsize=20)
        ax2.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
        ax2.grid(True)
        
        # ---------------- Polar plot for HU ----------------
        ax3 = plt.subplot(1, 3, 3, projection='polar')
        angles_hu, HU_means, HU_stds = [], [], []
        for radial_id in range(cmd.n_radial):
            angle = 2.0 * math.pi * radial_id / float(cmd.n_radial)
            angles_hu.append(angle)
            radial_values = cmd.radial_data[radial_id]
            HU_vals = [hu for _, _, _, hu, _, _ in radial_values 
                      if not (isinstance(hu, float) and math.isnan(hu))]
            if HU_vals:
                HU_means.append(np.mean(HU_vals))
                HU_stds.append(np.std(HU_vals))
            else:
                HU_means.append(0)
                HU_stds.append(0)
        
        angles_hu.append(angles_hu[0])
        HU_means.append(HU_means[0])
        HU_stds.append(HU_stds[0])
        angles_hu = np.array(angles_hu)
        HU_means = np.array(HU_means)
        HU_stds = np.array(HU_stds)
        
        ax3.plot(angles_hu, HU_means, 'o-', linewidth=2, markersize=6, color='darkgreen', label='Mean HU')
        ax3.fill_between(angles_hu, HU_means - 0.25*HU_stds, HU_means + 0.25*HU_stds,
                         alpha=0.3, color='green', label='±0.25 STD')
        ax3.fill(angles_hu, HU_means, alpha=0.15, color='green')
        ax3.set_theta_zero_location("E")
        ax3.set_theta_direction("counterclockwise")
        
        for i in range(cmd.n_radial):
            angle = 2.0 * math.pi * i / float(cmd.n_radial)
            angle_deg = i * 360.0 / float(cmd.n_radial)
            label_radius = max(HU_means) * 1.2 if max(HU_means) > 0 else 100
            ax3.text(angle, label_radius, f'{angle_deg:.0f}°',
                   ha='center', va='center', fontsize=9, fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
        
        ax3.set_title(f"HU Distribution\n(Radius = {cmd.radius_mm:.1f} mm)",
                     va='bottom', pad=20, fontsize=11, fontweight='bold')
        ax3.set_ylabel('HU Value', labelpad=30, fontsize=20)
        ax3.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
        ax3.grid(True)
        
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, "circle_polar_plot_with_CT.png"), dpi=150, bbox_inches='tight')
        plt.close()

    def _create_radial_profile_plots(self, cmd, save_dir):
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        distances_mm = np.linspace(0, cmd.radius_mm, cmd.radial_samples + 1)
        
        ax1 = axes[0, 0]
        for radial_id in range(cmd.n_radial):
            radial_values = cmd.radial_data[radial_id]
            hu_vals = [hu for _, _, _, hu, _, _ in radial_values]
            angle_deg = (radial_id / float(cmd.n_radial)) * 360
            ax1.plot(distances_mm, hu_vals, alpha=0.6, label=f'{angle_deg:.0f}°')
        ax1.set_xlabel('Distance from Center [mm]')
        ax1.set_ylabel('HU Value')
        ax1.set_title('HU Profiles Along All Radial Lines')
        ax1.grid(True, alpha=0.3)
        if cmd.n_radial <= 8:
            ax1.legend(fontsize=8)
        
        ax2 = axes[0, 1]
        hu_means, hu_stds = [], []
        for sample_idx in range(cmd.radial_samples + 1):
            sample_hus = []
            for radial_id in range(cmd.n_radial):
                hu = cmd.radial_data[radial_id][sample_idx][3]
                if not (isinstance(hu, float) and math.isnan(hu)):
                    sample_hus.append(hu)
            if sample_hus:
                hu_means.append(np.mean(sample_hus))
                hu_stds.append(np.std(sample_hus))
            else:
                hu_means.append(np.nan)
                hu_stds.append(0)
        ax2.plot(distances_mm, hu_means, 'b-', linewidth=2, label='Mean HU')
        ax2.fill_between(distances_mm, np.array(hu_means) - 0.25*np.array(hu_stds),
                         np.array(hu_means) + 0.25*np.array(hu_stds),
                         alpha=0.3, color='blue', label='±0.25 STD')
        ax2.set_xlabel('Distance from Center [mm]')
        ax2.set_ylabel('HU Value')
        ax2.set_title('Average HU Profile (All Directions)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        ax3 = axes[1, 0]
        for radial_id in range(cmd.n_radial):
            radial_values = cmd.radial_data[radial_id]
            E_vals = [E for _, _, _, _, _, E in radial_values]
            angle_deg = (radial_id / float(cmd.n_radial)) * 360
            ax3.plot(distances_mm, E_vals, alpha=0.6, label=f'{angle_deg:.0f}°')
        ax3.set_xlabel('Distance from Center [mm]')
        ax3.set_ylabel('Elastic Modulus [MPa]')
        ax3.set_title('Elastic Modulus Profiles Along All Radial Lines')
        ax3.grid(True, alpha=0.3)
        if cmd.n_radial <= 8:
            ax3.legend(fontsize=8)
        
        ax4 = axes[1, 1]
        E_means, E_stds = [], []
        for sample_idx in range(cmd.radial_samples + 1):
            sample_Es = []
            for radial_id in range(cmd.n_radial):
                E = cmd.radial_data[radial_id][sample_idx][5]
                if not (isinstance(E, float) and math.isnan(E)):
                    sample_Es.append(E)
            if sample_Es:
                E_means.append(np.mean(sample_Es))
                E_stds.append(np.std(sample_Es))
            else:
                E_means.append(np.nan)
                E_stds.append(0)
        ax4.plot(distances_mm, E_means, 'r-', linewidth=2, label='Mean E')
        ax4.fill_between(distances_mm, np.array(E_means) - 0.25*np.array(E_stds),
                         np.array(E_means) + 0.25*np.array(E_stds),
                         alpha=0.3, color='red', label='±0.25 STD')
        ax4.set_xlabel('Distance from Center [mm]')
        ax4.set_ylabel('Elastic Modulus [MPa]')
        ax4.set_title('Average Elastic Modulus Profile (All Directions)')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, "circle_radial_profiles.png"), dpi=150)
        plt.close()

# --------------- Multi-Line User Action -----------------
class ExtractHUUserAction(UserAction):
    def __init__(self):
        UserAction.__init__(self, "Measurements", "Multi-line measurements", "Multi-Line HU Tool")
        tool_icon = os.path.join(App.GetInstance().GetInstallLocation(), "Plugins", "GU Logo.png")
        if os.path.isfile(tool_icon):
            self.SetBitmapFileName(tool_icon)
            
    def OnActivated(self):  
        app = App.GetInstance()
        doc = app.GetActiveDocument()
        if doc is None: return
        background = doc.GetActiveBackground()
        if background is None: return
        clear_previous_annotations(doc)

        try:
            n_lines = InputDialog.GetInteger("How many main lines?", "1")
        except:
            return
        try:
            n_radial = InputDialog.GetInteger("Radial lines per point (0=none):", "0")
        except:
            n_radial = 0

        radial_length_mm = 10.0
        radial_samples = 20
        if n_radial > 0:
            try:
                radial_length_mm = float(InputDialog.GetInteger("Radius for sampling/circles (mm):", "10"))
                radial_samples = InputDialog.GetInteger("Samples per radial line:", "20")
            except:
                pass

        all_commands = []
        for line_id in range(1, n_lines+1):
            start = doc.PickPoint()
            if start is None: break
            end = doc.PickPoint()
            if end is None: break

            cmd = ExtractHUCommand(doc, background, start, end, line_id=line_id,
                                 n_radial=n_radial, radial_length_mm=radial_length_mm,
                                 radial_samples=radial_samples)
            doc.SubmitCommand(cmd)
            all_commands.append(cmd)

        if not all_commands:
            return

        save_dir = InputDialog.ChooseDirectory("Choose folder to save CSV and plots",
                                              os.path.expanduser("~"))
        if not save_dir:
            app.ShowMessage("Export cancelled.", "Cancelled")
            return

        self._export_data(all_commands, save_dir, n_radial > 0)
        app.ShowMessage(f"Exported to:\n{save_dir}", "Complete")

    def _export_data(self, commands, save_dir, has_radials):
        csv_path = os.path.join(save_dir, "hu_density_E_all_lines.csv")
        with open(csv_path, 'w') as f:
            f.write("LineID,SampleIdx,X,Y,Z,HU,Density[g/cm3],E[MPa],Valid\n")
            for cmd in commands:
                for line_id, x,y,z,hu,rho,E,sample_idx in cmd.baseline_values:
                    is_valid = "TRUE" if not (isinstance(hu, float) and math.isnan(hu)) else "FALSE"
                    hu_str = f"{hu:.3f}" if is_valid == "TRUE" else "NaN"
                    rho_str = f"{rho:.3f}" if is_valid == "TRUE" else "NaN"
                    E_str = f"{E:.2f}" if is_valid == "TRUE" else "NaN"
                    f.write(f"{line_id},{sample_idx},{x:.3f},{y:.3f},{z:.3f},{hu_str},{rho_str},{E_str},{is_valid}\n")

        if has_radials:
            radial_csv = os.path.join(save_dir, "radial_data_all_lines.csv")
            with open(radial_csv, 'w') as f:
                f.write("LineID,SampleIdx,RadialID,RadialSample,X,Y,Z,HU,Density[g/cm3],E[MPa],Valid\n")
                for cmd in commands:
                    for sample_idx, radial_dict in cmd.radial_data.items():
                        for radial_id, radial_values in radial_dict.items():
                            for r_sample, (rx,ry,rz,hu,rho,E) in enumerate(radial_values):
                                is_valid = "TRUE" if not (isinstance(hu, float) and math.isnan(hu)) else "FALSE"
                                hu_str = f"{hu:.3f}" if is_valid == "TRUE" else "NaN"
                                rho_str = f"{rho:.3f}" if is_valid == "TRUE" else "NaN"
                                E_str = f"{E:.2f}" if is_valid == "TRUE" else "NaN"
                                f.write(f"{cmd.line_id},{sample_idx},{radial_id},{r_sample},{rx:.3f},{ry:.3f},{rz:.3f},{hu_str},{rho_str},{E_str},{is_valid}\n")

        self._create_visualizations(commands, save_dir, has_radials)

    def _create_visualizations(self, commands, save_dir, has_radials):
        plt.figure(figsize=(12,8))
        plt.subplot(2,1,1)
        for cmd in commands:
            hu_vals = [hu for _,_,_,_,hu,_,_,_ in cmd.baseline_values if not (isinstance(hu, float) and math.isnan(hu))]
            xs = list(range(len(hu_vals)))
            c = get_colour_for_index(cmd.line_id-1)
            color = [c.GetRed()/255.0, c.GetGreen()/255.0, c.GetBlue()/255.0]
            plt.plot(xs, hu_vals, marker='o', linestyle='-', markersize=3,
                     label=f"Line {cmd.line_id}", color=color)
        plt.xlabel("Sample Number")
        plt.ylabel("HU Value")
        plt.title("HU Values Along Multiple Lines")
        plt.legend()
        plt.grid(True, alpha=0.3)

        plt.subplot(2,1,2)
        for cmd in commands:
            E_vals = [E for _,_,_,_,_,_,E,_ in cmd.baseline_values if not (isinstance(E, float) and math.isnan(E))]
            xs = list(range(len(E_vals)))
            c = get_colour_for_index(cmd.line_id-1)
            color = [c.GetRed()/255.0, c.GetGreen()/255.0, c.GetBlue()/255.0]
            plt.plot(xs, E_vals, marker='o', linestyle='-', markersize=3,
                     label=f"Line {cmd.line_id}", color=color)
        plt.xlabel("Sample Number")
        plt.ylabel("Elastic Modulus [MPa]")
        plt.title("Elastic Modulus (E) Along Multiple Lines")
        plt.legend()
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, "baseline_plots.png"), dpi=150)
        plt.close()

        if has_radials:
            self._create_radial_visualizations(commands, save_dir)

    def _create_radial_visualizations(self, commands, save_dir):
        for cmd in commands:
            if not cmd.radial_data:
                continue
            self._create_cylindrical_core_plot(cmd, save_dir)
            self._create_polar_plot(cmd, save_dir)

    def _create_cylindrical_core_plot(self, cmd, save_dir):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))
        sample_indices = sorted(cmd.radial_data.keys())
        if not sample_indices:
            plt.close()
            return

        # ---------- HU DATA ----------
        baseline_averages_hu, baseline_stds_hu = [], []
        for sample_idx in sample_indices:
            radial_dict = cmd.radial_data.get(sample_idx, {})
            all_hu = []
            for radial_id in radial_dict:
                for _, _, _, hu, _, _ in radial_dict[radial_id]:
                    if not (isinstance(hu, float) and math.isnan(hu)):
                        all_hu.append(hu)
            if all_hu:
                baseline_averages_hu.append(np.mean(all_hu))
                baseline_stds_hu.append(np.std(all_hu))
            else:
                baseline_averages_hu.append(np.nan)
                baseline_stds_hu.append(np.nan)

        baseline_hu = []
        for sample_idx in sample_indices:
            baseline_val = None
            for _, _, _, _, hu, _, _, s_idx in cmd.baseline_values:
                if s_idx == sample_idx and not (isinstance(hu, float) and math.isnan(hu)):
                    baseline_val = hu
                    break
            baseline_hu.append(baseline_val if baseline_val is not None else np.nan)

        # Plot HU
        ax1.plot(sample_indices, baseline_hu, 'b-o', label='Baseline HU', linewidth=2)
        ax1.plot(sample_indices, baseline_averages_hu, 'r-s', label='Radial Average', linewidth=2)
        ax1.fill_between(sample_indices, np.array(baseline_averages_hu) - np.array(baseline_stds_hu),
                        np.array(baseline_averages_hu) + np.array(baseline_stds_hu),
                        alpha=0.3, color='red', label='±1 STD')
        ax1.set_xlabel('Position Along Baseline')
        ax1.set_ylabel('HU Value')
        ax1.set_title('HU Variation Along Core')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # ---------- ELASTIC MODULUS DATA ----------
        baseline_averages_E, baseline_stds_E = [], []
        for sample_idx in sample_indices:
            radial_dict = cmd.radial_data.get(sample_idx, {})
            all_E = []
            for radial_id in radial_dict:
                for _, _, _, _, _, E in radial_dict[radial_id]:
                    if not (isinstance(E, float) and math.isnan(E)):
                        all_E.append(E)
            if all_E:
                baseline_averages_E.append(np.mean(all_E))
                baseline_stds_E.append(np.std(all_E))
            else:
                baseline_averages_E.append(np.nan)
                baseline_stds_E.append(np.nan)

        baseline_E = []
        for sample_idx in sample_indices:
            baseline_val = None
            for _, _, _, _, _, _, E, s_idx in cmd.baseline_values:
                if s_idx == sample_idx and not (isinstance(E, float) and math.isnan(E)):
                    baseline_val = E
                    break
            baseline_E.append(baseline_val if baseline_val is not None else np.nan)

        # Plot E
        ax2.plot(sample_indices, baseline_E, 'b-o', label='Baseline E', linewidth=2)
        ax2.plot(sample_indices, baseline_averages_E, 'r-s', label='Radial Average', linewidth=2)
        ax2.fill_between(sample_indices, np.array(baseline_averages_E) - np.array(baseline_stds_E),
                        np.array(baseline_averages_E) + np.array(baseline_stds_E),
                        alpha=0.3, color='red', label='±1 STD')
        ax2.set_xlabel('Position Along Baseline')
        ax2.set_ylabel('Elastic Modulus [MPa]')
        ax2.set_title('Elastic Modulus Variation Along Core')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, f"line_{cmd.line_id}_core.png"), dpi=150)
        plt.close()

    def _create_polar_plot(self, cmd, save_dir):
        if not cmd.radial_data:
            return
        first_sample = next(iter(cmd.radial_data))
        n_radials = len(cmd.radial_data[first_sample])
        if n_radials == 0:
            return

        angles = np.linspace(0, 2*np.pi, n_radials, endpoint=False)
        hu_means, hu_stds = [], []

        for radial_id in range(n_radials):
            all_hu_vals = []
            for sample_idx, radial_dict in cmd.radial_data.items():
                radial_values = radial_dict.get(radial_id, [])
                hu_vals = [v[3] for v in radial_values if not (isinstance(v[3], float) and math.isnan(v[3]))]
                if hu_vals:
                    all_hu_vals.extend(hu_vals)
            if all_hu_vals:
                hu_means.append(np.mean(all_hu_vals))
                hu_stds.append(np.std(all_hu_vals))
            else:
                hu_means.append(np.nan)
                hu_stds.append(0)

        angles = np.append(angles, angles[0])
        hu_means = np.append(hu_means, hu_means[0])
        hu_stds = np.append(hu_stds, hu_stds[0])

        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111, polar=True)
        ax.plot(angles, hu_means, 'o-', linewidth=2, label=f"Line {cmd.line_id}")
        ax.fill_between(angles, hu_means - hu_stds, hu_means + hu_stds, alpha=0.2, label='±1 STD')
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_title(f"Line {cmd.line_id} - HU Distribution")
        ax.legend()

        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, f"line_{cmd.line_id}_polar.png"), dpi=150)
        plt.close()

def main():
    app = App.GetInstance()
    
    action1 = ExtractHUUserAction()
    action2 = ExtractCircleHUUserAction()
    
    try:
        app.AddUserAction(action1)
        app.AddUserAction(action2)
        app.ShowMessage(
            "TWO HU Analysis Tools Registered:\n\n"
            "1. MULTI-LINE TOOL - Draw multiple measurement lines with optional radial sampling\n"
            "2. CIRCLE TOOL - Click single point for circular radial analysis\n\n"
            "Circle Tool Features:\n"
            "✓ Visible radial lines from center to edge\n"
            "✓ Combined CT + E + HU polar plots\n"
            "✓ Radial profile plots (HU and E vs distance)\n"
            "✓ Complete CSV export\n\n"
            "Both tools ready to use!",
            "Success"
        )
    except:
        pass

if __name__ == "__main__":
    main()

# --------------- Multi-Line Command -----------------
class ExtractHUCommand(Command):
    def __init__(self, doc, background, start, end, line_id=1, num_samples=100,
                 n_radial=0, radial_length_mm=10.0, radial_samples=20):
        Command.__init__(self)
        self.doc = doc
        self.background = background
        self.start = start
        self.end = end
        self.num_samples = num_samples
        self.line_id = line_id
        self.n_radial = n_radial
        self.radial_length_mm = float(radial_length_mm)
        self.radial_samples = radial_samples
        self.baseline_values = []
        self.radial_data = {}

    def GetName(self):
        return f"Extract HU+E with radials for line {self.line_id}"

    def Do(self):
        delete_old_measurements(self.doc, [
            f"HU_Line_{self.line_id}",
            f"HU_Start_{self.line_id}",
            f"HU_End_{self.line_id}",
            f"HU_Text_{self.line_id}"
        ])
        self._sample_baseline()
        if self.n_radial > 0:
            self._sample_radial_lines()
        self._create_annotations()
        return True

    def _sample_baseline(self):
        pixel_type, _ = get_background_info(self.background)
        nx, ny, nz = get_counts_from_any(self.background)
        sx, sy, sz = self.start.GetX(), self.start.GetY(), self.start.GetZ()
        ex, ey, ez = self.end.GetX(), self.end.GetY(), self.end.GetZ()

        for i in range(self.num_samples + 1):
            t = i / float(self.num_samples)
            x = sx + t * (ex - sx)
            y = sy + t * (ey - sy)
            z = sz + t * (ez - sz)
            ix, iy, iz = int(round(x)), int(round(y)), int(round(z))
            if nx is not None: ix = clamp(ix, 0, nx-1)
            if ny is not None: iy = clamp(iy, 0, ny-1)
            if nz is not None: iz = clamp(iz, 0, nz-1)

            try:
                if pixel_type and 'Float' in str(pixel_type):
                    hu = self.background.GetVoxelValueReal(ix, iy, iz)
                else:
                    try:
                        hu = self.background.GetVoxelValueReal(ix, iy, iz)
                    except:
                        hu = float(self.background.GetVoxelValue(ix, iy, iz))
                rho = hu_to_density(hu)
                E = density_to_E(rho)
                self.baseline_values.append((self.line_id, x, y, z, hu, rho, E, i))
            except:
                self.baseline_values.append((self.line_id, x, y, z, float('nan'), float('nan'), float('nan'), i))

    def _sample_radial_lines(self):
        sx, sy, sz = self.start.GetX(), self.start.GetY(), self.start.GetZ()
        ex, ey, ez = self.end.GetX(), self.end.GetY(), self.end.GetZ()
        sx_mm, sy_mm, sz_mm = get_spacing_mm(self.doc)
        dx_w = (ex - sx) * sx_mm
        dy_w = (ey - sy) * sy_mm
        dz_w = (ez - sz) * sz_mm
        main_w = _unit((dx_w, dy_w, dz_w))
        if _norm(main_w) < 1e-12:
            return

        up = (0.0, 0.0, 1.0) if abs(main_w[2]) < 0.9 else (0.0, 1.0, 0.0)
        perp1_w = _unit(_cross(main_w, up))
        if _norm(perp1_w) < 1e-12:
            up = (1.0, 0.0, 0.0)
            perp1_w = _unit(_cross(main_w, up))
        perp2_w = _unit(_cross(main_w, perp1_w))

        for sample_idx, base in enumerate(self.baseline_values):
            _, bx, by, bz, hu0, _, _, _ = base
            if hu0 is None or (isinstance(hu0, float) and math.isnan(hu0)):
                continue
            self.radial_data[sample_idx] = {}

            for radial_id in range(self.n_radial):
                angle = 2.0 * math.pi * radial_id / float(self.n_radial)
                cos_a, sin_a = math.cos(angle), math.sin(angle)
                rad_dir_w = (
                    cos_a * perp1_w[0] + sin_a * perp2_w[0],
                    cos_a * perp1_w[1] + sin_a * perp2_w[1],
                    cos_a * perp1_w[2] + sin_a * perp2_w[2],
                )
                radial_values = []
                half_len_mm = self.radial_length_mm

                for r_sample in range(self.radial_samples + 1):
                    t = (r_sample / float(self.radial_samples) - 0.5) * 2.0
                    off_wx = t * half_len_mm * rad_dir_w[0]
                    off_wy = t * half_len_mm * rad_dir_w[1]
                    off_wz = t * half_len_mm * rad_dir_w[2]
                    rx = bx + off_wx / sx_mm
                    ry = by + off_wy / sy_mm
                    rz = bz + off_wz / sz_mm

                    hu = self._sample_hu_at_point(rx, ry, rz)
                    if not (isinstance(hu, float) and math.isnan(hu)):
                        rho = hu_to_density(hu)
                        E = density_to_E(rho)
                        radial_values.append((rx, ry, rz, hu, rho, E))
                    else:
                        radial_values.append((rx, ry, rz, float('nan'), float('nan'), float('nan')))
                self.radial_data[sample_idx][radial_id] = radial_values

    def _sample_hu_at_point(self, x, y, z):
        try:
            nx, ny, nz = get_counts_from_any(self.background)
            ix, iy, iz = int(round(x)), int(round(y)), int(round(z))
            if nx is not None: ix = clamp(ix, 0, nx-1)
            if ny is not None: iy = clamp(iy, 0, ny-1)
            if nz is not None: iz = clamp(iz, 0, nz-1)
            pixel_type, _ = get_background_info(self.background)
            if pixel_type and 'Float' in str(pixel_type):
                return self.background.GetVoxelValueReal(ix, iy, iz)
            else:
                try:
                    return self.background.GetVoxelValueReal(ix, iy, iz)
                except:
                    return float(self.background.GetVoxelValue(ix, iy, iz))
        except:
            return float('nan')

    def _create_annotations(self):
        annotations = self.doc.GetAnnotations()
        orientation = self.doc.GetActiveSliceView().GetOrientation()
        slices = self.doc.GetSliceIndices(orientation)
        colour = get_colour_for_index(self.line_id - 1)
        line_ann = annotations.AddLine(
            f"HU_Line_{self.line_id}", orientation, slices,
            to_real_point3d(self.start), to_real_point3d(self.end), True, True)
        line_ann.SetColour(colour)
        line_ann.SetWidth(2.0)

        if self.n_radial > 0 and self.radial_data:
            self._add_radial_circles(annotations, orientation, colour)

        try:
            dims = self.doc.GetDimensions()
            dx_mm = (self.end.GetX() - self.start.GetX()) * dims.GetSpacingX()
            dy_mm = (self.end.GetY() - self.start.GetY()) * dims.GetSpacingY()
            dz_mm = (self.end.GetZ() - self.start.GetZ()) * dims.GetSpacingZ()
            length_mm = math.sqrt(dx_mm*dx_mm + dy_mm*dy_mm + dz_mm*dz_mm)
            midx = (self.start.GetX() + self.end.GetX()) / 2.0
            midy = (self.start.GetY() + self.end.GetY()) / 2.0
            midz = (self.start.GetZ() + self.end.GetZ()) / 2.0
            txt = f"{length_mm:.2f} mm"
            if self.n_radial > 0:
                txt += f" ({self.n_radial} radials, radius {self.radial_length_mm:.1f} mm)"
            sv = IntVector()
            sv.append(self.doc.GetActiveSliceView().GetActiveSlice())
            annotations.AddTextBox(f"HU_Text_{self.line_id}", orientation, sv,
                                 RealPoint3D(midx, midy, midz), RealPoint3D(midx+30, midy+15, midz),
                                 RealSize(80, 20), txt, 1.0, True, False)
        except:
            pass

    def _add_radial_circles(self, annotations, orientation, base_colour):
        try:
            all_slices = IntVector()
            try:
                total = self.doc.GetSliceCount(orientation)
                for i in range(total):
                    all_slices.append(i)
            except:
                all_slices.append(self.doc.GetActiveSliceView().GetActiveSlice())

            sx_mm, sy_mm, sz_mm = get_spacing_mm(self.doc)
            spacing = (sx_mm, sy_mm, sz_mm)
            radius_mm = float(self.radial_length_mm)
            axA, axB, axN = _plane_axes_from_orientation(orientation)
            r_vox = [0.0, 0.0, 0.0]
            r_vox[axA] = radius_mm / spacing[axA]
            r_vox[axB] = radius_mm / spacing[axB]

            total_samples = len(self.baseline_values) - 1
            if total_samples < 0:
                return
            stations = [0, max(0, total_samples // 2), max(0, total_samples)]
            station_colour = {
                stations[0]: Colour(255, 0, 0),
                stations[1]: Colour(255, 255, 0),
                stations[2]: Colour(0, 255, 0),
            }

            if radius_mm < 10:
                num_segments = 12
            elif radius_mm < 20:
                num_segments = 24
            else:
                num_segments = 36

            for s_idx in stations:
                center = None
                for _, x, y, z, hu, _, _, idx in self.baseline_values:
                    if idx == s_idx:
                        center = [x, y, z]
                        break
                if center is None:
                    continue
                col = station_colour.get(s_idx, base_colour)

                for seg in range(num_segments):
                    a1 = 2.0 * math.pi * seg / num_segments
                    a2 = 2.0 * math.pi * (seg + 1.0) / num_segments
                    p1 = center[:]
                    p2 = center[:]
                    p1[axA] = center[axA] + r_vox[axA] * math.cos(a1)
                    p1[axB] = center[axB] + r_vox[axB] * math.sin(a1)
                    p2[axA] = center[axA] + r_vox[axA] * math.cos(a2)
                    p2[axB] = center[axB] + r_vox[axB] * math.sin(a2)
                    seg_name = f"HU_RadialCircle_{self.line_id}_{s_idx}_seg_{seg}"
                    line = annotations.AddLine(seg_name, orientation, all_slices,
                                             RealPoint3D(p1[0], p1[1], p1[2]),
                                             RealPoint3D(p2[0], p2[1], p2[2]), True, True)
                    line.SetColour(col)
                    line.SetWidth(2.0)
        except:
            pass

    def CanUndo(self):
        return False
    
    def OnNativeDelete(self):
        try:
            self.center_data = None
            self.radial_data = None
            self.doc = None
            self.background = None
        except:
            pass