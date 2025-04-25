import trimesh
import numpy as np
import math
from vedo import Plotter, Mesh
from vedo import Sphere, show as vshow, Mesh
import time
import colorsys
import imageio.v2 as imageio
import argparse
import trimesh.voxel as vox
import numpy as np
from collections import defaultdict 
import skeletor as sk
import csv
from vedo import Points, Lines, show, Mesh, write
import networkx as nx
from scipy.spatial.distance import euclidean
import pickle

def radius_tapering_metric(skel, bp_coords):
    G = nx.Graph()
    radius_list = skel.swc['radius'].tolist()
    coords = skel.skeleton.vertices

    # Build graph
    for e in skel.skeleton.entities:
        v1, v2 = e.points
        G.add_edge(v1, v2, weight=euclidean(coords[v1], coords[v2]))

    degrees = dict(G.degree())
    branch_indices = [idx for idx in degrees if degrees[idx] > 2]
    endpoint_indices = [idx for idx in degrees if degrees[idx] == 1]

    taper_metrics = []
    taper_profiles = []

    for i,bp in enumerate(branch_indices):
        for ep in branch_indices[i+1:]+endpoint_indices:
            try:
                path = nx.shortest_path(G, source=bp, target=ep, weight='weight')
                radii_path = [radius_list[i] for i in path]
                coords_path = [coords[i] for i in path]

                if len(radii_path) < 2:
                    continue

                start_r = radii_path[0]
                end_r = radii_path[-1]
                rel_drop = (start_r - end_r) / (start_r + 1e-8)  # avoid divide by 0
                abs_drop = start_r - end_r

                # Path length
                path_len = sum(euclidean(coords_path[i], coords_path[i+1]) for i in range(len(coords_path)-1))
                slope = abs_drop / (path_len + 1e-8)

                taper_metrics.append({
                    "bp_index": bp,
                    "ep_index": ep,
                    "start_radius": start_r,
                    "end_radius": end_r,
                    "abs_drop": abs_drop,
                    "rel_drop": rel_drop,
                    "path_len": path_len,
                    "slope": slope,
                    "path": path
                })
                profile = []
                cum_dist = 0.0
                for j in range(len(path)):
                    idx = path[j]
                    rad = radius_list[idx]

                    if j > 0:
                        dist = euclidean(coords[path[j]], coords[path[j-1]])
                        cum_dist += dist
                    else:
                        cum_dist = 0.0

                    profile.append((cum_dist, rad))

                taper_profiles.append({
                    "start_bp": bp,
                    "end_bp": ep,
                    "profile": profile,
                    "path": path
                })
            except nx.NetworkXNoPath:
                continue

    return taper_metrics, taper_profiles

def skeleton_and_branchpoints(mesh, visualize=True):
    import skeletor as sk

    # Step 1: Skeletonize the mesh
    fixed = sk.pre.fix_mesh(mesh, remove_disconnected=1, inplace=False)
    skel = sk.skeletonize.by_wavefront(fixed, waves=1, progress=True)
    skel = sk.post.clean_up(skel)
    skel = sk.post.remove_bristles(skel, fixed)
    sk.post.radii(skel, fixed)
    skeleton = skel.skeleton  # trimesh.path.Path3D

    # Step 2: Build adjacency graph to find branch points
    adj = {}
    for entity in skeleton.entities:
        v1_idx, v2_idx = entity.points
        for a, b in [(v1_idx, v2_idx), (v2_idx, v1_idx)]:
            adj.setdefault(a, set()).add(b)

    # Step 3: Identify branch points (degree > 2)metric in taper_info:
            
    branch_indices = [idx for idx, neighbors in adj.items() if len(neighbors) > 2]
    branch_coords = skeleton.vertices[branch_indices]

    if visualize:
        # Original mesh in translucent sky blue
        vedo_mesh = trimesh_to_vedo(mesh)

        # Skeleton lines
        line_start = [skeleton.vertices[e.points[0]] for e in skeleton.entities]
        line_end = [skeleton.vertices[e.points[1]] for e in skeleton.entities]
        lines = Lines(line_start, line_end, c='green', lw=2)

        # Branch points
        branch_pts = Points(branch_coords, r=10, c='red')

        # Display everything
        show(vedo_mesh, lines, branch_pts, axes=1, title="Skeleton and Branch Points")
    return skel, branch_coords

def generate_distinct_colors(n):
    colors = []
    for i in range(n):
        h = i / n
        r, g, b = colorsys.hsv_to_rgb(h, 1.0, 1.0)
        colors.append([int(255*r), int(255*g), int(255*b), 255])
    return colors

def split_mesh_into_connected_components(mesh):
    labels = trimesh.graph.connected_component_labels(mesh.face_adjacency)
    components = []
    for label in np.unique(labels):
        faces = np.where(labels == label)[0]
        component = mesh.submesh([faces], only_watertight=False)[0]
        components.append(component)
    return components

def trimesh_to_vedo(tri):
    return Mesh([tri.vertices, tri.faces], c='lightblue', alpha=0.5)

def recenter_and_scale_mesh(mesh, target_size=1.0):
    centroid = mesh.center_mass
    mesh.vertices -= centroid
    bbox = mesh.bounding_box.extents
    scale = target_size / np.max(bbox)
    mesh.vertices *= scale
    return mesh

def tiled_plot_with_rotation(components, cols=7, rot_frames=100, gif_path="rotation.gif"):
    rows = math.ceil(len(components) / cols)
    plotter = Plotter(shape=(rows, cols), title="Rotating Components", axes=0, interactive=False)

    actors = []
    for i, comp in enumerate(components):
        comp = recenter_and_scale_mesh(comp)
        vm = trimesh_to_vedo(comp)
        plotter.show(vm, at=i, resetcam=True)
        actors.append(vm)

    plotter.interactive = True

    frames = []
    for _ in range(rot_frames):
        for actor in actors:
            actor.rotate_y(5)
        plotter.show(resetcam=False)
        screenshot = plotter.screenshot(asarray=True)
        frames.append(screenshot)
        time.sleep(0.05)

    imageio.mimsave(gif_path, frames, duration=0.05)
    print(f"GIF saved to: {gif_path}")

    plotter.show(interactive=True)



def parse_arguments():
    parser = argparse.ArgumentParser(description='Vein mesh processing')
    parser.add_argument('--mesh_path', type=str, help='Path to the mesh file to read', default=None)
    parser.add_argument('--analyze', action='store_true', help='Analyze various metrics of brain', default=False)
    parser.add_argument('--gif_tiles', action='store_true', help='Visualize rotating components and save as GIF', default=False)
    parser.add_argument('--landmark_select', action='store_true', help='Visualize and select important landmarks on the mesh. Save the points and mesh', default=False)
    parser.add_argument('--read_saved_comp', action='store_true', help='Read and visualize data saved in landmark_select', default=False)
    return parser.parse_args()

if __name__ == '__main__':
    a = parse_arguments()
    if a.mesh_path:
        mesh_path = a.mesh_path
        mesh = trimesh.load_mesh(mesh_path)

        components = split_mesh_into_connected_components(mesh)
                
        print(f"Total connected components: {len(components)}")

        total_area = sum([c.area for c in components])
        total_volume = sum([c.volume for c in components if c.is_volume])

        filtered = []
        filtered_c = []
        area_thresh = 0.0005
        volume_thresh = 0.005

        for comp in components:
            area_ratio = comp.area / total_area
            volume_ratio = comp.volume / total_volume if comp.is_volume else 0

            if area_ratio > area_thresh or volume_ratio > volume_thresh:
                filtered.append(comp)

        distinct_colors = generate_distinct_colors(len(filtered))
        for comp, color in zip(filtered, distinct_colors):
            comp.visual.face_colors = np.tile(color, (len(comp.faces), 1))
            filtered_c.append(comp)

        print(f"Total filtered components: {len(filtered_c)}")
    
    if a.analyze:
        csv_path = "radius_metrics.csv"
        csv_path_prof = "radius_profile.csv"
        metric_per_branch  = [] ## per component
        profile_per_branch = [] ## per component
        for i, comp in enumerate(filtered_c):
            skel, bp_coords = skeleton_and_branchpoints(comp, visualize=True)
            [taper_info, taper_profile] = radius_tapering_metric(skel, bp_coords)
            for metric in taper_info:
                metric_per_branch.append(metric)
            for profile in taper_profile:
                profile_per_branch.append(profile)
        print(f"Total metrics: {len(metric_per_branch)}")
        print(f"Total metrics: {len(profile_per_branch)}")
        with open(csv_path, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=[
                "bp_index", "ep_index", "start_radius", "end_radius",
                "abs_drop", "rel_drop", "path_len", "slope", "path"
            ])
            writer.writeheader()
            for row in metric_per_branch:
                row_copy = row.copy()
                row_copy['path'] = str(row_copy['path'])  # Convert list to string for CSV
                writer.writerow(row_copy)

        print(f"Tapering metrics saved to {csv_path}")
        with open(csv_path_prof, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=[
                "start_bp", "end_bp", "profile", "path"
            ])
            writer.writeheader()
            for row in profile_per_branch:
                row_copy = row.copy()
                row_copy['path'] = str(row_copy['path'])  # Convert list to string for CSV
                writer.writerow(row_copy)

        print(f"Tapering profiles saved to {csv_path}")


                
    
    if a.gif_tiles:
        tiled_plot_with_rotation(filtered_c[:70], gif_path="vein_rotation.gif")
    
    if a.landmark_select:
        from vedo import Points, show as vshow, Mesh

        selected_results = []
        esc_pressed = [False]
        # Create a Plotter for handling interactive visualization
        plotter = Plotter(title="Click on Each Mesh to Select a Point")

        # Function to handle mouse click events
        def on_click(event):
            if event.actor and event.picked3d is not None:
                clicked_point = np.array(event.picked3d)
                print("Selected point:", clicked_point)

                # Highlight the clicked point with a black sphere
                highlight = Sphere(pos=clicked_point, r=0.00001, c='black')
                plotter.add(highlight)
                selected_results.append((clicked_point, event.actor))

        def on_key(event):
            if event.keypress.lower() == "escape":
                print("ESC pressed. Exiting visualization.")
                esc_pressed[0] = True
                plotter.close()

        plotter.add_callback("left mouse click", on_click)
        plotter.add_callback("key press", on_key)
        # Show meshes iteratively
        for mesh in filtered_c:
            if esc_pressed[0]:
                break
            # Convert mesh to vedo object
            m = trimesh_to_vedo(mesh)

            # Add mesh to plotter
            plotter.add(m)

            # Show the mesh and wait for clicks
            plotter.show(interactive=True)

        mesh_point_map = {}
        def add_clicked_point(vedo_actor, clicked_point):
            mesh_id = id(vedo_actor)  # use memory address to identify unique meshes
            
            if mesh_id not in mesh_point_map:
                mesh_point_map[mesh_id] = {
                    "mesh": vedo_actor,
                    "points": []
                }
            
            mesh_point_map[mesh_id]["points"].append(np.array(clicked_point))
        for i, (point, vedo_actor) in enumerate(selected_results):
        #     print(f"Selected point {point} from mesh {vedo_actor}")
            add_clicked_point(vedo_actor, point)



        all_meshes_data = {}

        for mesh_id, data in mesh_point_map.items():
            vedo_mesh = data["mesh"]
            clicked_points = np.array(data["points"])
            print(vedo_mesh)
            vertices = vedo_mesh.points
            print(vertices.shape)
            faces = vedo_mesh.cells
            print(len(faces))
            
            all_meshes_data[str(mesh_id)] = {
                "vertices": vertices,
                "faces": faces,
                "selected_points": clicked_points
            }

        # Save to file
        with open("all_meshes_with_points.pkl", "wb") as f:
            pickle.dump(all_meshes_data, f)

        print(f"Saved {len(all_meshes_data)} unique mesh(es) with clicked points.")

    if a.read_saved_comp:
        with open("all_meshes_with_points.pkl", "rb") as f:
            data = pickle.load(f)

        for mesh_id, content in data.items():
            mesh = Mesh([content["vertices"], content["faces"]])
            points = Points(content["selected_points"], r=10, c="red")
            show(mesh, points, axes=1, viewup="z", interactive=True)

