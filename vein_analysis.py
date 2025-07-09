import trimesh
import numpy as np
import math
from vedo import Plotter, Mesh
from vedo import Sphere, show as vshow, Mesh, save
import vedo
import time
from scipy.interpolate import splprep, splev
import colorsys
import imageio.v2 as imageio
import argparse
import trimesh.voxel as vox
from collections import defaultdict 
import skeletor as sk
from skimage.morphology import skeletonize
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from vedo import Points, Lines, show, Mesh, write
import networkx as nx
from scipy.spatial.distance import euclidean
import pickle
import skeletor as sk
import trimesh


from vedo import Mesh, Points, Lines, show

def skeleton_and_branchpoints(mesh, visualize=True):
    import skeletor as sk
    import trimesh
    import numpy as np
    import networkx as nx
    from sklearn.neighbors import kneighbors_graph
    from vedo import Mesh, Points, Lines, show
    # fixed = sk.pre.fix_mesh(mesh, remove_disconnected=1, inplace=False)
    fixed = mesh
    skel = sk.skeletonize.by_wavefront(fixed, waves=3, progress=True)
    skel = sk.post.clean_up(skel)
    skel = sk.post.remove_bristles(skel, fixed, los_only=True)
    points = skel.vertices
    print(len(points))

    A = kneighbors_graph(points, n_neighbors=7, mode='distance', include_self=False)
    G = nx.from_scipy_sparse_array(A)
    tree = nx.minimum_spanning_tree(G)

    nodes_to_remove = []
    for node in tree.nodes():
        if tree.degree(node) == 1:
            neighbor = next(tree.neighbors(node))
            if tree.degree(neighbor) == 3:
                nodes_to_remove.append((node, neighbor))

    tree.remove_edges_from(nodes_to_remove)
    tree.remove_nodes_from([int(n) for n, _ in nodes_to_remove])

    print(f"Remaining nodes in tree: {len(tree.nodes())}")


    remaining_nodes = list(tree.nodes())

    old_to_new_index = {old_idx: new_idx for new_idx, old_idx in enumerate(remaining_nodes)}
    new_points = np.array([points[i] for i in remaining_nodes])
    edges = [
        [new_points[old_to_new_index[u]], new_points[old_to_new_index[v]]]
        for u, v in tree.edges()
    ]

    vedo_mesh = Mesh([mesh.vertices, mesh.faces], c="lightblue", alpha=0.3)
    nodes = Points(new_points, r=8, c='red')
    edge_lines = Lines(edges, c='blue', lw=2)


    from scipy.interpolate import splprep, splev
    import numpy as np
    import networkx as nx
    from vedo import Points, Lines


    def find_valid_branches(graph):
        branches = []
        visited = set()

        for node in graph.nodes():
            if graph.degree[node] >= 3:
                for neighbor in graph.neighbors(node):
                    if (node, neighbor) in visited or (neighbor, node) in visited:
                        continue

                    branch = [node, neighbor]
                    visited.add((node, neighbor))
                    visited.add((neighbor, node))
                    current = neighbor
                    prev = node

                    while graph.degree[current] == 2:
                        next_nodes = [n for n in graph.neighbors(current) if n != prev]
                        if not next_nodes:
                            break
                        next_node = next_nodes[0]
                        branch.append(next_node)
                        visited.add((current, next_node))
                        visited.add((next_node, current))
                        prev, current = current, next_node

                    if graph.degree[current] >= 3 and current != node:
                        branches.append(branch)

        return branches

    smoothed_lines = []
    for branch in find_valid_branches(tree):
        coords = np.array([points[i] for i in branch])
        if len(coords) < 3:
            continue  # Not enough points for a spline

        elif len(coords) == 3:
            smoothed_lines.append(coords)  # Keep it straight

        else:
            tck, u = splprep(coords.T, s=0.5)
            u_fine = np.linspace(0, 1, 50)
            smooth_coords = np.array(splev(u_fine, tck)).T
            smoothed_lines.append(smooth_coords)




    from vedo import Points as VedoPoints, Line, show
    from scipy.spatial import cKDTree

    # Step 1: All skeleton point indices from the tree
    all_tree_nodes = set(tree.nodes())

    # Step 2: Identify the ones used in splines
    spline_point_indices = set()
    for branch in find_valid_branches(tree):
        spline_point_indices.update(branch)

    if len(spline_point_indices) != 0 and len(smoothed_lines) != 0:
        print("Valid spline points found in the tree.")
        # Step 3: Get coordinates of spline-fitting skeleton points (orange)
        spline_fit_points = np.array([points[i] for i in spline_point_indices])
        all_spline_coords = np.vstack(smoothed_lines)
        spline_tree = cKDTree(all_spline_coords)

        # Step 4: Get other skeleton points (not used in splines)
        non_spline_indices = list(all_tree_nodes - spline_point_indices)
        non_spline_points = np.array([points[i] for i in non_spline_indices])

        
        # Step 6: Remove non-spline points that are too close to any spline
        radius = 0.0001  # Adjust depending on scale
        mask = np.ones(len(non_spline_points), dtype=bool)

        for i, pt in enumerate(non_spline_points):
            if spline_tree.query(pt)[0] < radius:
                mask[i] = False  # Too close to spline — mark for removal

        filtered_non_spline_points = non_spline_points[mask]
    
    else:
        print("No valid spline points found in the tree.")
        spline_fit_points = np.empty((0, 3))
        filtered_non_spline_points = np.array([points[i] for i in all_tree_nodes])


    combined_points = np.vstack([spline_fit_points, filtered_non_spline_points])
    from sklearn.neighbors import kneighbors_graph
    import networkx as nx

    A_combined = kneighbors_graph(combined_points, n_neighbors=7, mode='distance', include_self=False)
    G_combined = nx.from_scipy_sparse_array(A_combined)

    tree_combined = nx.minimum_spanning_tree(G_combined)

    import numpy as np

    # Set your desired threshold (tune this based on your scale)
    length_threshold = 0.00015  # for example

    # Find and remove long edges
    edges_to_remove = []
    for u, v in tree_combined.edges():
        p1 = combined_points[u]
        p2 = combined_points[v]
        dist = np.linalg.norm(p1 - p2)
        if dist > length_threshold:
            edges_to_remove.append((u, v))

    tree_combined.remove_edges_from(edges_to_remove)
    print(f"Removed {len(edges_to_remove)} long edges.")

    # Rebuild edges and draw
    edges_combined = [
        [combined_points[u], combined_points[v]] for u, v in tree_combined.edges()
    ]
    tree_lines = Lines(edges_combined, c='purple', lw=2)

    # Step 1: Combined points
    combined_points = np.vstack([spline_fit_points, filtered_non_spline_points])

    # Step 1: Build kNN graph
    from sklearn.neighbors import kneighbors_graph
    A_combined = kneighbors_graph(combined_points, n_neighbors=5, mode='distance', include_self=False)
    G_combined = nx.from_scipy_sparse_array(A_combined)

    # Step 2: Build MST
    tree_combined = nx.minimum_spanning_tree(G_combined)

    # Step 3: Attempt fallback connections if line is outside mesh
    kdtree = cKDTree(combined_points)
    final_edges = []
    failed_edges = []

    for u, v in tree_combined.edges():
        p1 = combined_points[u]
        p2 = combined_points[v]
        segment = np.linspace(p1, p2, num=20)

        if np.all(mesh.contains(segment)):
            final_edges.append([p1, p2])
        else:
            # Try alternative connections for u
            distances, indices = kdtree.query(p1, k=7)
            found = False
            for alt_idx in indices:
                if alt_idx == u or alt_idx == v:
                    continue
                alt_p = combined_points[alt_idx]
                alt_segment = np.linspace(p1, alt_p, num=20)
                if np.all(mesh.contains(alt_segment)):
                    final_edges.append([p1, alt_p])
                    found = True
                    break
            if not found:
                failed_edges.append([p1, p2])

    # Step 4: Visualize
    vedo_mesh = Mesh([mesh.vertices, mesh.faces], c="lightblue", alpha=0.3)
          

    import networkx as nx
    import numpy as np
    from vedo import Points as VedoPoints, Lines, Mesh, show

    # 1. Create an empty tree (not a general graph)
    tree_final = nx.Graph()

    # 2. Map points to indices
    point_to_index = {tuple(p): i for i, p in enumerate(combined_points)}

    # 3. Add edges only from final_edges
    for edge in final_edges:
        u = point_to_index[tuple(edge[0])]
        v = point_to_index[tuple(edge[1])]
        tree_final.add_edge(u, v)

    # 4. Extract 1-degree nodes (tips)
    tip_indices = [n for n in tree_final.nodes if tree_final.degree[n] == 1]
    tip_coords = np.array([combined_points[i] for i in tip_indices])

    # 5. Visualize
    vedo_mesh = Mesh([mesh.vertices, mesh.faces], c="lightblue", alpha=0.3)

    from networkx import connected_components

    components = list(nx.connected_components(tree_final))
    print(f"Found {len(components)} disconnected components.")
    component_points = []
    for comp in components:
        comp_indices = list(comp)
        comp_coords = np.array([combined_points[i] for i in comp_indices])
        component_points.append(comp_coords)
    from scipy.spatial import cKDTree

    new_connections = []
    while len(component_points) > 1:
        base = component_points[0]
        base_kdtree = cKDTree(base)
        best_pair = None
        best_dist = np.inf
        best_idx = -1

        for i in range(1, len(component_points)):
            target = component_points[i]
            distances, base_indices = base_kdtree.query(target)

            for j, dist in enumerate(distances):
                if dist < best_dist:
                    p1 = base[base_indices[j]]
                    p2 = target[j]
                    segment = np.linspace(p1, p2, num=20)
                    if np.all(mesh.contains(segment)):
                        best_dist = dist
                        best_pair = (p1, p2)
                        best_idx = i

        if best_pair:
            new_connections.append(best_pair)
            # merge this component into base
            component_points[0] = np.vstack([component_points[0], component_points[best_idx]])
            del component_points[best_idx]
        else:
            print("No valid connection found inside mesh. Stopping.")
            break
    for p1, p2 in new_connections:
        u = np.where((combined_points == p1).all(axis=1))[0][0]
        v = np.where((combined_points == p2).all(axis=1))[0][0]
        tree_final.add_edge(u, v)
    all_tree_edges = [
        [combined_points[u], combined_points[v]] for u, v in tree_final.edges()
    ]
    final_tree_lines = Lines(all_tree_edges, c='green', lw=3)
    
    # if visualize:
    #     show(vedo_mesh, final_tree_lines, Points(combined_points, r=6, c='red'), axes=1)
    components = list(nx.connected_components(tree_final))
    print(f"Found {len(components)} disconnected components.", type(combined_points), len(combined_points))
    return tree_final, combined_points, mesh, vedo_mesh 

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


def global_tortuosity(tree, points):
    tortuosities = []
    for comp in nx.connected_components(tree):
        subgraph = tree.subgraph(comp)

        # Step 1: Find longest shortest path (diameter)
        lengths = dict(nx.all_pairs_dijkstra_path_length(subgraph))
        max_len = 0
        max_pair = (None, None)
        for u in lengths:
            for v in lengths[u]:
                if lengths[u][v] > max_len:
                    max_len = lengths[u][v]
                    max_pair = (u, v)

        if max_pair[0] is None or max_pair[1] is None:
            continue

        # Step 2: Get coordinates of path
        path_nodes = nx.shortest_path(subgraph, source=max_pair[0], target=max_pair[1])
        coords = points[path_nodes]

        # Step 3: Compute actual path length
        segment_lengths = np.linalg.norm(np.diff(coords, axis=0), axis=1)
        path_length = np.sum(segment_lengths)

        # Step 4: Euclidean distance between ends
        euclidean = np.linalg.norm(coords[0] - coords[-1])

        # Step 5: Tortuosity
        if euclidean > 0:
            tort = path_length / euclidean
            tortuosities.append(tort)

    return tortuosities


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

def compute_complexity(comp):
    try:
        tree_final ,combined_points, _,_ = skeleton_and_branchpoints(comp, visualize=False)
        return len(combined_points)
    except Exception as e:
        print(f"Failed to compute skeleton for a component: {e}")
        return -1


def parse_arguments():
    parser = argparse.ArgumentParser(description='Vein mesh processing')
    parser.add_argument('--mesh_path', type=str, help='Path to the mesh file to read', default=None)
    parser.add_argument('--analyze', type = int, help='Analyze various metrics of brain. Radial Profile - 0, Tortuosity - 1 ', default=None)
    parser.add_argument('--gif_tiles', action='store_true', help='Visualize rotating components and save as GIF', default=False)
    parser.add_argument('--landmark_select', action='store_true', help='Visualize and select important landmarks on the mesh. Save the points and mesh', default=False)
    parser.add_argument('--read_saved_comp', action='store_true', help='Read and visualize data saved in landmark_select', default=False)
    parser.add_argument('--test', action='store_true', help='Run the test for new functionality', default=False)
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
        area_thresh = 0.001
        volume_thresh = 0.01

        for comp in components:
            area_ratio = comp.area / total_area
            volume_ratio = comp.volume / total_volume if comp.is_volume else 0

            if area_ratio > area_thresh or volume_ratio > volume_thresh:
                filtered.append(comp)
                filtered_c.append(comp)
        filtered_c = sorted(filtered_c, key=compute_complexity, reverse=True)
        print(f"Total filtered components: {len(filtered_c)}")


    if a.analyze == 0:
        meshes = []
        global_sampled_pts = []
        global_radii = []
        for i, comp in enumerate(filtered_c):
            tree_final, combined_points, mesh, vedo_mesh = skeleton_and_branchpoints(comp, visualize=True)
            from scipy.spatial import cKDTree
            import numpy as np

            # Build a KD-tree on the mesh surface points
            surface_points = mesh.sample(100000)  # Sample dense enough for accuracy
            surface_tree = cKDTree(surface_points)
            sampled_pts = []
            radii = []
            for u, v in tree_final.edges():
                p1 = combined_points[u]
                p2 = combined_points[v]
                segment = np.linspace(p1, p2, num=10)

                for pt in segment:
                    dist, _ = surface_tree.query(pt)
                    global_sampled_pts.append(pt)
                    global_radii.append(dist)
                    sampled_pts.append(pt)
                    radii.append(dist)
            sampled_pts = np.array(sampled_pts)
            radii = np.array(radii)
            meshes.append((comp, sampled_pts, radii))

        global_radii = np.array(global_radii)
        global_sampled_pts = np.array(global_sampled_pts)
        min_gr =  np.percentile(global_radii, 0.1)
        max_gr = np.percentile(global_radii, 99.9)
        norm = mcolors.Normalize(vmin=min_gr, vmax=max_gr)
        cmap = plt.get_cmap('jet')

        vedo_meshes = []
        for comp, sampled_pts, radii in meshes:
            vmesh = Mesh([comp.vertices, comp.faces])
            face_centers = vmesh.cell_centers().points

            sample_tree = cKDTree(sampled_pts)
            dists, idxs = sample_tree.query(face_centers)
            face_radii = radii[idxs]

            face_colors = (cmap(norm(face_radii))[:, :3] * 255).astype(np.uint8)
            vmesh.cellcolors = face_colors
            vmesh.alpha(0.6)
            vedo_meshes.append(vmesh)

        import time
        from vedo import Plotter, Text2D, Mesh
        import matplotlib.pyplot as plt
        import numpy as np

        plt_ = Plotter(title="Click a Component to Show Radius Stats", axes=1)
        info = Text2D("", pos="top-left", c='black')
        highlight_actor = None

        def on_click(evt):
            global highlight_actor
            if not evt.actor:
                return

            try:
                idx = vedo_meshes.index(evt.actor)
                radii = meshes[idx][2]

                # Remove previous highlight
                if highlight_actor:
                    plt_.remove(highlight_actor)

                # Highlight current component
                highlight_actor = evt.actor.clone()
                highlight_actor.c('red').alpha(0.3).wireframe(True).linewidth(5)
                plt_.add(highlight_actor)
                plt_.render()  # Force the highlight to render now

                # Short delay so user sees the highlight before popup
                time.sleep(0.3)

                # Now update text info
                info.text(
                    f"Component #{idx} → Radius Stats:\n"
                    f"Min:  {radii.min():.7f}\n"
                    f"Mean: {radii.mean():.7f}\n"
                    f"Max:  {radii.max():.7f}"
                )
                plt_.render()  # Update text immediately too

                # Popup histogram
                plt.figure(figsize=(6, 3))
                plt.hist(radii, bins=50, color='skyblue', edgecolor='black')
                plt.title(f"Radius Histogram - Component {idx}")
                plt.xlabel("Radius")
                plt.ylabel("Frequency")
                plt.tight_layout()
                plt.show()

            except ValueError:
                info.text("Component not found.")
                plt_.render()

        # Add scene elements
        plt_ += vedo_meshes
        plt_ += info
        plt_.add_callback("mouse click", on_click)
        plt_.show(viewup='z')

    # if a.analyze == 1:
    #     for i, comp in enumerate(filtered_c):
    #         tree_final, combined_points, mesh, vedo_mesh = skeleton_and_branchpoints(comp, visualize=True)
    #         import numpy as np
    #         import networkx as nx

    #         global_torts = global_tortuosity(tree_final, combined_points)
    #         print(f"Found {len(global_torts)} components.")
    #         print(f"Mean global tortuosity: {np.mean(global_torts):.3f}")
    #         print(f"Min: {np.min(global_torts):.3f}, Max: {np.max(global_torts):.3f}")

    #         from vedo import Lines, Plotter, Text2D
    #         import numpy as np
    #         import networkx as nx

    #         def build_component_lines_with_tortuosity(tree, points):
    #             components = list(nx.connected_components(tree))
    #             actors = []
    #             tortuosities = []

    #             for i, comp in enumerate(components):
    #                 sub = tree.subgraph(comp)
    #                 edges = list(sub.edges())
    #                 edge_lines = [[points[u], points[v]] for u, v in edges]
    #                 actor = Lines(edge_lines, c='gray', lw=2)
    #                 actor.name = f"component_{i}"
    #                 actors.append(actor)

    #                 # Tortuosity via longest path in component
    #                 lengths = dict(nx.all_pairs_dijkstra_path_length(sub))
    #                 max_len = 0
    #                 max_pair = (None, None)
    #                 for u in lengths:
    #                     for v in lengths[u]:
    #                         if lengths[u][v] > max_len:
    #                             max_len = lengths[u][v]
    #                             max_pair = (u, v)

    #                 if max_pair[0] is not None:
    #                     path_nodes = nx.shortest_path(sub, source=max_pair[0], target=max_pair[1])
    #                     coords = points[path_nodes]
    #                     path_len = np.sum(np.linalg.norm(np.diff(coords, axis=0), axis=1))
    #                     euclidean = np.linalg.norm(coords[0] - coords[-1])
    #                     tort = path_len / euclidean if euclidean > 0 else 1.0
    #                 else:
    #                     tort = 1.0

    #                 tortuosities.append(tort)

    #             return actors, tortuosities

    #         # Generate actors and tort values
    #         actors, tortuosities = build_component_lines_with_tortuosity(tree_final, combined_points)

    #         # Vedo interactive plotter
    #         plt = Plotter(title="Click a Component to Show Tortuosity", axes=1)

    #         # Text overlay
    #         info = Text2D("", pos="top-left", c='black')

    #         def on_click(evt):
    #             if not evt.actor:
    #                 return
    #             for a in actors:
    #                 a.color("gray").lw(2)  # reset others
    #             evt.actor.color("red").lw(4)
    #             idx = actors.index(evt.actor)
    #             info.text(f"Component #{idx} → Tortuosity: {tortuosities[idx]:.3f}")

    #         plt += actors
    #         plt += info
    #         plt.add_callback("mouse click", on_click)
    #         plt.show(viewup='z')

    # if a.analyze == 1:
    #     for i, comp in enumerate(filtered_c):
    #         tree_final, combined_points, mesh, vedo_mesh = skeleton_and_branchpoints(comp, visualize=False)
    #         import numpy as np
    #         import networkx as nx
    #         from vedo import Lines, Plotter, Text2D

    #         def extract_branches_between_junctions(graph):
    #             deg = dict(graph.degree())
    #             junctions = [n for n, d in deg.items() if d >= 3 or d == 1]
    #             visited_edges = set()
    #             branches = []

    #             for node in junctions:
    #                 for neighbor in graph.neighbors(node):
    #                     if (node, neighbor) in visited_edges or (neighbor, node) in visited_edges:
    #                         continue

    #                     path = [node, neighbor]
    #                     visited_edges.add((node, neighbor))
    #                     visited_edges.add((neighbor, node))
    #                     current = neighbor
    #                     prev = node

    #                     while True:
    #                         next_nodes = [n for n in graph.neighbors(current) if n != prev]
    #                         if len(next_nodes) != 1:
    #                             break
    #                         prev, current = current, next_nodes[0]
    #                         path.append(current)
    #                         visited_edges.add((prev, current))
    #                         visited_edges.add((current, prev))

    #                         if graph.degree[current] == 1 or graph.degree[current] >= 3:
    #                             break

    #                     branches.append(path)

    #             return branches

    #         branch_paths = extract_branches_between_junctions(tree_final)
    #         branch_torts = []
    #         branch_actors = []

    #         for idx, path in enumerate(branch_paths):
    #             coords = combined_points[path]
    #             if len(coords) < 2:
    #                 continue
    #             euclidean = np.linalg.norm(coords[0] - coords[-1])
    #             path_len = np.sum(np.linalg.norm(np.diff(coords, axis=0), axis=1))
    #             tort = path_len / euclidean if euclidean > 0 else 1.0
    #             branch_torts.append(tort)

    #             lines = [[coords[i], coords[i+1]] for i in range(len(coords)-1)]
    #             actor = Lines(lines, c='gray', lw=2)
    #             actor.name = f"branch_{idx}"
    #             branch_actors.append(actor)

    #         print(f"Total branches: {len(branch_actors)}")
    #         print(f"Mean tortuosity: {np.mean(branch_torts):.3f}")
    #         print(f"Min: {np.min(branch_torts):.3f}, Max: {np.max(branch_torts):.3f}")

    #         # Vedo interactive plotter
    #         plt = Plotter(title="Click Branch to Show Tortuosity", axes=1)
    #         info = Text2D("", pos="top-left", c='black')

    #         def on_click(evt):
    #             if not evt.actor:
    #                 return
    #             for a in branch_actors:
    #                 a.color("gray").lw(2)
    #             evt.actor.color("red").lw(4)
    #             idx = branch_actors.index(evt.actor)
    #             info.text(f"Branch #{idx} → Tortuosity: {branch_torts[idx]:.3f}")
    #             plt.render()

    #         plt += branch_actors
    #         plt += info
    #         plt.add_callback("mouse click", on_click)
    #         plt.show(viewup='z')

    if a.analyze == 1:
        all_branch_actors = []
        all_tortuosities = []

        mesh_id_map = {}  # maps each actor to its tortuosity index

        for comp_idx, comp in enumerate(filtered_c):
            tree_final, combined_points, mesh, vedo_mesh = skeleton_and_branchpoints(comp, visualize=False)

            def extract_branches_between_junctions(graph):
                deg = dict(graph.degree())
                junctions = [n for n, d in deg.items() if d >= 3 or d == 1]
                visited_edges = set()
                branches = []

                for node in junctions:
                    for neighbor in graph.neighbors(node):
                        if (node, neighbor) in visited_edges or (neighbor, node) in visited_edges:
                            continue

                        path = [node, neighbor]
                        visited_edges.add((node, neighbor))
                        visited_edges.add((neighbor, node))
                        current = neighbor
                        prev = node

                        while True:
                            next_nodes = [n for n in graph.neighbors(current) if n != prev]
                            if len(next_nodes) != 1:
                                break
                            prev, current = current, next_nodes[0]
                            path.append(current)
                            visited_edges.add((prev, current))
                            visited_edges.add((current, prev))

                            if graph.degree[current] == 1 or graph.degree[current] >= 3:
                                break

                        branches.append(path)

                return branches

            branches = extract_branches_between_junctions(tree_final)

            for b_idx, path in enumerate(branches):
                coords = combined_points[path]
                if len(coords) < 2:
                    continue

                euclidean = np.linalg.norm(coords[0] - coords[-1])
                path_len = np.sum(np.linalg.norm(np.diff(coords, axis=0), axis=1))
                tort = path_len / euclidean if euclidean > 0 else 1.0

                lines = [[coords[i], coords[i+1]] for i in range(len(coords)-1)]
                actor = Lines(lines, c='gray', lw=2)
                actor.name = f"comp{comp_idx}_branch{b_idx}"

                all_branch_actors.append(actor)
                all_tortuosities.append(tort)
                mesh_id_map[id(actor)] = len(all_tortuosities) - 1

        print(f"Total branches across all components: {len(all_branch_actors)}")

        # Show all in one scene
        from vedo import Plotter, Text2D

        plt = Plotter(title="Click Branch to Show Tortuosity", axes=1)
        info = Text2D("", pos="top-left", c='black')

        def on_click(evt):
            if not evt.actor:
                return
            for a in all_branch_actors:
                a.color("gray").lw(2)
            evt.actor.color("red").lw(4)

            idx = mesh_id_map.get(id(evt.actor), None)
            if idx is not None:
                tort = all_tortuosities[idx]
                info.text(f"Branch #{idx} → Tortuosity: {tort:.3f}")
                plt.render()

        plt += all_branch_actors
        plt += info
        plt.add_callback("mouse click", on_click)
        plt.show(viewup='z')

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

