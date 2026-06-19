# Interactive Feature Analysis on 3D Brain Vascular Meshes

This Python script provides a comprehensive pipeline for processing, analyzing, and visualizing 3D vein-like meshes using `trimesh`, `vedo`, `skeletor`, and `networkx`.
The primary applications include connected component analysis, skeletonization, branch point detection, and radius tapering metric computation. Optional functionalities like landmark selection and rotating component visualization (GIF) are also provided.

---

## Features

- **Mesh Component Splitting**: Automatically segments a mesh into its connected components.
- **Skeletonization**: Uses wavefront-based skeletonization (`skeletor`) for mesh contaction and uses B-Spline method for filtering and centerline creation.
- **Radial  Analysis**: Computes radial profile of the mesh.
- **Tortuosity Analysis**: Computed tortuosity for different branches of the brain vasculature.
- **Rotating GIF Output**: Creates a tiled GIF of all components rotating in 3D.
- **Interactive Landmark Selection**: Allows manual point selection and mesh annotation using mouse clicks.

----

## Requirements

Install all dependencies:

```bash
pip install trimesh vedo skeletor networkx scipy imageio numpy pandas plotly
```
---

## Basic Component Processing

This will perform mesh contaction, spline fitting and centerline extraction to provide mesh wise radial and tortuosiiity features.

The profiles are displayed in an interactive fashion for faster analysis and inference. Click on the components to see the global features.

### For radial profile
```bash
python vein_analysis.py --mesh_path <path_to_your_mesh_file> --analyze 0 
```

### For tortuosity profile
```bash
python vein_analysis.py --mesh_path <path_to_your_mesh_file> --analyze 1
```

### To save structural information
```bash
python vein_analysis.py --mesh_path <path_to_your_mesh_file> --structural
```
It saves the,
sampled_points: points along skeleton segments
and radius values: local radius (distance to surface) at those sample points.

It also saves the mesh with the face_radii --> saves in recomputing the face mapping with radial values.

## Landmark selection for mesh manipulation

While selecting landmark, after selection is done, instead of traversing through all the mesh components, press 'Esc' to exit the plotter and save the results. Saved landmarks with the associated mesh is available as 'all_saved_comp.pkl' file. These files can be used for mesh manipulation and analysis in a collaborative fashion.

```bash
python vein_analysis.py --mesh_path <path_to_your_mesh_file> --landmark_select
```
## Visualize different components of the brain vasculature in tiled gif
```bash
python vein_analysis.py --mesh_path <path_to_your_mesh_file> --gif_tiles
```
## Visualize components saved during landmark selection
Landmark selection saves the component on which landmark is selected as a seperate mesh along with the selected point. Visualize for further analysis

```bash
python vein_analysis.py --mesh_path <path_to_your_mesh_file> --landmark_select --read_saved_comp
```
If .pkl file is already generated with components and landmarks:

```bash
python vein_analysis.py --mesh_path <path_to_your_mesh_file> --read_saved_comp
```

## Acknowledgements
Skeletor: https://github.com/navis-org/skeletor/tree/master




