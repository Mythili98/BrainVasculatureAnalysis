# Vein Mesh Processing and Analysis

This Python script provides a comprehensive pipeline for processing, analyzing, and visualizing 3D vein-like meshes using `trimesh`, `vedo`, `skeletor`, and `networkx`.
The primary applications include connected component analysis, skeletonization, branch point detection, and radius tapering metric computation. Optional functionalities like landmark selection and rotating component visualization (GIF) are also provided.

---

## Features

- **Mesh Component Splitting**: Automatically segments a mesh into its connected components.
- **Skeletonization**: Uses wavefront-based skeletonization (`skeletor`) to extract the mesh's medial axis.
- **Branch Point Detection**: Identifies nodes with degree > 2 on the skeleton graph.
- **Radius Tapering Analysis**: Computes absolute/relative radius drop, tapering slope, and cumulative profiles along skeleton paths.
- **Rotating GIF Output**: Creates a tiled GIF of all components rotating in 3D.
- **Interactive Landmark Selection**: Allows manual point selection and mesh annotation using mouse clicks.
- **CSV Export**: Saves tapering metrics and radius profiles to `.csv` files for downstream analysis.

----

## Requirements

Install all dependencies:

```bash
pip install trimesh vedo skeletor networkx scipy imageio numpy pandas
```
---

## Basic Component Processing

This will provide results of radii tapering.
This will perform skeletonization and calculates: 

Absolute and relative radius drop,

Radius tapering slope

Cumulative length and tapering profile

Saves the data to:

radius_metrics.csv and radius_profile.csv


```bash
python vein_analysis.py --mesh_path <path_to_your_mesh_file> --analyze
```
# Other metrics coming soon..
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




