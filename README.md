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

        Absolute and relative radius drop

        Radius tapering slope

        Cumulative length and tapering profile

    Saves:

        radius_metrics.csv

        radius_profile.csv

```bash
python main.py --mesh_path <path_to_your_mesh_file> --analyze
```
## Landmark selection for mesh manipulation

```bash
python main.py --mesh_path <path_to_your_mesh_file> --landmark_select
```
## Visualize different components of the brain vasculature in tiled gif
```bash
python main.py --mesh_path <path_to_your_mesh_file> --gif_tiles
```
## Visualize components saved during landmark selection
Landmark selection saves the component on which landmark is selected as a seperate mesh along with the selected point. Visualize for further analysis

```bash
python main.py --mesh_path <path_to_your_mesh_file> --landmark_select --read_saved_comp
```
If .pkl file is already generated with components and landmarks:

```bash
python main.py --mesh_path <path_to_your_mesh_file> --read_saved_comp
```






