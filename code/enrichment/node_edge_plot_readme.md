# Creating the BrainNetViewer plot docs

## Creating node and edge files

The file `02-oneway_righties.R` creates:

- a list of nodes to keep (*n = 44*), and writes it to a file: `enrichment_nodes_to_keep.tsv`.
  - `node`: Node name
  - `degree`: Node degree (>1)
- a (symmetric) matrix of edges (*44x44*), where each cell:
  - `0`: no edge
  - `!=0`: weight value

The file `create_node_file.R` reads in all of the Glasser nodes and creates a
`.node` file for all of them: `glasser_all.node`. For this file, both `color`
and `size` are set to 1.

Then, `create_node_file.R` filters the whole list down to the 44 nodes, which
are part of nine networks (including the excluded OAN), and uses the `degree`
value as size.

Final file: `enrichment_result.node` and `enrichment_abs_result.edge`.

## Using BrainNetViewer

**Use MATLAB R2024a**, the changes in R2025 have broken BNV.

### Networks/colors

Using color labels from Ji et al. (2019), except OAN since it was dropped.
Uses `colors.txt`, a TSV with one row per ROI.

1. Auditory - `#F98BFC`
2. Cingulo Opercular - `#C350BF`
3. Default - `#FF4E31`
4. Dorsal Attn. - `#1EEE3A`
5. Frontoparietal - `#FDF66D`
6. Language - `#05A0A1`
7. Orbito-Affective - `#7f7f7f`
8. Somatomotor - `#62FAFD`
9. Visual2 - `#425BFB`

- Use surface: `BrainMesh_ICBM152.nv`.
- To load colors > Options > Color > Modular > Load Custom Color

View angle: Azimuth: -100, Elevation: 20
Edge size: Scale: 0.80 / abs value.
Edze color: Neg: blue #425bfb; Pos: red #FB425B

Settings saved to `bnv_final.mat`.

## References

- Ji et al. (2019). https://doi.org/10.1016/j.neuroimage.2018.10.006.
