# LinkerForge schematic assembly -- NOT a predicted structure.
load sfgfp_mcherry_flexible_non_cleavable.pdb
bg_color white
hide everything
show cartoon, chain A or chain B
color forest, chain A
color tv_red, chain B
show spheres, chain L
set sphere_scale, 0.45, chain L
color magenta, chain L
show spheres, chain A and resi 232 and name CA
show spheres, chain B and resi 4 and name CA
color gold, (chain A and resi 232 and name CA) or (chain B and resi 4 and name CA)
orient
