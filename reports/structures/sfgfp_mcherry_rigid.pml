# LinkerForge PyMOL Visualization Script
# Construct: sfGFP_2B3P (1-238) — Linker (239-256) — mCherry_2H5Q (257-492)
bg_color white
hide everything, all
show cartoon, all
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
select domain_a, resi 1-238
select port_a, resi 236-238
select designed_linker, resi 239-256
color magenta, designed_linker
show sticks, designed_linker and (name CA+C+N+O)
select port_b, resi 257-259
select domain_b, resi 257-492
color forest, domain_a
color tv_red, domain_b
color gold, port_a or port_b
show spheres, (port_a or port_b) and name CA
set sphere_scale, 0.65, (port_a or port_b) and name CA
deselect
orient
zoom all, 4
