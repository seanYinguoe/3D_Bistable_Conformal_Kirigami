# Selected target meshes

These files are unchanged copies of the author's target geometries from the existing repository. The selection includes `quarter_dome`, `double_dome` and `hemisphere`; the planar `circle.obj` is retained for the original 2D driver.

Each `*_flat.obj` stores the curved target vertices (`v`), planar UV coordinates (`vt`), and face indices. The UV maps were supplied with the research project and are described by the original driver as BFF outputs. `vertice_sort` aligns UV ordering with target connectivity.

Do not assume that OBJ coordinates have an intrinsic unit. Use a consistent length convention across grid size, ligament geometry and material parameters. The source manifest records hashes and the exact original paths within the repository.
