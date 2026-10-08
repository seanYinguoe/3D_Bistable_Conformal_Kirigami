# Target meshes and supplied UV maps

The retained examples are `quarter_dome`, `double_dome` and `hemisphere`. Each `*_flat.obj` already contains both curved target vertices (`v`) and planar UV coordinates (`vt`), together with face indices. Separate raw OBJ copies are unnecessary for the retained workflows.

These are unchanged inputs from `Xiaoyuan_bistable_kirigami_unit/Input_model/Reference_model/<name>/` at source commit `afb8ffa`. The original driver describes the UV maps as Boundary First Flattening outputs; `vertice_sort` aligns their indexing. BFF is not bundled.

OBJ files do not specify physical units. Keep grid size, ligament geometry, SVG dimensions and material units consistent.

| File | SHA-256 |
|---|---|
| `double_dome/double_dome_flat.obj` | `f22d84bc883434a0d6403c7ffd044bf696ba7e5645e9ad0ee30ddd5102fd1070` |
| `hemisphere/hemisphere_flat.obj` | `5ffefc8cb6c32563d7239114e319ab015b6d28ce59af77816e1abf5bab4601bc` |
| `quarter_dome/quarter_dome_flat.obj` | `25c394b83e34e388159147ff4a63d6ee4937ae1b9202c16ed007a66ca53a4d7a` |
