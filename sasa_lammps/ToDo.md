# ToDo

- Currently, the test_against_vmd_reference_data tests are failing for all parameters. I already checked: The newly added vdW-radii do not make a significant difference. This basically leaves only our new C extension implementation. My approach is to first generate more test data from the original vmd-python using the generate_vmd_reference.py. I could for example use a different enzyome other than lysozyme and also various other molecules, like amino acids. Of course I need to first create the xyz data for those. I then check for optimization potential in the C extension by stringly orienting myself on all the tests (not only the test_against_vmd_reference_data, but also all the others). 
- Generalize the postprocessing: More adaptable graphs, less hardcoded  labels etc...

## Refactoring
- Think about a more object-oriented approach for better readability. Maybe profile the whole package to avoid performance issues. Thinks to consider here:
  - Use of `__slots__` in objects
  - Define interfaces using `abc`
  - Pre-allocate memory in collection-object, like so `l = [0] * 1000000`
  - Think about where methods and functions lie within objects: local vs. global
  - If I want to profile the whole thing, I should proably mock the lmp-executable to make the profiling more efficient
- Add type-hints an docstrings everywhere they are still missing
- Re-distribute the methods in `conversion.py`, `gro2lammps.py` and `helper.py` and maybe change the naming
- Create a submodule `helper` (or `misc` or similar) and move several files here
- IMO the `sasa_lammps/` should only contain `sasa_core.py` and `sasa_main.py` 
- Make the postprocessing optional. To do so construct a proper result object that, at least holds information about the result files, but at best has the results as objects. Then the postprocessing can be called on that object directly. This allows 1.) the postprocessing to be optional 2.) the later addition of user-configuration for the postprocessing