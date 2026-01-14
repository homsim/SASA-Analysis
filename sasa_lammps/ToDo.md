# ToDo

- Currently, the test_against_vmd_reference_data tests are failing for all parameters. I already checked: The newly added vdW-radii do not make a significant difference. This basically leaves only our new C extension implementation. My approach is to first generate more test data from the original vmd-python using the generate_vmd_reference.py. I could for example use a different enzyome other than lysozyme and also various other molecules, like amino acids. Of course I need to first create the xyz data for those. I then check for optimization potential in the C extension by stringly orienting myself on all the tests (not only the test_against_vmd_reference_data, but also all the others). 
- Generalize the postprocessing: More adaptable graphs, less hardcoded  labels etc...
- Provide a proper warning message in case the `lammps_exe` variable is provided but points to a non-existent file
