# ToDo

- Generalize the postprocessing: More adaptable graphs, less hardcoded  labels etc... (this may not be my decision though)
  - Go over the necessary files and check if an invocation of `postprocessing` is possible without a previous `compute` (in case someone wants to do postprocessing without re-computing the whole thing). 
  - Throw a proper exception in case files are missing for the postprocessing
- Adapt a unified coding style (especially in the namings)
- Go over all the tests
  - Refactor the fixtures and share paths, the resource-path in particular
  - Sort out unnecessary tests and maybe add new ones 