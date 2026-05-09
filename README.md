# Microwave Cathode Full-Fluid Plasma Code

This repository contains a full-fluid plasma simulation code for a microwave cathode based on a non-neutral drift-diffusion (NDD) model.

The code is intended for numerical analysis of plasma transport, electrostatic potential, microwave energy coupling, and related plasma-fluid behavior in a microwave cathode system.

## Overview

This code solves a full-fluid plasma model for a microwave cathode.

The model is based on a non-neutral drift-diffusion formulation and includes fluid equations for charged and neutral species, electrostatic field calculation, and microwave-related source terms.

The current implementation is written in C++ and uses external numerical libraries for linear solvers.

## Features

- Full-fluid plasma simulation for a microwave cathode
- Non-neutral drift-diffusion (NDD) model
- Ion, electron, and neutral-species transport modeling
- Electrostatic potential calculation
- Microwave power coupling parameters
- Secondary electron emission settings
- Restart file support
- YAML-based setup file
- Output data generation for post-processing

## Directory Structure

```text
.
├── README.md
├── Makefile
├── include/
│   └── Header files
├── src/
│   └── Source files
├── input/
│   └── setup.yaml
├── restart/
│   └── Restart files
├── results/
│   └── Output files
└── scripts/
    └── Post-processing or utility scripts
```

## Requirements

The code requires the following libraries:

- C++ compiler with C++17 support
- MPI
- HYPRE
- yaml-cpp

On macOS with Homebrew, these libraries may be installed using:

```bash
brew install open-mpi hypre yaml-cpp
```

## Build

Compile the code using:

```bash
make
```

or, depending on the Makefile target:

```bash
make Normal
```

The executable file is generated as:

```text
a.out
```

TODO: Update this section if the executable name is changed, for example to `build/solver`.

## Run

Run the simulation using:

```bash
./a.out
```

The setup file is read from:

```text
input/setup.yaml
```

Before running the code, make sure that `input/setup.yaml` exists and contains all required parameters.

## Setup File

The simulation parameters are defined in:

```text
input/setup.yaml
```

The setup file is written in YAML format. Each parameter is grouped by category, while the parameter names are kept identical to the corresponding C++ variable names.

Example:

```yaml
microwave:
  Pmw: 8.0                # (double) Microwave power [W]
  omegam: 26703537556e0   # (double) Microwave angular frequency

scheme:
  icon_PC: 1              # (int) Whether to use the PC method
  dt_ini: 1.0e-10         # (double) Initial time step
```

The basic rule is:

```text
category:
  C++_variable_name: value
```

The category is used only for readability and organization. The parameter name itself should correspond one-to-one with the variable name used in the C++ code.

## Main Parameter Categories

The setup file contains the following categories:

- `microwave`: microwave power, frequency, and related parameters
- `plasma`: plasma and neutral particle parameters
- `material`: dielectric material properties
- `bias`: bias voltage settings
- `neutral_inlet`: neutral gas inlet conditions
- `transport_model`: transport model switches and coefficients
- `see`: secondary electron emission parameters
- `scheme`: numerical scheme and time-step settings
- `solver`: solver convergence criteria and iteration settings
- `microwave_coupling`: microwave calculation interval
- `relaxation`: under-relaxation factors
- `simulation`: global simulation settings
- `output`: output settings

## Output

Simulation results are written to the output directory defined by the code and setup parameters.

Typical output files may include:

- plasma density
- electric potential
- electric field
- velocity or flux data
- convergence information
- restart data

TODO: Add the exact output file names and their meanings.

## Restart Files

Restart files are used to resume simulations from previously saved states.

The `restart/` directory is intended for restart data. Since restart files can become large, they are usually not tracked by Git.

Recommended `.gitignore` setting:

```gitignore
restart/*
!restart/.gitkeep
```

## Results and Generated Files

Generated result files should usually not be committed to GitHub.

Recommended `.gitignore` setting:

```gitignore
results/*
!results/.gitkeep

*.o
*.out
a.out
build/
```

Small sample input files may be tracked, but large simulation outputs should be excluded from Git.

## Post-processing

Post-processing scripts can be placed in:

```text
scripts/
```

For example:

```text
scripts/
├── plot_result.py
└── plot_contour.py
```

TODO: Add instructions for visualization and analysis scripts.

## Notes for Developers

When adding new setup parameters:

1. Add the parameter to `input/setup.yaml`.
2. Use the same name as the corresponding C++ variable.
3. Add the parameter-reading line in `inputParam()`.
4. Specify the type in the YAML comment, such as `(int)` or `(double)`.
5. Check that the parameter is correctly loaded at runtime.

Example:

```yaml
scheme:
  new_parameter: 1.0   # (double) Description of the new parameter
```

Corresponding C++ side:

```cpp
category = "scheme";
READ_DOUBLE(category, new_parameter);
```

## Version Control Policy

The following files should be tracked by Git:

- source files in `src/`
- header files in `include/`
- setup files in `input/`
- Makefile
- README.md
- small sample input files
- utility scripts

The following files should usually not be tracked:

- executable files
- object files
- build directories
- large result files
- restart files
- temporary files

## License

TODO: Add license information.

Commercial use is permitted if the selected license allows it.

Recommended common options:

- MIT License
- BSD 3-Clause License
- Apache License 2.0

For a simple permissive license, the MIT License is commonly used.

## Citation

TODO: Add citation information.

If this code is used in academic work, please cite the related paper, thesis, or repository.

Example:

```text
Ryo Shirakawa, "TODO: Title of related work", TODO: Journal/Conference/Thesis, TODO: Year.
```

## Author

Ryo Shirakawa

TODO: Add affiliation, laboratory, or contact information if needed.

## Acknowledgements

TODO: Add acknowledgements for supervisors, collaborators, funding, or related software libraries.
