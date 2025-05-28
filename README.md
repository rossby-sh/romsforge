# romsforge (Unified Python + Fortran ROMS Preprocessing System)

> ⚠️ **Experimental Notice**  
> This project is in an **experimental stage** and intended for research and development purposes only.  
> Interfaces and functionalities are subject to change. 

- **Author**: Seonghyun Jo  
- **Contact**: shjo9536@gmail.com / birostris36@gmail.com

---

## Structure

### `py/`
- Python-based preprocessing tools for ROMS (forcing, ini, bry, obs)
- Real-time development in `dev_individual/`, integrated into `src_stable/`

### `ft/` *(study-only)*
- Fortran routines for remapping, NetCDF I/O, diagnostics

### `workflows/` *(planned)*
- Shell-based orchestration of Python and Fortran processes

### `logs/` *(planned)*
- Development history and debugging records

---

## Build *(planned)*

```bash
cd ft/
make

romsforge/
├── py/                        # Python modules
│   ├── src_stable/            # Stable pipelines
│   ├── src_experimental/      # In-progress integrated workflows
│   ├── dev_individual/        # Prototypes, isolated experiments
│   ├── roms_to_zstd/          # Postprocessing and analysis tools
│   └── env_settings/          # Common config and path settings

├── ft/                        # Fortran modules
│   ├── src_stable/            # Stable remapping, I/O routines
│   ├── src_experimental/      # Work-in-progress Fortran code
│   ├── utils/                 # Shared submodules
│   ├── test/                  # Unit tests
│   └── Makefile               # Build system

├── workflows/                 # Orchestration scripts
│   ├── run_preprocess.sh
│   ├── run_model.sh
│   └── run_postprocess.sh

├── logs/                      # Development logs
│   ├── py/
│   └── ft/

├── docs/                      # Diagrams, explanations
│   └── structure.md

└── README.md

```
