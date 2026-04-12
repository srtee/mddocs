# Constant Potential Molecular Dynamics in LAMMPS — Learning Taxonomy

## 0. Prerequisites: What You Need Before Starting

### 0.1 Programming & Scripting Foundations
- **Bash command-line familiarity** — navigating directories, running programs, reading log output
- **Text editor proficiency** — writing and editing input scripts (no IDE required; VS Code or Nano suffice)
- **Basic Python (optional but helpful)** — for post-processing and data analysis
- **Understanding of file paths** — absolute vs. relative, working directories in HPC environments

### 0.2 Computer Setup for LAMMPS
- **Installing LAMMPS** — pre-built binaries vs. building from source; CMake vs. make
- **Running LAMMPS** — command-line invocation, input scripts, output log files
- **Verifying your installation** — running a simple test case (e.g., in.lj from examples)

---

## 1. Molecular Dynamics Foundations

### 1.1 What Molecular Dynamics Does
- Atoms have positions, velocities, and forces; Newton's laws produce trajectories
- An interatomic potential (force field) computes forces from atom positions
- Time integration advances positions/velocities stepwise

### 1.2 Key LAMMPS Concepts
- **Units** — LAMMPS has configurable unit systems (real, metal, si, etc.); consistency matters
- **Atom styles** — what properties each particle carries (charge, mass, radius, etc.)
- **Pair styles** — the interatomic potential function (Lennard-Jones, Coulombic, etc.)
- **Bonds, angles, dihedrals** — bonded interactions (not always needed)
- **Fixes** — operations applied each timestep (integrators, constraints, thermostats, etc.)
- **Computes** — calculate quantities on the fly
- **Groups** — named subsets of atoms for targeted commands

### 1.3 A Simple MD Simulation in LAMMPS
- Initialization commands (`units`, `atom_style`, `boundary`, `region`, `create_box`, `create_atoms`)
- Force field setup (`pair_style`, `pair_coeff`, `bond_style`, etc.)
- Energy minimization before production
- Fixes for time integration (`fix nve`, `fix nvt`, etc.)
- `run` command and thermodynamic output (`thermo`, `thermo_style`)

### 1.4 Thermostats and Barostats
- NVE (microcanonical) — energy conserved, no thermostat
- NVT (canonical) — constant particle number, volume, temperature
- NPT (isothermal-isobaric) — constant particle number, pressure, temperature
- Why temperature/pressure control matters for meaningful simulations

---

## 2. Electrostatics in Molecular Simulations

### 2.1 Why Electrostatics Is Hard
- Coulombic interactions are long-range: \( V \propto 1/r \)
- Direct summation scales as \( N^2 \) — prohibitive for large systems
- Periodic boundary conditions complicate the math

### 2.2 Cutoff and Short-Range Interactions
- Truncating Coulomb at a cutoff saves computation but introduces error
- Real-space vs. k-space (reciprocal space) splitting
- When cutoff-based methods fail (high ionic strength, long-range correlations)

### 2.3 Ewald Summation
- The fundamental idea: split Coulomb into real-space, reciprocal-space, and self terms
- The reciprocal-space sum is solved via FFT — faster than direct sum for large N
- Accuracy controlled by Ewald parameter (`alpha`) and cutoff

### 2.4 PPPM (Particle-Particle Particle-Mesh)
- PPPM / PME is the practical descendant of Ewald: maps charges to a mesh, uses FFT
- Lower computational scaling than naive Ewald: \( N \log N \) instead of \( N^{3/2} \)
- K-space solvers in LAMMPS: `kspace_style ewald`, `kspace_style pppm`

### 2.5 Long-Range Solvers and Pair Styles
- `pair_style coul/long` pairs with `kspace_style ewald` or `pppm`
- Must use matching pair + kspace styles for consistency
- Energy and forces computed consistently across real and k-space contributions

---

## 3. The Physical Chemistry of Electrodes and Electrolytes

### 3.1 Conductors vs. Insulators in MD
- In a conductor, charges rearrange to equalize electrostatic potential
- In an insulator, charges are fixed; no redistribution within the material
- Electrode materials (metals) behave as conductors

### 3.2 The Electrode-Electrolyte Interface
- Ions in solution accumulate near charged surfaces (double layer)
- Applied potential changes the free energy of ion adsorption
- Structure of solvent (water) near electrode surfaces matters

### 3.3 Fixed Charge vs. Polarizable Models
- Fixed-charge force fields (e.g., AMBER, CHARMM) — charges do not change
- Polarizable force fields — allow charge redistribution (Drude oscillators, core-shell, etc.)
- Constant potential method — the electrode itself is treated as polarizable/conductive

### 3.4 What "Constant Potential" Means Physically
- The electrostatic potential on each electrode surface is held fixed
- Charges on electrode atoms fluctuate to maintain that potential
- This is distinct from constant charge (fixed q on each atom)

---

## 4. The Constant Potential Method (CPM)

### 4.1 Historical Origins
- Siepmann & Sprik (1995) — original constant potential approach
- Reed et al. (2007) — modern reformulation enabling practical MD
- The method minimizes electrostatic energy with respect to electrode charges under the constraint of fixed potential

### 4.2 Core Idea: Energy Minimization Under Constraint
- The total electrostatic energy depends on all charges (electrode + electrolyte)
- We want each electrode to be at a specified potential — not a specified charge
- Mathematically: minimize \( E(\mathbf{q}) \) subject to \( \mathbf{A}\mathbf{q} = \mathbf{V} \) (potential constraints)
- This produces induced surface charges, not predetermined ones

### 4.3 Gaussian Charge Smearing
- Point charges lead to singular matrix at short distances (Coulomb catastrophe)
- Smearing each electrode charge as a Gaussian distribution makes the matrix invertible
- The `eta` parameter controls the width of the Gaussian smearing
- Larger eta = narrower distribution = more like point charge; smaller eta = wider distribution

### 4.4 Capacitance Matrix
- The \( N \times N \) capacitance matrix \( \mathbf{C} \) relates electrode charges to potentials:
  \[
  \mathbf{Q} = \mathbf{Q}_{0V} + \mathbf{C} \cdot \mathbf{V}
  \]
- \( \mathbf{Q}_{0V} \) is the charge configuration at zero applied potential (influenced by electrolyte)
- The matrix is computed from electrode geometry and dielectric environment

### 4.5 Constant Potential vs. Constant Charge
- **Constant charge** — each atom has a fixed charge; potential varies freely
- **Constant potential** — each electrode group has a fixed potential; charges fluctuate
- CPM is essential for studying electrochemical interfaces where potential is the controlled variable

### 4.6 Related Methods
- **Charge equilibration (QEq)** — redistributes atomic charges within a single分子的
- **Drude oscillators** — adds extra particles to model electronic polarizability
- **Core-shell models** — another polarizability approach
- CPM is distinct: it controls boundary conditions on the electrostatic potential, not internal polarizability

---

## 5. The ELECTRODE Package in LAMMPS

### 5.1 Package Overview
- Authors: Ahrens-Iwers, Tee, Meissner (since LAMMPS 4May2022)
- Implements three constant-potential fix styles: `electrode/conp`, `electrode/conq`, `electrode/thermo`
- Requires KSPACE package and LAPACK/BLAS

### 5.2 Installing the ELECTRODE Package
- CMake build: `-D PKG_ELECTRODE=yes -D PKG_KSPACE=yes`
- May need `-D USE_INTERNAL_LINALG=yes` if LAPACK linking fails
- No traditional make build support (only CMake)

### 5.3 `fix electrode/conp` — Constant Potential
- Sets the **potential** on each electrode group; charges respond
- Syntax: `fix ID group electrode/conp potential eta`
- `potential` — electrode potential in volts (or equal-style variable for dynamic control)
- `eta` — Gaussian width parameter; smearing determines matrix invertibility
- `couple` keyword — add additional electrode groups

### 5.4 `fix electrode/conq` — Constant Charge
- Sets the **total charge** on each electrode; potentials respond
- Syntax: `fix ID group electrode/conq charge eta`
- Useful when you know the charge injection/extraction rather than the voltage
- The potentials computed can be monitored to understand the system

### 5.5 `fix electrode/thermo` — Thermopotentiostat
- Implements a thermodynamically consistent thermostat for electrode charge/potential
- Adds thermal fluctuations to both potential and charge with correct statistics
- Syntax: `fix ID group electrode/thermo potential eta temp T_v tau_v`
- For studying fluctuations and non-equilibrium dynamics

### 5.6 Key Fix Keywords
- `algo` — algorithm selection:
  - `mat_inv` — precompute capacitance matrix, fast but memory-intensive
  - `mat_cg` — precompute elastance matrix, use conjugate gradient
  - `cg` — no precomputation, on-the-fly CG each step
- `symm` — charge neutrality constraint across all electrodes (`on` or `off`)
- `ffield` — finite-field mode: use periodic z-direction and internal E-field instead of slab geometry
- `etypes` — type-based optimized neighbor lists (faster if electrode/electrolyte types don't overlap)
- `write_mat` / `write_inv` / `read_mat` / `read_inv` — save/load capacitance matrix for restart or reuse

### 5.7 Output from ELECTRODE Fixes
- **Global scalar** — energy added to system by the fix (negative of total electrode charge × potential)
- **Global vector** — current potential on each electrode (useful for `conq` and `thermo`)
- **Global array** — capacitance matrix rows, elastance matrix rows, and charge-at-0V values

---

## 6. K-space Solvers for ELECTRODE

### 6.1 Why Special K-space Styles Are Needed
- Standard `pppm` assumes charges are fixed point charges
- ELECTRODE needs to provide electrode-electrode interaction matrix and electrode-electrolyte interaction vector to the fix
- The electrode charges are not fixed — they are solved for each timestep

### 6.2 Available Styles
- `kspace_style ewald/electrode`
- `kspace_style pppm/electrode`
- `kspace_style pppm/electrode/intel` (accelerated variant)
- These behave like their standard counterparts but support the ELECTRODE fix interface

### 6.3 Kspace Modify Options for Electrode Systems
- Slab Ewald / 2D Ewald for non-periodic z-direction: `kspace_modify slab ew2d`
- Wire boundary conditions: `kspace_modify wire`
- PPPM `amat` option: `onestep` vs. `twostep` for elastance matrix calculation (memory vs. speed tradeoff)

---

## 7. Building a Constant Potential Simulation

### 7.1 System Design Decisions
- **Electrode geometry** — planar surfaces, nanostructured electrodes (nanopores, graphene), tip/substrate
- **Electrolyte** — water model (TIP4P, SPC/E), ionic liquid, aqueous salt solution
- **Boundary conditions** — non-periodic in z (slab geometry) for default mode; or periodic with `ffield on`
- **Box size** — large enough to avoid electrode-electrode interactions through periodic images

### 7.2 Initialization for Electrode Systems
```
units metal          # or real — must be consistent
atom_style charge    # charge is required
boundary p p f       # non-periodic z; f = fixed
region electrode_bot  # define electrode region
region electrolyte    # define electrolyte region
create_box            # with multiple sub-regions
create_atoms          # electrode atoms + electrolyte atoms
```
- Electrode atoms must be in their own group(s)
- Electrolyte atoms are the remaining particles

### 7.3 Force Field Setup
```
pair_style lj/cut/coul/long  # short-range Coulomb + LJ
pair_coeff * * 0.0            # no LJ for coulomb-only electrolyte
pair_coeff 1 1 0.01 3.4      # electrode LJ parameters
pair_coeff 1 2 0.01 3.4      # electrolyte-electrode LJ
kspace_style pppm/electrode 1e-5  # long-range solver
```
- Electrode atoms typically have large LJ epsilon to prevent penetration
- Electrode-electrolyte interactions carefully tuned

### 7.4 Applying the ELECTRODE Fix
```
fix fxconp bot electrode/conp -1.0 1.805 couple top 1.0 symm on
```
- Bottom electrode at -1.0 V; top electrode at +1.0 V
- `symm on` constrains total electrode charge to zero
- The fix manages charge equilibration each timestep

### 7.5 Adding Thermostats
- Electrode atoms should typically be thermostatted differently from electrolyte
- Example: `fix frozen all nve` for frozen electrodes
- Or use `fix langevin` on electrolyte only, `fix nve` on frozen electrodes
- `fix electrode/thermo` includes built-in thermostatting

### 7.6 A Minimal Working Input Script
```
units metal
atom_style charge
boundary p p f

region box block 0 20 0 20 0 60
create_box 2 box

# electrode atoms
region bot block 0 20 0 20 0 5
region top block 0 20 0 20 55 60
create_atoms 1 region bot
create_atoms 2 region top

# electrolyte
region elyte block 0 20 0 20 5 55
create_atoms 3 region elyte

group bot region bot
group top region top
group elyte region elyte

pair_style lj/cut/coul/long 10.0
pair_coeff 1 1 0.01 3.0
pair_coeff 2 2 0.01 3.0
pair_coeff 1 2 0.005 3.0
pair_coeff 1 3 0.005 3.0
pair_coeff 2 3 0.005 3.0
pair_coeff 3 3 0.15 3.2
kspace_style pppm/electrode 1e-5

fix fxconp bot electrode/conp -1.0 1.805 couple top 1.0 symm on

fix frozen all nve
fix mytemp elyte langevin 300 300 100 12345

run 10000
```

---

## 8. Practical Considerations and Common Issues

### 8.1 Convergence
- Electrode charges must converge each timestep
- Tolerance in the fix (`algo mat_inv` tolerance is implicit in precomputation)
- If charges don't equilibrate: increase `maxiter`, check eta value, check geometry

### 8.2 Electrode Immobilization
- Matrix-based algorithms (`mat_inv`, `mat_cg`) require electrode positions to be fixed
- `algo cg` allows moving electrodes but is slower
- For dynamic electrodes (charging/discharging), consider `ffield on` with `algo cg`

### 8.3 Potential Control and Variables
- Any potential/charge parameter can be an equal-style variable
- Ramp potentials over simulation time: `variable V equal ramp(0.0, 2.0)`
- Use `fix modify` with `tf` option for Thomas-Fermi metallicity model (quantum corrections for real metals)

### 8.4 Parallel Performance
- Electrode atoms should be evenly distributed across processors
- `processors * * 2` maps 2D decomposition for slab geometries
- Matrix storage can exceed 0.5 GiB per MPI process — use `mat_cg` or `algo cg` for large systems

### 8.5 Restarting and Reproducibility
- No restart data written yet (as of current LAMMPS version)
- Save/restore capacitance matrix with `write_inv` / `read_inv`
- Always set random seed explicitly for thermostat (`rng_v` in `electrode/thermo`)

### 8.6 Units Gotchas
- Potentials are always in volts regardless of `units` setting
- eta is in inverse length units; check that the value is appropriate for your unit system
- Charges in `electrode/conq` must be in the same units as the rest of the simulation

---

## 9. Example Systems to Study

### 9.1 Planar Electrode with Aqueous Electrolyte
- Classic double-layer capacitor setup
- Two planar metal electrodes, water + salt in between
- Measure capacitance, ion density profiles, potential distribution

### 9.2 Nanostructured Electrodes
- Graphene sheets, nanopores, electrode tips
- Non-planar geometry is where CPM shines vs. fixed-charge models
- Study charging dynamics, ion selectivity, field enhancement

### 9.3 Graphite-Ionic Liquid Interface
- Room-temperature ionic liquids as electrolytes
- High potential differences; fixed-charge models fail at these conditions
- Study energy storage mechanisms, electrodecreening

### 9.4 Electrode with Finite Field (`ffield on`)
- Periodic in z-direction; potential difference applied via internal E-field
- Allows smaller box in z; better for some geometries
- Requires `symm on`; electrode charge neutrality is imposed

---

## 10. Further Reading

### 10.1 Foundational Papers
- Siepmann & Sprik, J. Chem. Phys. 102, 511 (1995) — original CPM
- Reed et al., J. Chem. Phys. 126, 084704 (2007) — modern reformulation
- Ahrens-Iwers & Meissner, J. Chem. Phys. 155, 104104 (2021) — LAMMPS ELECTRODE package
- Ahrens-Iwers et al., J. Chem. Phys. 157, 084801 (2022) — ELECTRODE package validation

### 10.2 Related Methods Papers
- Gingrich, MSc thesis (2010) — Gaussian smearing for CPM
- Deissenbeck et al., Phys. Rev. Letters 126, 136803 (2021) — thermopotentiostat
- Dufils et al., Phys. Rev. Letters 123, 195501 (2019) — finite-field CPM
- Scalfi et al., J. Chem. Phys. 153, 174704 (2020) — Thomas-Fermi model for metals

### 10.3 LAMMPS Documentation
- [ELECTRODE package details](https://docs.lammps.org/Packages_details.html#electrode-package)
- [fix electrode/conp](https://docs.lammps.org/fix_electrode.html)
- [kspace_style](https://docs.lammps.org/kspace_style.html) — electrode variants
- [LAMMPS examples/PACKAGES/electrode](https://github.com/lammps/lammps/tree/develop/examples/PACKAGES/electrode)

---

## Scope Notes

This outline assumes a typical undergraduate in chemistry, physics, or engineering who has completed:
- Introductory physics (mechanics, electromagnetism basics)
- One semester of general chemistry (understanding of ions, Coulomb's law, electrostatic energy)
- Basic programming experience (any language)

The progression moves from **what MD does** → **how electrostatics is handled computationally** → **why electrode surfaces need special treatment** → **the CPM method** → **LAMMPS implementation specifics** → **practical simulation design**.

Articles at the leaf level of this taxonomy should be written as **standalone explainers** with:
- A conceptual introduction (why this concept matters)
- Concrete, working example code where applicable
- Expected output / interpretation guidance
- Common pitfalls and how to troubleshoot them
