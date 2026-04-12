# Constant Potential Molecular Dynamics in LAMMPS — Learning Taxonomy

## 0. Orientation: What You Need to Know

### 0.1 Prior Knowledge
- Introductory mechanics and electromagnetism (Coulomb's law, electrostatic energy)
- Basic chemistry (ions, water structure, what an electrode is)
- One programming language (any) — for post-processing
- No prior MD experience required

### 0.2 Your Starting Point
- You have LAMMPS installed (or access to an HPC cluster with it)
- You can run an input script and see log output
- You want to understand what the commands *mean*, not just what they *do*

### 0.3 How to Use This Guide
- The goal is to connect LAMMPS input lines → physical meaning → simulation output
- Theory sections explain the physics; Implementation sections show the corresponding LAMMPS syntax
- Every major concept links to an observable you can measure

---

## 1. Molecular Dynamics: What the Simulation Actually Does

### 1.1 The Core Loop
- Newton determines atomic motion; forces come from interatomic potentials
- A timestep Δt advances positions and velocities: `x(t+Δt) = x(t) + v(t)Δt`
- Forces are recalculated at every step — this is what makes MD computationally expensive
- The trajectory is a deterministic sequence; statistics over many steps give thermodynamic properties

### 1.2 Key LAMMPS Concepts (Syntax → Meaning)
- **Units** — define the physical scale (real = Å, kcal/mol, fs; metal = Å, eV, ps). Mismatch causes garbage.
- **Atom style** — what properties each particle carries. `charge` is mandatory for electrostatics.
- **Boundary** — `p p f` means periodic in x,y, fixed (non-periodic) in z. Determines slab geometry for electrode simulations.
- **Pair styles** — the functional form of the interatomic potential. Determines what forces exist.
- **Fixes** — operations applied every timestep (integration, thermostats, constraints). The verb of LAMMPS.
- **Computes** — calculate quantities on the fly (temperature, stress, charge density)
- **Groups** — named subsets of atoms for targeted commands. Essential for distinguishing electrode from electrolyte.

### 1.3 Reading LAMMPS Output
- `Thermodynamic output` — energy (total, kinetic, potential), temperature, pressure, volume
- Energy conservation indicates a well-behaved simulation (NVE without thermostat)
- Potential energy components (pair, bond, kspace) reveal where computation time goes
- Drift in conserved quantity signals integration problems

---

## 2. Force Fields — What Interactions Exist and How to Configure Them

### 2.1 The Role of a Force Field
- A force field is a set of mathematical functions that describe interatomic forces
- It is an approximation — no force field is "correct", only appropriate for certain systems
- Choice of force field determines what physics is captured (and what is omitted)
- All force field parameters are fitted to experimental data or quantum calculations

### 2.2 Non-Bonded Interactions

#### 2.2.1 Van der Waals (Dispersion/Repulsion)
- **Lennard-Jones (LJ)**: `pair_style lj/cut` — the workhorse of MD
  - Formula: `V(r) = 4ε[(σ/r)^12 - (σ/r)^6]`
  - ε = depth of the energy well (binding strength); σ = distance at which energy = zero
  - The r⁻¹² term is steep repulsion at short range; r⁻⁶ is softer attraction at longer range
  - Physical meaning: Pauli repulsion at short distance, induced-dipole attraction at longer range
  - Cutoff: LJ is truncated at some `rcut` — introduces error; `rcut` typically 10–12 Å for water
  - **Why cut?** Computation scales as N²; truncation saves time with acceptable error
- **Other forms**: Buckingham, Morse, 12-6-4 (for metals), Gaussian — each suited to different materials

#### 2.2.2 Electrostatic (Coulombic)
- **Formula**: `V(r) = k q₁q₂/r` — long-range, never truly zero
- **The problem**: Direct summation scales as N²; prohibitive for large systems
- **Solutions**: Ewald summation, PPPM (Particle-Particle Particle-Mesh) — see Section 3

#### 2.2.3 Combined Short-Range Potentials
- `pair_style lj/cut/coul/long` — LJ plus real-space Coulomb; k-space handles the long-range Coulomb
- This is the standard choice for electrolyte simulations: LJ for repulsion/attraction, Coulomb for ionic interactions

### 2.3 Bonded Interactions (When They Matter)
- **Bonds**: `bond_style harmonic` — spring-like; `V = k(r - r₀)²`
  - Water models (TIP3P, TIP4P) have explicit O-H bonds
  - Electrolyte ions may have internal bonds (e.g., carbonate)
- **Angles**: angle between three atoms; important for water orientation
- **Dihedrals**: rotations around bonds; relevant for organic electrolytes

### 2.4 Cross-Terms and Combining Rules
- **Mixing rules**: Lorentz-Berthelot (`ε_ij = √(ε_i ε_j)`, `σ_ij = (σ_i + σ_j)/2`)
- **pair_coeff** syntax: `ID1 ID2 epsilon sigma` for each unique pair type
- **All pairs must be defined**: LAMMPS will error if any pair type is unspecified

### 2.5 Force Field Parameters — Where They Come From
- Water models: TIP3P, TIP4P, SPC/E — each is a compromise between accuracy and cost
- Ion parameters: from literature (e.g., Joung & Cheatham, 2008)
- Electrode atoms: typically large ε for hard-wall repulsion, no Coulomb (metallic screening)
- **Physical meaning of LJ on electrodes**: prevents electrolyte atoms from penetrating the surface; does not model bonding

### 2.6 Force Field → Output Connection
| Input | Output Effect |
|-------|---------------|
| `pair_coeff * * ε σ` | Potential energy scale; particle density; diffusion coefficient |
| `dielectric` | Strength of Coulomb interactions (vacuum = 1, water ≈ 80) |
| `bond_style` + `bond_coeff` | Vibrational spectra; bulk modulus; compressibility |
| Cutoff distance | Energy drift; simulation runtime (larger cut = slower but more accurate) |

---

## 3. Physical Conditions — How the Simulation Enforces Temperature, Pressure, and Time

### 3.1 What Is Being Held Constant?
- Real experiments control T, P (or V), and composition
- In MD, we control what we simulate: NVE (microcanonical), NVT (canonical), NPT (isothermal-isobaric)
- The choice is not just convention — it determines which statistical ensemble produces your output

### 3.2 Time Integration and Timestep

#### 3.2.1 Newton's Equations and the Timestep
- Position update: `x(t+Δt) = 2x(t) - x(t-Δt) + (F/m)Δt²` (Verlet integration)
- **Velocity Verlet**: `v(t) = v(t-Δt) + 0.5(F(t) + F(t+Δt))Δt` — more stable
- `Δt` must be small enough to resolve the fastest motion (typically C-H bond vibration)
- **Typical values**: `Δt = 1 fs` for explicit water; `Δt = 0.5 fs` if bonds are flexible; `Δt = 2 fs` with bond constraints

#### 3.2.2 Why Timestep Matters for Output
- Too large Δt → energy drift, simulation "blows up" (atoms overlap catastrophically)
- Too small Δt → wasted computation; results unchanged
- **Constraint**: `fix shake` or `fix rattle` removes fast degrees of freedom, allowing larger Δt

#### 3.2.3 Verlet → Input Connection
```
fix 1 all nve          # basic NVE integration
fix 1 all nvt          # NVE + thermostat (see below)
```
- `nve` = no thermostat; energy should be conserved (check `E_total` drift in log)
- `nvt` = NVE modified by thermostat; total energy is NOT conserved

### 3.3 Thermostats — Controlling Temperature

#### 3.3.1 Temperature in MD
- Temperature is a statistical property: `T = (2/3Nk) × (1/2 mv²)_avg`
- Instantaneous T fluctuates; we average over time to compare to experiment
- Thermostats manipulate velocities to sample the canonical ensemble

#### 3.3.2 Thermostat Algorithms (Meaning → LAMMPS)

| Method | LAMMPS Fix | How It Works | Output Implication |
|--------|-----------|--------------|-------------------|
| Velocity-rescale | `fix temp/rescale` | Scale velocities to match target T | Simple but not rigorous; good for equilibration |
| Berendsen | `fix berendsen` | Weak coupling to heat bath | Produces incorrect fluctuations |
| Nosé-Hoover | `fix nvt` or `fix npt` | Extended Lagrangian; ergodic sampling | Standard for production; `temp/press` output reflects ensemble |
| Langevin | `fix langevin` | Stochastic damping + random forces | Good for non-equilibrium; adds noise to dynamics |
| Canonical (Nosé-Hoover chain) | `fix nvt` with `tchain` | Multiple thermostats for better ergodicity | More expensive but better for some systems |

#### 3.3.3 What Temperature Controls Physically
- Kinetic energy → particle velocities → diffusion coefficient
- Too high T → bond breaking, chemical reactions (LAMMPS not suited for this)
- Too low T → slow dynamics, metastable states

### 3.4 Barostats — Controlling Pressure

#### 3.4.1 Pressure in MD
- Pressure from virial theorem: `P = (NkT/V) + (1/3V)⟨Σ r·F⟩`
- Mechanical pressure from forces + kinetic contribution
- NPT adjusts box size to achieve target P; NPH holds P fixed without T control

#### 3.4.2 Barostat Algorithms (Meaning → LAMMPS)

| Method | LAMMPS Fix | Output Implication |
|--------|-----------|-------------------|
| Nose-Hoover | `fix npt` | P fluctuates; volume adjusts to target; good for density equilibration |
| Berendsen | `fix barostat` | Weak coupling; incorrect fluctuations |
| MTK | `fix npt` with `iso` or `aniso` | Martyna-Tobias-Klein — more reliable for anisotropic stress |
| Bussi | `fix npt` with `temp` Bussi | Stochastic barostat; better for non-equilibrium |

#### 3.4.3 Barostat → Physical Output
- **Volume** → density (g/cm³) — compare to experiment
- **Stress tensor** → mechanical properties; anisotropic systems (like electrode interfaces) need careful choice
- Isotropic (`iso`) vs. anisotropic (`aniso`) — slab geometry uses `aniso` or `semiiso`

### 3.5 Integrators in LAMMPS (Summary)

| Fix | Ensemble | Use Case |
|-----|----------|----------|
| `fix nve` | NVE | Microcanonical; equilibration check; frozen atoms |
| `fix nvt` | NVT | Canonical ensemble; most production runs |
| `fix npt` | NPT | Isothermal-isobaric; density-sensitive properties |
| `fix nph` | NPH | Isobaric without thermostat; rarely used alone |

### 3.6 Observables — What You Actually Measure

#### 3.6.1 Thermodynamic Output
- `thermo` N — print every N steps
- `thermo_style custom` — select variables: `temp`, `press`, `pe`, `ke`, `vol`, `density`

#### 3.6.2 Computed Quantities
- `compute temp all temp` — kinetic temperature (available for groups)
- `compute stress/atom` — per-atom stress tensor (local pressure)
- `compute charge` — per-atom charge (if needed for analysis)

#### 3.6.3 Fix-Specific Output
- `fix electrode/conp` outputs: electrode potentials (global vector), energy (scalar), capacitance matrix rows (array)
- Use `fix ID group_name electrode/conp ...` then `variable Vx equal f_ID[2]` to access electrode potential

#### 3.6.4 Trajectory Files
- `dump 1 all custom 100 traj.xyz id type x y z q` — write trajectory
- `dump 2 all dcd 100 traj.dcd` — binary format; smaller files
- Post-processing: Python (MDAnalysis, pytim), VMD, Ovito

---

## 4. Electrostatics in Molecular Simulations — How LAMMPS Handles Long-Range Interactions

### 4.1 Why Electrostatics Is the Hard Part
- Coulomb: `V ∝ 1/r` — long-range even at moderate distances
- Direct summation: N(N-1)/2 pairs → computationally prohibitive for N > 1000
- Periodic boundaries make the math non-trivial (Ewald, slab corrections)

### 4.2 Ewald Summation — The Conceptual Basis

#### 4.2.1 The Split
- Real space: short-range Coulomb interactions summed directly (cutoff-based)
- Reciprocal space: long-range interactions handled via Fourier transform
- Self term: corrects for double-counting of interactions in periodic systems

#### 4.2.2 Ewald Parameter (α)
- Controls the split between real and reciprocal space
- Large α → narrow real-space Gaussian → faster real-space convergence → slower k-space
- Small α → opposite
- LAMMPS default is usually reasonable; for charged systems with small box, tuning may help

### 4.3 PPPM (Particle-Particle Particle-Mesh)

#### 4.3.1 Why PPPM Is Used Instead of Naive Ewald
- Ewald: O(N³/²) scaling — still bad for very large N
- PPPM: maps charge to mesh, uses FFT → O(N log N)
- Standard for systems with >10,000 atoms

#### 4.3.2 PPPM → LAMMPS Syntax
```
kspace_style pppm 1e-5
```
- `1e-5` = desired RMS force accuracy (relative to force magnitude)
- Tighter tolerance = more accurate but slower
- For electrode systems: `kspace_style pppm/electrode 1e-5`

### 4.4 Cutoff-Based Methods (When PPPM Is Unnecessary)
- Small systems (N < 1000) or very short simulations: direct Coulomb with cutoff may suffice
- `pair_style coul/cut 10.0` — cutoff at 10 Å; no k-space calculation
- **Warning**: ignoring long-range electrostatics introduces serious error; do not use for production

### 4.5 Boundary Conditions and Electrostatics

#### 4.5.1 Slab Geometry (Non-Periodic Z)
- `boundary p p f` — z is non-periodic
- Electrostatics requires 2D Ewald correction: `kspace_modify slab 3`
- Without this, forces in z-direction are wrong (artificial dipole-dipole interactions)

#### 4.5.2 Periodic Z with Finite Field
- `boundary p p p` with `ffield yes` in ELECTRODE fix
- Internal E-field applied; simulates infinite stack of identical cells
- No slab correction needed; better for some geometries

### 4.6 Pair Style — Kspace Consistency
- `pair_style lj/cut/coul/long` MUST pair with `kspace_style pppm` (or ewald)
- `pair_style lj/cut/coul/cut` pairs with NO kspace (direct Coulomb only)
- Mismatch causes LAMMPS to error or produce wrong forces

---

## 5. Conductors, Electrodes, and the Constant Potential Method

### 5.1 Conductors vs. Insulators in MD
- **Insulator**: charges are fixed; electrostatic potential varies freely
- **Conductor**: charges redistribute to equalize potential within the material
- **Electrode (metal)**: behaves as conductor; charges respond to applied voltage

### 5.2 Fixed Charge vs. Polarizable Models

#### 5.2.1 Fixed Charge (Standard Force Fields)
- AMBER, CHARMM, OPLS — atomic charges set at start, never change
- Appropriate for insulators and vacuum
- **Problem for electrodes**: charge cannot respond to applied potential

#### 5.2.2 Polarizable Models
- **Drude oscillators**: extra particle attached by harmonic spring; models electron cloud response
- **Core-shell**: massless shell encloses ion; responds to local field
- **QEq (charge equilibration)**: redistribute charges within a single molecule
- **CPM (Constant Potential Method)**: controls boundary condition on electrostatic potential

### 5.3 The Electrode-Electrolyte Interface

#### 5.3.1 Structure
- Ions accumulate near charged surfaces → electric double layer
- Water molecules orient near electrode (hydrogen up/down)
- Capacitance depends on ion size, hydration, electrode structure

#### 5.3.2 Applied Potential
- Real experiments: potentiostat holds electrode at fixed voltage
- Fixed-charge MD: cannot enforce voltage; can only set atom charges (wrong)
- CPM: enforces voltage directly → correct electrochemical interface

### 5.4 The Constant Potential Method — Theory

#### 5.4.1 The Problem
- We want: electrostatic potential on electrode surface = V (fixed)
- We have: electrode charges, electrolyte charges, geometry
- Energy: `E = (1/2) Σ q_i φ_i` — depends on all charges
- Constraint: electrode potentials must equal specified values

#### 5.4.2 Energy Minimization Under Constraint
- Minimize `E(q)` subject to `Aq = V` (potential constraints)
- Solution: induced surface charges that depend on geometry and electrolyte
- These are NOT predetermined — they emerge from the minimization

#### 5.4.3 Gaussian Charge Smearing
- Point charges at surface → singular matrix (diverges as r → 0)
- Smear each charge as a Gaussian distribution (width η)
- `η` controls smearing: larger η = narrower (more point-like); smaller η = wider
- Physical meaning: metals have finite electron density at surface; smearing mimics this

#### 5.4.4 The Capacitance Matrix
- N×N matrix relates electrode potentials to charges:
  ```
  Q = Q₀ᵥ + C · V
  ```
- `Q₀ᵥ` = charge at zero applied potential (influenced by electrolyte structure)
- `C` = capacitance matrix (depends on electrode geometry and dielectric)
- Physical meaning: Cij = how much charge on electrode i changes when electrode j's potential changes

### 5.5 Constant Potential vs. Constant Charge

| Aspect | Constant Potential (electrode/conp) | Constant Charge (electrode/conq) |
|--------|-------------------------------------|----------------------------------|
| What is fixed | Electrode voltage | Total electrode charge |
| What varies | Electrode atom charges | Electrode potentials |
| Use case | Experiment-like (potentiostat) | Known charge injection |
| Output | Fluctuating charges | Fluctuating potentials |

- CPM is essential for electrochemistry where potential is the controlled variable

---

## 6. The ELECTRODE Package — LAMMPS Implementation

### 6.1 Package Overview
- Authors: Ahrens-Iwers, Tee, Meissner (LAMMPS 4May2022+)
- Implements constant potential, constant charge, and thermopotentiostat methods
- Requires: KSPACE package, LAPACK/BLAS

### 6.2 Installing the ELECTRODE Package
- CMake: `-D PKG_ELECTRODE=yes -D PKG_KSPACE=yes`
- If LAPACK linking fails: `-D USE_INTERNAL_LINALG=yes`
- No make-based build (CMake only)

### 6.3 The Three Fix Styles

#### 6.3.1 `fix electrode/conp` — Constant Potential
```
fix fxconp bot electrode/conp -1.0 1.805 couple top 1.0 symm on
```
- **Purpose**: hold electrode groups at specified potentials; charges respond
- **Arguments**: `potential eta` (potential in volts; eta in inverse length units)
- **Key options**:
  - `couple groupName V` — add additional electrode group at potential V
  - `symm on/off` — enforce total electrode charge neutrality
  - `algo mat_inv|mat_cg|cg` — algorithm choice (see below)
  - `ffield on/off` — finite-field mode for periodic z
  - `etypes` — optimize neighbor lists when electrode/electrolyte types don't overlap

#### 6.3.2 `fix electrode/conq` — Constant Charge
```
fix fxconq bot electrode/conq -0.5 1.805 couple top 0.5 symm on
```
- **Purpose**: set total charge on electrode groups; potentials respond
- **Physical meaning**: model charge injection/extraction (battery discharge)
- **Output**: potentials on electrodes (can be monitored to understand response)

#### 6.3.3 `fix electrode/thermo` — Thermopotentiostat
```
fix fxthermo bot electrode/thermo -1.0 1.805 300 12345 50.0
```
- **Purpose**: adds thermal fluctuations to electrode potential (and charge)
- **Arguments**: `potential eta temp seed tau_v`
- **Physical meaning**: models Johnson-Nyquist noise in electrochemical systems
- **Use case**: non-equilibrium dynamics, fluctuation-dissipation studies

### 6.4 Algorithm Choices (`algo`)

| Algorithm | Precomputation | Memory | Speed | Use Case |
|-----------|----------------|--------|-------|----------|
| `mat_inv` | Capacitance matrix | High (N²) | Fast | < 10,000 electrode atoms; static electrodes |
| `mat_cg` | Elastance matrix | Medium | Medium | Medium systems; limited electrode motion |
| `cg` | None | Low | Slow | Dynamic electrodes; large systems |

- `mat_inv`: invert matrix once; fast on-the-fly charge calculation
- `mat_cg`: solve linear system with conjugate gradient; less memory than `mat_inv`
- `cg`: no precomputation; solve full problem each step; needed for moving electrodes

### 6.5 K-space Compatibility
- Use electrode-specific kspace styles:
  - `kspace_style ewald/electrode`
  - `kspace_style pppm/electrode`
  - `kspace_style pppm/electrode/intel` (optimized)
- These provide electrode-electrode interaction matrices to the fix

### 6.6 Slab Geometry Modifications
```
kspace_modify slab 3        # 2D Ewald correction (default for non-periodic z)
kspace_modify wire          # for wire boundary conditions
kspace_modify one-step      # for faster but less accurate elastance matrix
```

### 6.7 Matrix Save/Load
- `write_inv` / `read_inv` — save elastance matrix for restart
- Avoids recomputation for large systems
- No restart file support yet (as of current LAMMPS); matrix save is the workaround

### 6.8 Output from ELECTRODE Fixes

| Output | Access | Meaning |
|--------|--------|---------|
| Energy (scalar) | `f_ID` | `-Σ q_i V_i` — work done by electrodes on system |
| Electrode potentials (vector) | `f_ID[2], f_ID[3]...` | Potential on each electrode group |
| Capacitance matrix rows (array) | `f_ID[4+]` | C_ij rows; useful for analysis |

---

## 7. Building a Constant Potential Simulation — Complete Example

### 7.1 Design Decisions Before Writing Input

#### 7.1.1 Electrode Geometry
- Planar surfaces: simplest; good for double-layer studies
- Nanostructured: graphene, nanopores, tips — CPM shines here (charge responds to complex shape)
- Cell size: must be large enough that periodic images don't interact through vacuum

#### 7.1.2 Electrolyte Choice
- Water: SPC/E or TIP4P (not TIP3P — poor dielectric)
- Salt: NaCl, KCl — matched ion parameters to water model
- Ionic liquids: heavy; slower dynamics but interesting electrostatics

#### 7.1.3 Boundary Conditions
- Slab (`p p f`): standard for open surfaces; needs `kspace_modify slab`
- Periodic with `ffield`: smaller z-possible; requires `symm on`

### 7.2 Initialization Commands (Meaning)

```lammps
units metal                # eV, Å, ps — match your force field
atom_style charge         # charge needed; no bonds for simple electrolytes
boundary p p f            # non-periodic z for slab geometry

region box block 0 20 0 20 0 60
create_box 2 box          # 2 atom types (bottom, top electrode)

region bot block 0 20 0 20 0 5
region top block 0 20 0 20 55 60
create_atoms 1 region bot
create_atoms 2 region top
```

- Electrode atoms get distinct types (1, 2) to distinguish groups
- `group bot region bot` — defines the group for the electrode fix

### 7.3 Force Field Setup

```lammps
pair_style lj/cut/coul/long 10.0    # LJ + real-space Coulomb; cutoff 10 Å
pair_coeff 1 1 0.01 3.0              # electrode-electrode LJ (soft repulsion)
pair_coeff 2 2 0.01 3.0              # top electrode
pair_coeff 1 2 0.005 3.0            # cross-interaction (smaller)
pair_coeff 3 3 0.15 3.2              # water LJ (from water model)
pair_coeff 1 3 0.005 3.0            # electrode-water LJ
pair_coeff 2 3 0.005 3.0

kspace_style pppm/electrode 1e-5    # long-range solver; electrode-aware
```

- Electrode LJ ε is small (0.01 eV) — creates soft repulsion, not bonding
- Water parameters depend on the water model used (SPC/E, TIP4P, etc.)
- kspace tolerance `1e-5` is standard; `1e-6` for production

### 7.4 Applying the ELECTRODE Fix

```lammps
fix fxconp bot electrode/conp -1.0 1.805 couple top 1.0 symm on
```

- **Bottom electrode**: -1.0 V
- **Top electrode**: +1.0 V (via `couple`)
- **symm on**: total electrode charge = 0 (neutral system)
- **eta = 1.805 Å⁻¹**: standard value (check literature for your system)

### 7.5 Thermostating Strategy

```lammps
fix frozen all nve                  # frozen electrodes (default for static electrodes)
fix mytemp elyte langevin 300 300 100 12345  # Langevin on electrolyte only
```

- Electrodes are frozen (no dynamics) — common for studying equilibrium double layer
- If electrodes should move: use `fix nve` on them and thermostat appropriately
- Langevin adds stochastic noise; good for temperature control in non-equilibrium

### 7.6 Energy Minimization Before Production

```lammps
minimize 1e-6 1e-9 1000 10000    # before running dynamics
```

- Remove overlaps; establish reasonable starting configuration
- Without minimization: first steps may have huge forces → simulation destabilizes

### 7.7 Production Run and Output

```lammps
thermo 100
thermo_style custom step temp pe ke vol density

dump 1 all custom 100 traj.lammpstrj id type x y z q
dump 2 all dcd 100 traj.dcd

run 10000
```

- `thermo` prints every 100 steps; check for energy conservation
- Trajectory file for post-processing: density profiles, RDF, MSD

### 7.8 Expected Output (And What It Means)

| Quantity | Expected Behavior | What It Tells You |
|----------|------------------|-------------------|
| `Temp` | Fluctuates around 300 K | Thermostat working |
| `PotEng` | Stable (small drift) | Force field reasonable |
| `Density` | ~1 g/cm³ for water | Water model correct |
| `f_fxconp[1]` | ~0 (energy conserved) | Electrode fix correct |
| Electrode charge | Fluctuates around ±Q | Double layer formation |

---

## 8. Practical Considerations and Troubleshooting

### 8.1 Convergence Issues

#### 8.1.1 Symptoms
- Charge oscillates wildly between steps
- `nan` or `inf` in thermodynamic output
- Energy increases without bound

#### 8.1.2 Causes and Fixes
| Cause | Fix |
|-------|-----|
| Overlapping atoms | Run minimization; check `pair_coeff` |
| eta too small | Increase eta (wider Gaussian = better-conditioned matrix) |
| eta too large | Decrease eta (narrower = more accurate, but harder to converge) |
| Tolerance too tight | `kspace_style pppm 1e-4` (less accurate but easier) |
| Poor geometry | Check region definitions; ensure no vacuum gaps in electrode |

### 8.2 Electrode Immobilization
- `mat_inv` and `mat_cg` precompute based on electrode positions
- If electrodes move: either use `algo cg` (slow but允许 movement) or recompute matrix
- For charging dynamics: `ffield on` with `algo cg` is the standard approach

### 8.3 Parallel Performance
- Electrode atoms should be evenly distributed: `processors * * 2` (2D slab decomposition)
- Large electrode matrix: can exceed 0.5 GiB per MPI rank
- For >50,000 electrode atoms: use `mat_cg` or `cg`

### 8.4 Units Summary (Critical Reference)

| Quantity | real units | metal units |
|----------|-----------|-------------|
| Distance | Å | Å |
| Energy | kcal/mol | eV |
| Time | fs | fs |
| Potential | kcal/mol·e | V (always!) |
| eta | Å⁻¹ | Å⁻¹ |

- **Potentials are ALWAYS in volts** regardless of `units` setting
- Verify eta is appropriate for your unit system (inverse length)

### 8.5 Reproducibility
- Always set thermostat seeds explicitly: `12345` in `fix langevin`
- Set `random/philox` seed if using GPU (LAMMPS 2024+)
- Matrix save/load ensures same capacitance matrix on restart

### 8.6 Common Pitfalls Summary

| Pitfall | Consequence | Prevention |
|---------|-------------|-------------|
| Wrong units in pair_coeff | Wrong energies; crash | Check units matching `units` setting |
| Mismatched pair + kspace | Error or wrong forces | `coul/long` → `pppm`; `coul/cut` → no kspace |
| Forgot `kspace_modify slab` | Artificially high z-pressures | Add for non-periodic z |
| Electrode atoms in wrong group | Fix applies to wrong atoms | Verify with `group bot` command |
| Cutoff too small | Energy drift; bad electrostatics | Use ≥10 Å for water |

---

## 9. Example Systems

### 9.1 Planar Electrode with Aqueous NaCl
- Two planar gold electrodes, 1M NaCl water
- Study: double-layer capacitance, ion density profiles, potential distribution
- Expected output: capacitance ~10–20 μF/cm² for gold

### 9.2 Graphene Supercapacitor
- Two graphene layers (porous electrode)
- Ionic liquid electrolyte
- Study: charging dynamics, ion insertion, energy density
- CPM essential: charge localizes at sharp curvature

### 9.3 Tip-Substrate Junction
- STM tip near metal surface
- Study: tunneling current, field emission
- Complex geometry: CPM captures tip charge distribution correctly

### 9.4 Finite-Field Charging Simulation
- `boundary p p p` with `ffield on`
- Smaller z-dimension possible (no vacuum needed)
- Model: battery discharge curve (potential vs. charge)
- Use `electrode/conq` to set charge and measure potential response

---

## 10. Further Reading

### 10.1 Foundational Papers
- Siepmann & Sprik, J. Chem. Phys. 102, 511 (1995) — original CPM
- Reed et al., J. Chem. Phys. 126, 084704 (2007) — modern reformulation
- Ahrens-Iwers & Meissner, J. Chem. Phys. 155, 104104 (2021) — ELECTRODE package
- Ahrens-Iwers et al., J. Chem. Phys. 157, 084801 (2022) — package validation

### 10.2 Related Methods
- Gingrich, MSc thesis (2010) — Gaussian smearing for CPM
- Deissenbeck et al., Phys. Rev. Lett. 126, 136803 (2021) — thermopotentiostat
- Dufils et al., Phys. Rev. Lett. 123, 195501 (2019) — finite-field CPM
- Scalfi et al., J. Chem. Phys. 153, 174704 (2020) — Thomas-Fermi model

### 10.3 LAMMPS Documentation
- [ELECTRODE package details](https://docs.lammps.org/Packages_details.html#electrode-package)
- [fix electrode/conp](https://docs.lammps.org/fix_electrode.html)
- [kspace_style](https://docs.lammps.org/kspace_style.html) — electrode variants
- [Examples/PACKAGES/electrode](https://github.com/lammps/lammps/tree/develop/examples/PACKAGES/electrode)

---

## Design Notes

This outline targets a student who has completed introductory physics and chemistry, has access to a working LAMMPS installation, and wants to understand the physics behind the commands.

The progression is:
1. **What MD does** → trajectory, forces, ensemble
2. **What interacts** → force fields (LJ, Coulomb, bonded)
3. **What is controlled** → thermostats, barostats, timestep
4. **How electrostatics is solved** → Ewald, PPPM, slab corrections
5. **Why electrodes are different** → conductor vs. insulator
6. **What CPM does** → energy minimization under potential constraint
7. **How to run it** → ELECTRODE package syntax
8. **How to build a simulation** → full example with output interpretation

Every section connects input syntax to physical meaning and output observables.