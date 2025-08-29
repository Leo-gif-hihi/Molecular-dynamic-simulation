#!/usr/bin/env python3
"""
openmm_proteinwater.py

A flexible script to run protein-in-water molecular dynamics simulations with OpenMM.

Automatically runs successive production chunks of length --checkpoint-ps
until --total-ns is reached. Safe to interrupt and resume.
"""

import os, sys, argparse
from pdbfixer import PDBFixer
from openmm.app import (Modeller, ForceField, Simulation,
                        DCDReporter, StateDataReporter, CheckpointReporter,
                        PME, HBonds, PDBFile)
from openmm import XmlSerializer, unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator

def parse_args():
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pdbid", help="Input PDB ID or filename.")
    p.add_argument("--total-ns",      type=float, default=1.0,   help="Total production simulation time in nanoseconds. Default: 1.0")
    p.add_argument("--traj-interval", type=float, default=1.0,   help="Trajectory snapshot interval in picoseconds. Default: 1.0")
    p.add_argument("--equil-time",    type=float, default=10.0,  help="Equilibration time in picoseconds. Default: 10.0")
    p.add_argument("--checkpoint-ps", type=float, default=10.0,  help="Checkpoint interval in picoseconds. Default: 10.0")

    # Arguments for physical parameters
    p.add_argument("--temp",          type=float, default=300.0, help="Temperature in Kelvin. Default: 300.0")
    p.add_argument("--pressure",      type=float, default=1.0,   help="Pressure in atmospheres. Default: 1.0")
    p.add_argument("--friction",      type=float, default=1.0,   help="Friction coefficient for Langevin integrator in 1/ps. Default: 1.0")

    # Arguments for force field
    p.add_argument("--ff",            type=str, default='amber14-all.xml', help="Force field XML for the protein. Default: 'amber14-all.xml'")
    p.add_argument("--water",         type=str, default='amber14/tip3p.xml', help="Force field XML for the water model. Default: 'amber14/tip3p.xml'")

    # Argument for platform selection
    p.add_argument("--platform",      type=str, default='auto', help="Computation platform: 'auto', 'CUDA', 'OpenCL', 'CPU'. Default: 'auto'")
    return p.parse_args()

def get_platform(name='auto'):
    """Selects the fastest available platform or a specific one."""
    if name.lower() == 'auto':
        platforms = [Platform.getPlatform(i) for i in range(Platform.getNumPlatforms())]
        # Speed hierarchy: CUDA > OpenCL > CPU
        ranked_platforms = sorted(platforms, key=lambda p: {'CUDA': 0, 'OpenCL': 1, 'CPU': 2}.get(p.getName(), 3))
        selected_platform = ranked_platforms[0]
        print(f"  • Auto-selected platform: {selected_platform.getName()}")
        return selected_platform
    else:
        try:
            platform = Platform.getPlatformByName(name)
            print(f"  • Using specified platform: {name}")
            return platform
        except Exception:
            sys.exit(f"Error: Platform '{name}' not found or not available.")

def make_sim(top, sys_, dt, temp, friction, platform_name):
    """Creates an OpenMM Simulation object."""
    integ = LangevinMiddleIntegrator(temp*unit.kelvin, friction/unit.picosecond, dt)
    plat = get_platform(platform_name)
    return Simulation(top, sys_, integ, plat)

def setup_and_equilibrate(cleaned_pdb, dt, eq_steps, xml_sys, pdb_saved, chk_file,
                          ff_protein, ff_water, temp, pressure, friction, platform_name):
    """Sets up a new simulation system from scratch and equilibrates it."""
    fixer = PDBFixer(filename=cleaned_pdb)
    fixer.removeHeterogens(keepWater=True)
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    fixer.addMissingAtoms(); fixer.addMissingHydrogens(pH=7.0)

    modeller = Modeller(fixer.topology, fixer.positions)
    try:
        ff = ForceField(ff_protein, ff_water)
    except Exception as e:
        sys.exit(f"Error creating ForceField: {e}\nCheck if force field files '{ff_protein}' and '{ff_water}' are available.")

    modeller.addSolvent(ff, model='tip3p', padding=1.0*unit.nanometer, ionicStrength=0.15*unit.molar)

    system = ff.createSystem(modeller.topology,
                             nonbondedMethod=PME,
                             nonbondedCutoff=1.0*unit.nanometer,
                             constraints=HBonds)
    system.addForce(MonteCarloBarostat(pressure*unit.atmosphere, temp*unit.kelvin, 25))

    sim = make_sim(modeller.topology, system, dt, temp, friction, platform_name)
    sim.context.setPositions(modeller.positions)

    print("  • Minimizing…"); sim.minimizeEnergy()
    print(f"  • Equilibrating NVT at {temp} K…")
    sim.context.setVelocitiesToTemperature(temp*unit.kelvin)
    sim.reporters.append(StateDataReporter('nvt.log', eq_steps//10, step=True, temperature=True))
    sim.step(eq_steps); sim.reporters.clear()

    print(f"  • Equilibrating NPT at {temp} K and {pressure} atm…")
    sim.reporters.append(StateDataReporter('npt.log', eq_steps//10, step=True, temperature=True, volume=True))
    sim.step(eq_steps); sim.reporters.clear()

    with open(xml_sys,'w') as f:
        f.write(XmlSerializer.serializeSystem(sim.context.getSystem()))
    PDBFile.writeFile(sim.topology,
                      sim.context.getState(getPositions=True).getPositions(),
                      open(pdb_saved,'w'))
    sim.saveCheckpoint(chk_file)
    return sim, 0

def resume_sim(xml_sys, pdb_saved, chk_file, dt, temp, friction, platform_name):
    """Resumes a simulation from a checkpoint."""
    with open(xml_sys) as f:
        system = XmlSerializer.deserializeSystem(f.read())
    pdb = PDBFile(pdb_saved)
    sim = make_sim(pdb.topology, system, dt, temp, friction, platform_name)
    with open(chk_file,'rb') as f:
        sim.loadCheckpoint(f)
    steps = sim.context.getState().getStepCount()
    print(f"  • Resumed at {steps} steps")
    return sim, steps

def main():
    """Main execution function."""
    args = parse_args()
    pdbid = args.pdbid.lower()
    workdir = os.path.join(os.getcwd(), pdbid)
    cleaned_pdb = os.path.join(workdir, f"{pdbid}_cleaned.pdb")
    if not os.path.exists(cleaned_pdb):
        sys.exit(f"Error: Cleaned PDB not found at {cleaned_pdb}")
    os.chdir(workdir)

    # MD setup
    dt = 2.0*unit.femtoseconds
    dt_ps = dt.value_in_unit(unit.picoseconds)
    total_steps = int((args.total_ns*unit.nanoseconds)/dt)
    snapshot_steps = int((args.traj_interval*unit.picoseconds)/dt)
    eq_steps = int((args.equil_time*unit.picoseconds)/dt)
    chk_steps = int((args.checkpoint_ps*unit.picoseconds)/dt)

    xml_system = 'system.xml'
    pdb_saved = 'solvated.pdb'
    chk_file = 'prod.chk'

    # Setup or resume
    if os.path.isfile(chk_file) and os.path.isfile(xml_system) and os.path.isfile(pdb_saved):
        print("▶ Resuming simulation…")
        sim, steps_done = resume_sim(xml_system, pdb_saved, chk_file, dt,
                                     args.temp, args.friction, args.platform)
    else:
        print("▶ Setting up new simulation…")
        sim, steps_done = setup_and_equilibrate(cleaned_pdb, dt, eq_steps,
                                                xml_system, pdb_saved, chk_file,
                                                args.ff, args.water, args.temp,
                                                args.pressure, args.friction,
                                                args.platform)

    # Production loop
    while steps_done < total_steps:
        steps_to_run = min(chk_steps, total_steps - steps_done)
        ps_start = int(steps_done * dt_ps)
        ps_end = int(ps_start + steps_to_run * dt_ps)
        tag = f"{ps_start}to{ps_end}ps"
        print(f"▶ Running chunk {tag}…")

        # Attach reporters
        sim.reporters = []
        sim.reporters.append(DCDReporter(f'prod_{tag}.dcd', snapshot_steps))
        sim.reporters.append(StateDataReporter(f'prod_{tag}.log', snapshot_steps,
                                               step=True, time=True,
                                               potentialEnergy=True, temperature=True,
                                               volume=True, speed=True))
        sim.reporters.append(CheckpointReporter(chk_file, chk_steps))

        try:
            sim.step(steps_to_run)
        except KeyboardInterrupt:
            print("⚠ Interrupted — checkpoint saved.")
            # No need for finally block as sim.step() handles checkpointing on exceptions
            # when a CheckpointReporter is active. Saving again is redundant but safe.
            sim.saveCheckpoint(chk_file)
            sys.exit(0) # Exit cleanly

        # Update steps_done. Note: OpenMM's step count is already updated internally.
        steps_done = sim.context.getState().getStepCount()
        print(f"  ✔ Completed chunk. Total steps: {steps_done}")

    print("✔︎ Full production simulation complete.")

if __name__ == "__main__":
    main()