import MDAnalysis as mda
import mdtraj as md
import os
import sys
import argparse
import time
import logging
import re

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

start_time = time.time()

# Argument parsing
parser = argparse.ArgumentParser(description="Remove hydrogens from reference and trajectory files.")
parser.add_argument('-r', '--ref', dest='RefFile', required=True, help="Reference file (.gro, .pdb, .psf)")
parser.add_argument('-t', '--traj', dest='TrajFile', required=True, help="Trajectory file (.xtc, .trr, .dcd)")
parser.add_argument('--gro-out', dest='gro_output', required=True, help="Output .gro filename (without extension)")
parser.add_argument('--traj-out', dest='traj_output', required=True, help="Output .xtc filename (without extension)")
args = parser.parse_args()

RefFile, TrajFile = args.RefFile, args.TrajFile
gro_output_filename = args.gro_output.strip()
traj_output_filename = args.traj_output.strip()

# Hardcoded script directory
script_dir = "/home/alessia.guadagnin/propre/lib" 
sys.path.append(script_dir)
logger.info("Using hardcoded script directory: %s", script_dir)

# Import custom modules
from inp_out import *
from check_errors import *

# Validate output filenames
def validate_filename(name):
    if not re.match(r'^[\w\-.]+$', name):
        logger.error("ERROR: Invalid filename '%s'. Use only alphanumeric characters, dashes, and underscores.", name)
        sys.exit(1)

validate_filename(gro_output_filename)
validate_filename(traj_output_filename)

# Validate input files
mandatory_files_present_removeH(RefFile, TrajFile)
checking_file_found(RefFile)
check_empty_file(RefFile)
checking_file_found(TrajFile)
check_empty_file(TrajFile)

# Load molecular structure & trajectory
logger.info("Loading reference and trajectory files...")
u = mda.Universe(RefFile, TrajFile)

# Check if the trajectory already lacks hydrogens
if not u.select_atoms("name H*").n_atoms:
    logger.info("Trajectory already lacks hydrogen atoms. Skipping processing.")
    processed_traj_file = TrajFile
    processed_gro_file = RefFile
else:
    Ref_noH = u.select_atoms("not name H*")

    # Write new reference & trajectory
    Ref_noH.write(f"{gro_output_filename}.gro")
    logger.info("Hydrogen-free reference file saved as %s.gro", gro_output_filename)

    with mda.Writer(f"{traj_output_filename}.xtc", Ref_noH.n_atoms) as W:
        for ts in u.trajectory:
            W.write(Ref_noH)
    logger.info("Hydrogen-free trajectory file saved as %s.xtc", traj_output_filename)
    
    processed_traj_file = f"{traj_output_filename}.xtc"
    processed_gro_file = f"{gro_output_filename}.gro"

# Convert XTC to XYZ efficiently
logger.info("Converting trajectory to XYZ format...")
with md.formats.XYZTrajectoryFile(f"{traj_output_filename}.xyz", 'w') as f:
    for chunk in md.iterload(processed_traj_file, top=processed_gro_file, chunk=100):
        f.write(chunk.xyz * 10)
logger.info("XYZ trajectory saved as %s.xyz", traj_output_filename)

end_time = time.time()
logger.info("✅ Process completed in %.2f seconds.", end_time - start_time)
