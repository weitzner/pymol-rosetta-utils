"""PyMOL utility scripts that use Rosetta as a calculator."""
from pymol import cmd
import pyrosetta
import pyrosetta.distributed.io as io
import pyrosetta.distributed.packed_pose as packed_pose


def pymol_to_rosetta(pymol_obj, state=-1):
    """Convert a pymol selection to a rosetta selection."""
    return io.pose_from_pdbstring(cmd.get_pdbstr(pymol_obj, state=state))


def rosetta_to_pymol(pack_or_pose, pymol_obj_name):
    """Convert a rosetta pose to a pymol object."""
    wpose = packed_pose.to_pose(pack_or_pose)
    pdb_string_stream = pyrosetta.rosetta.std.stringstream()
    pyrosetta.rosetta.core.io.pdb.dump_pdb(wpose, pdb_string_stream)
    cmd.read_pdbstr(pdb_string_stream.str(), pymol_obj_name)


def really_silly_test(pymol_obj_name):
    """Make a pose from sequence and show it in PyMOL."""
    raw_input_pose = io.pose_from_sequence("TESTSEATTLETEST")
    rosetta_to_pymol(raw_input_pose, pymol_obj_name)
