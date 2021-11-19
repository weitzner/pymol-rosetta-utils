"""PyMOL utility scripts that use Rosetta as a calculator."""
from pymol import cmd
import pyrosetta
import pyrosetta.distributed.io as io
import pyrosetta.distributed.packed_pose as packed_pose


def pymol_to_rosetta(pymol_obj, state=-1):
    """Convert a pymol selection to a rosetta selection."""
    pdb_string = cmd.get_pdbstr(pymol_obj, state=state)
    pose = io.pose_from_pdbstring(pdb_string)
    return pose


def rosetta_to_pymol(pack_or_pose, pymol_obj_name):
    """Convert a rosetta pose to a pymol object."""
    wpose = packed_pose.to_pose(pack_or_pose)
    string_stream = pyrosetta.rosetta.std.stringstream()
    pyrosetta.rosetta.core.io.pdb.dump_pdb(wpose, string_stream)
    pdb_string = string_stream.str()
    cmd.read_pdbstr(pdb_string, pymol_obj_name)


def really_silly_test(pymol_obj_name):
    """Make a pose from sequence and show it in PyMOL."""
    import pyrosetta.distributed.tasks.score as score

    raw_input_pose = score.ScorePoseTask()(io.pose_from_sequence("TESTSEATTLETEST"))
    rosetta_to_pymol(raw_input_pose, pymol_obj_name)
