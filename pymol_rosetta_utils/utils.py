"""PyMOL utility scripts that use Rosetta as a calculator."""
from collections import namedtuple
import logging
from pymol import cmd
import pyrosetta
import pyrosetta.distributed.io as io
import pyrosetta.distributed.packed_pose as packed_pose

logger = logging.getLogger("pymol_rosetta_utils.utils")
ResidueID = namedtuple("ResidueID", "chain res icode")


def pymol_to_rosetta(pymol_obj, state=-1):
    """Convert a pymol selection to a rosetta selection."""
    return io.pose_from_pdbstring(cmd.get_pdbstr(pymol_obj, state=state))


def rosetta_to_pymol(pack_or_pose, pymol_obj_name="pose"):
    """Convert a rosetta pose to a pymol object."""
    wpose = packed_pose.to_pose(pack_or_pose)
    pdb_string_stream = pyrosetta.rosetta.std.stringstream()
    pyrosetta.rosetta.core.io.pdb.dump_pdb(wpose, pdb_string_stream)
    cmd.read_pdbstr(pdb_string_stream.str(), pymol_obj_name)


def pymol_selection_to_pdb_resnums(selection="all"):
    """Convert a pymol selection to a list of ResidueIDs and return it."""
    local_space = {"raw_resnums": []}
    cmd.iterate(selection, "raw_resnums.append((chain, resv, resi))", space=local_space)

    pdb_resnums = []
    # filter duplicate entries while preserving order (py3.7+ only)
    for chain, resno_int, resno_str in dict.fromkeys(local_space["raw_resnums"]).keys():
        try:
            int(resno_str)
        except ValueError as val_err:
            logger.warning(
                "Detected non-integer value for residue number - %s - treating as inserion code",
                resno_str,
            )
            icode = resno_str[-1]
            if int(resno_str[:-1]) != resno_int:
                raise ValueError from val_err
            pdb_resnums.append(ResidueID(chain, resno_int, icode))
        else:
            pdb_resnums.append(ResidueID(chain, resno_int, " "))

    return pdb_resnums


def pymol_selection_to_residue_selector(
    selection="all", name="index_selector", pose=None
):
    """Create a Rosetta ResidueIndexSelector XML tag from a PyMOL selection and return it."""
    if pose is None:
        pose = packed_pose.to_pose(pymol_to_rosetta(selection))
    pdb_resnums = pymol_selection_to_pdb_resnums(selection)

    pose_resnums = ",".join(
        str(pose.pdb_info().pdb2pose(rsd.chain, rsd.res, rsd.icode))
        for rsd in pdb_resnums
    )
    return f"""<Index name="{name}" resnums="{pose_resnums}" />"""


def really_silly_test(pymol_obj_name):
    """Make a pose from sequence and show it in PyMOL."""
    raw_input_pose = io.pose_from_sequence("TESTSEATTLETEST")
    rosetta_to_pymol(raw_input_pose, pymol_obj_name)
