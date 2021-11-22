"""Functions that display the results of Rosetta ResidueSelectors in PyMOL"""
import logging
from pymol import cmd
import pyrosetta.distributed.packed_pose as packed_pose
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from .utils import pymol_to_rosetta, pymol_selection_to_residue_selector

logger = logging.getLogger("pymol_rosetta_utils.selectors")


def _show_residues_from_residue_selector(selector, pose, label="selected_residues"):
    """Apply a residue selector to a pose and show the selected residues in sticks."""
    residues = selector.apply(pose)
    sele = []
    for pos, val in enumerate(residues, start=1):
        if not val:
            continue
        # Note: this will fail on insertion codes (rosetta issue)
        resno, chain = pose.pdb_info().pose2pdb(pos).split()
        sele.append("(resi {} and chain {})".format(resno, chain))
    cmd.select(label, "+".join(sele))
    cmd.show("sticks", label)


def show_interface_residues(pymol_obj_name, grp1, grp2):
    """Show interface residues between two selections."""
    wpose = packed_pose.to_pose(pymol_to_rosetta(pymol_obj_name))
    objs = XmlObjects.create_from_string(
        f"""
<RESIDUE_SELECTORS>
    {pymol_selection_to_residue_selector(grp1, name="grp1", pose=wpose)}
    {pymol_selection_to_residue_selector(grp2, name="grp2", pose=wpose)}
    <InterfaceByVector name="interface" grp1_selector="grp1" grp2_selector="grp2"/>
</RESIDUE_SELECTORS>
"""
    )
    logging.debug(objs)
    _show_residues_from_residue_selector(
        objs.get_residue_selector("interface"), wpose, label="interface_residues"
    )
