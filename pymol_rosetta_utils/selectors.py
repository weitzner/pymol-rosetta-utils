"""Functions that display the results of Rosetta ResidueSelectors in PyMOL"""
from pymol import cmd
import pyrosetta.distributed.packed_pose as packed_pose
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

from .utils import pymol_to_rosetta


def _show_residues_from_residue_selector(selector, pose, label="selected_residues"):
    """Apply a residue selector to a pose and show the selected residues in sticks."""
    residues = selector.apply(pose)
    sele = []
    for pos, val in enumerate(residues, start=1):
        if not val:
            continue
        # note: this will fail on insertion codes
        resno, chain = pose.pdb_info().pose2pdb(pos).split()
        sele.append("(resi {} and chain {})".format(resno, chain))
    cmd.select(label, "+".join(sele))
    cmd.show("sticks", label)


def show_interface_residues(pymol_obj_name, selection_1, selection_2):
    """Show interface residues between two selections."""
    ppose = pymol_to_rosetta(pymol_obj_name)
    wpose = packed_pose.to_pose(ppose)
    objs = XmlObjects.create_from_string(
        f"""
<RESIDUE_SELECTORS>
    <Chain name="chain1" chains="{selection_1}" />
    <Chain name="chain2" chains="{selection_2}" />
    <InterfaceByVector name="interface" grp1_selector="chain1" grp2_selector="chain2"/>
</RESIDUE_SELECTORS>
"""
    )

    interface_selector = objs.get_residue_selector("interface")
    _show_residues_from_residue_selector(
        interface_selector, wpose, label="interface_residues"
    )
