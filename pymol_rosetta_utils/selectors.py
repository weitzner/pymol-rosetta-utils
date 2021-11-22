"""Functions that display the results of Rosetta ResidueSelectors in PyMOL"""
from collections import namedtuple
import logging
from pymol import cmd
import pyrosetta.distributed.packed_pose as packed_pose
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from .utils import (
    pymol_to_rosetta,
    rosetta_to_pymol,
    pymol_selection_to_residue_selector,
)

logger = logging.getLogger("pymol_rosetta_utils.selectors")
PDBAtomID = namedtuple("PDBAtomID", "resi chain name")


def _apply_residue_selector_and_show_as_sticks(
    selector, pose, label="selected_residues"
):
    """Apply a residue selector to a pose and show the selected residues in sticks."""
    residues = selector.apply(pose)
    sele = []
    for pos, val in enumerate(residues, start=1):
        if not val:
            continue
        # this will fail on insertion codes because rosetta only returns
        # residue numbers and chainIDs
        resno, chain = pose.pdb_info().pose2pdb(pos).split()
        sele.append("(resi {} and chain {})".format(resno, chain))
    cmd.select(label, "+".join(sele))
    cmd.show("sticks", label)
    return residues


def _hbond_pdb(hbond, pose):
    """Convert donor and acceptor atom ids to PDB numbering and return a tuple
    in the form of (donor, acceptor).
    """
    (don_res, don_chain), don_atm = (
        pose.pdb_info().pose2pdb(hbond.don_res()).split(),
        pose.residue(hbond.don_res()).atom_name(hbond.don_hatm()).strip(),
    )

    (acc_res, acc_chain), acc_atm = (
        pose.pdb_info().pose2pdb(hbond.acc_res()).split(),
        pose.residue(hbond.acc_res()).atom_name(hbond.acc_atm()).strip(),
    )

    return PDBAtomID(don_res, don_chain, don_atm), PDBAtomID(
        acc_res, acc_chain, acc_atm
    )


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
    logger.debug(objs)
    _apply_residue_selector_and_show_as_sticks(
        objs.get_residue_selector("interface"), wpose, label="interface_residues"
    )


def show_hydrogen_bonds(pymol_obj_name, include_bb_bb=True, selection="all"):
    """Show hydrogen bonds between two selections."""
    wpose = packed_pose.to_pose(pymol_to_rosetta(pymol_obj_name))
    # update the PyMOL object to include all hydrogens
    rosetta_to_pymol(wpose, pymol_obj_name)

    objs = XmlObjects.create_from_string(
        f"""
<RESIDUE_SELECTORS>
    {pymol_selection_to_residue_selector(selection, name="selector", pose=wpose)}
    <HBond name="hbonds" residue_selector="selector" include_bb_bb="{include_bb_bb}" scorefxn="REF2015"/>
</RESIDUE_SELECTORS>
"""
    )
    logger.debug(objs)

    # show sticks for selected residues
    selected_residues = _apply_residue_selector_and_show_as_sticks(
        objs.get_residue_selector("hbonds"), wpose, label="hbond_residues"
    )

    # hide non-polar hydrogens
    cmd.hide("everything", "(h. and (e. c extend 1))")

    # get all hbonds in pose
    for i, hbond in enumerate(
        wpose.get_hbonds(exclude_bb=(not include_bb_bb)).hbonds(), start=1
    ):
        if (
            not selected_residues[hbond.don_res()]
            or not selected_residues[hbond.acc_res()]
        ):
            continue
        # get the residues and atom names involved in the hbond
        don, acc = _hbond_pdb(hbond, wpose)
        # draw dashed line to indicate hbond
        cmd.distance(
            f"hb_{i}",
            f"(resi {don.resi} and chain {don.chain} and name {don.name})",
            f"(resi {acc.resi} and chain {acc.chain} and name {acc.name})",
        )
        cmd.hide("labels", f"hb_{i}")
    # group all of the hbonds together in a collapsable list
    cmd.group("hbonds", "hb_*")
