"""
PyMOL visualization shortcuts that use Rosetta.
"""
__author__ = "Brian D. Weitzner"
__email__ = "brian.weitzner@gmail.com"

import logging
from .utils import rosetta_to_pymol, pymol_to_rosetta
from .selectors import show_interface_residues

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger("pymol_rosetta_utils")
