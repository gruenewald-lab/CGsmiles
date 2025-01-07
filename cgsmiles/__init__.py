"""
CGsmiles: Coarse-Grained Smiles (CGsmiles) for representing abitrarily complex molecules using line notation.
"""
import pbr.version
__version__ = pbr.version.VersionInfo('cgsmiles').release_string()

from .read_cgsmiles import read_cgsmiles
from .read_fragments import read_fragments
from .resolve import MoleculeResolver
from .sample import MoleculeSampler
