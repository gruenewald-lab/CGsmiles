"""
CGsmiles: Coarse-Grained Smiles (CGsmiles) for representing abitrarily complex molecules using line notation.
"""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version('cgsmiles')
except PackageNotFoundError:
    __version__ = '0+unknown'

from .read_cgsmiles import read_cgsmiles
from .read_fragments import read_fragments
from .resolve import MoleculeResolver
from .sample import MoleculeSampler
