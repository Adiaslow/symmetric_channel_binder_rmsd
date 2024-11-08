from Bio import PDB
from Bio.PDB.Superimposer import Superimposer
import numpy as np
from pathlib import Path
import argparse
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Union
import sys
import logging

# Set up logging
logging.basicConfig(level=logging.INFO,
                   format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class AnalysisResult:
    binder_rmsd: float
    channel_rmsd: float
    binder_length_ref: int
    binder_length_model: int
    channel_length: int

class ChannelBinderAnalyzer:
    def __init__(self, quiet: bool = True):
        self.parser = PDB.PDBParser(QUIET=quiet)
        self.mmcif_parser = PDB.MMCIFParser(QUIET=quiet)
        self.sup = Superimposer()

    def load_structure(self, filepath: Union[str, Path]) -> PDB.Structure.Structure:
        """Load structure from either PDB or mmCIF file."""
        filepath = str(filepath)
        logger.info(f"Loading structure from: {filepath}")
        if filepath.endswith('.pdb'):
            return self.parser.get_structure('structure', filepath)
        elif filepath.endswith('.cif'):
            return self.mmcif_parser.get_structure('structure', filepath)
        raise ValueError(f"Unsupported file format: {filepath}")

    def list_chains_and_sequences(self, structure: PDB.Structure.Structure) -> Dict[str, str]:
        """List all chains and their sequences in the structure."""
        three_to_one = {
            'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F',
            'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
            'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R',
            'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
        }

        sequences = {}
        for chain in structure[0]:
            sequence = ''.join(three_to_one.get(res.get_resname().strip(), 'X')
                             for res in chain.get_residues() if res.id[0] == ' ')
            sequences[chain.id] = sequence
            logger.info(f"Chain {chain.id}: {sequence}")
        return sequences

    def find_chain_with_sequence(self, structure: PDB.Structure.Structure,
                               target_sequence: str) -> Optional[str]:
        """Find the chain containing the target sequence."""
        sequences = self.list_chains_and_sequences(structure)
        for chain_id, sequence in sequences.items():
            if target_sequence in sequence:
                return chain_id
        return None

    def get_ca_atoms(self, structure: PDB.Structure.Structure, chain_id: str,
                     target_sequence: Optional[str] = None) -> List[PDB.Atom.Atom]:
        """Get CA atoms for specified chain, optionally filtering by target sequence."""
        try:
            chain = structure[0][chain_id]
            residues = list(chain.get_residues())

            if target_sequence:
                three_to_one = {
                    'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F',
                    'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
                    'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R',
                    'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
                }

                sequence = ''.join(three_to_one.get(res.get_resname().strip(), 'X')
                                 for res in residues if res.id[0] == ' ')
                start_idx = sequence.find(target_sequence)

                if start_idx == -1:
                    logger.error(f"Target sequence not found in chain {chain_id}")
                    logger.error(f"Chain sequence: {sequence}")
                    logger.error(f"Looking for: {target_sequence}")
                    return []

                target_residues = residues[start_idx:start_idx + len(target_sequence)]
            else:
                target_residues = [res for res in residues if res.id[0] == ' ']

            ca_atoms = [residue['CA'] for residue in target_residues if 'CA' in residue]
            return ca_atoms

        except Exception as e:
            logger.error(f"Error getting CA atoms for chain {chain_id}: {str(e)}")
            return []

    def analyze_structures(self,
                         ref_structure: PDB.Structure.Structure,
                         test_structure: PDB.Structure.Structure,
                         target_sequence: str) -> Optional[AnalysisResult]:
        """Analyze a pair of structures and return their RMSD values."""
        try:
            # Find chains containing target sequence
            ref_channel_chain = self.find_chain_with_sequence(ref_structure, target_sequence)
            test_channel_chain = self.find_chain_with_sequence(test_structure, target_sequence)

            if not ref_channel_chain or not test_channel_chain:
                logger.error("Could not find target sequence in structures")
                return None

            # Get binder chains (assuming it's the other chain)
            ref_chains = set(chain.id for chain in ref_structure[0])
            test_chains = set(chain.id for chain in test_structure[0])

            ref_binder_chain = (ref_chains - {ref_channel_chain}).pop()
            test_binder_chain = (test_chains - {test_channel_chain}).pop()

            # Get CA atoms for channel alignment
            fixed_ca = self.get_ca_atoms(ref_structure, ref_channel_chain, target_sequence)
            moving_ca = self.get_ca_atoms(test_structure, test_channel_chain, target_sequence)

            if not fixed_ca or not moving_ca:
                logger.error("Failed to get CA atoms for channel alignment")
                return None

            # Calculate and apply superposition
            self.sup.set_atoms(fixed_ca, moving_ca)
            for atom in test_structure.get_atoms():
                atom.set_coord(np.dot(atom.get_coord(), self.sup.rotran[0]) + self.sup.rotran[1])

            # Get binder atoms and calculate RMSD
            binder1_atoms = self.get_ca_atoms(ref_structure, ref_binder_chain)
            binder2_atoms = self.get_ca_atoms(test_structure, test_binder_chain)

            if not binder1_atoms or not binder2_atoms:
                logger.error("Failed to get CA atoms for binder comparison")
                return None

            # Calculate RMSDs
            min_length = min(len(binder1_atoms), len(binder2_atoms))
            binder1_coords = np.array([atom.get_coord() for atom in binder1_atoms])[:min_length]
            binder2_coords = np.array([atom.get_coord() for atom in binder2_atoms])[:min_length]

            binder_rmsd = np.sqrt(np.mean(np.sum((binder1_coords - binder2_coords) ** 2, axis=1)))
            channel_rmsd = self.sup.rms

            logger.info(f"Calculated RMSDs - Binder: {binder_rmsd:.2f}, Channel: {channel_rmsd:.2f}")

            return AnalysisResult(
                binder_rmsd=binder_rmsd,
                channel_rmsd=channel_rmsd,
                binder_length_ref=len(binder1_atoms),
                binder_length_model=len(binder2_atoms),
                channel_length=len(fixed_ca)
            )

        except Exception as e:
            logger.error(f"Error analyzing structures: {str(e)}")
            return None

def main():
    parser = argparse.ArgumentParser(description='Calculate RMSD between reference and test structures')
    parser.add_argument('-refpdb', required=True, help='Path to reference PDB/CIF file')
    parser.add_argument('-testpdb', required=True, help='Path to test PDB/CIF file')
    parser.add_argument('-targetseq', required=True, help='Target sequence for channel identification')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')
    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    analyzer = ChannelBinderAnalyzer()
    try:
        # Load structures
        ref_structure = analyzer.load_structure(args.refpdb)
        test_structure = analyzer.load_structure(args.testpdb)

        # Analyze structures
        result = analyzer.analyze_structures(
            ref_structure=ref_structure,
            test_structure=test_structure,
            target_sequence=args.targetseq
        )

        if result:
            print("\nAnalysis Results:")
            print(f"Binder RMSD: {result.binder_rmsd:.2f}")
            print(f"Channel RMSD: {result.channel_rmsd:.2f}")
            print(f"Binder length (reference): {result.binder_length_ref}")
            print(f"Binder length (test): {result.binder_length_model}")
            print(f"Channel length: {result.channel_length}")
        else:
            sys.exit(1)

    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

if __name__ == "__main__":
    main()
