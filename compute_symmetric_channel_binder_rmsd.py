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
        try:
            if filepath.endswith('.pdb'):
                structure = self.parser.get_structure('structure', filepath)
            elif filepath.endswith('.cif'):
                structure = self.mmcif_parser.get_structure('structure', filepath)
            else:
                raise ValueError(f"Unsupported file format: {filepath}")

            # Verify structure loaded successfully
            chain_count = len(list(structure[0]))
            logger.info(f"Successfully loaded structure with {chain_count} chains")
            return structure
        except Exception as e:
            logger.error(f"Failed to load structure {filepath}: {str(e)}")
            raise

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
                             for res in chain.get_residues()
                             if res.id[0] == ' ' and res.get_resname().strip() in three_to_one)
            sequences[chain.id] = sequence
            logger.info(f"Chain {chain.id} sequence length: {len(sequence)}")
        return sequences

    def verify_chain_ids(self, structure: PDB.Structure.Structure, expected_chains: List[str]) -> bool:
        """Verify that all expected chains exist in the structure."""
        available_chains = set(chain.id for chain in structure[0])
        logger.info(f"Available chains: {available_chains}")
        logger.info(f"Expected chains: {expected_chains}")

        missing_chains = set(expected_chains) - available_chains
        if missing_chains:
            logger.error(f"Missing chains: {missing_chains}")
            return False
        return True

    def find_matching_chain(self, structure: PDB.Structure.Structure,
                          target_sequence: str,
                          required_chain: Optional[str] = None) -> Optional[str]:
        """Find chain containing target sequence, optionally verifying against required chain."""
        sequences = self.list_chains_and_sequences(structure)

        if required_chain:
            if required_chain not in sequences:
                logger.error(f"Required chain {required_chain} not found in structure")
                return None

            if target_sequence in sequences[required_chain]:
                return required_chain
            else:
                logger.error(f"Target sequence not found in required chain {required_chain}")
                logger.error(f"Chain {required_chain} sequence: {sequences[required_chain]}")
                logger.error(f"Target sequence: {target_sequence}")
                return None

        # If no required chain, search all chains
        for chain_id, sequence in sequences.items():
            if target_sequence in sequence:
                logger.info(f"Found target sequence in chain {chain_id}")
                return chain_id

        logger.error("Target sequence not found in any chain")
        return None

    def get_ca_atoms(self, structure: PDB.Structure.Structure, chain_id: str,
                     target_sequence: Optional[str] = None) -> List[PDB.Atom.Atom]:
        """Get CA atoms for specified chain, optionally filtering by target sequence."""
        logger.info(f"Getting CA atoms for chain {chain_id}")

        try:
            if chain_id not in structure[0]:
                raise ValueError(f"Chain {chain_id} not found in structure")

            chain = structure[0][chain_id]
            residues = list(filter(lambda r: r.id[0] == ' ', chain.get_residues()))
            logger.info(f"Found {len(residues)} standard residues in chain")

            if target_sequence:
                three_to_one = {
                    'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F',
                    'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
                    'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R',
                    'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
                }

                sequence = ''.join(three_to_one.get(res.get_resname().strip(), 'X')
                                 for res in residues)

                logger.info(f"Chain sequence: {sequence}")
                logger.info(f"Target sequence: {target_sequence}")

                start_idx = sequence.find(target_sequence)
                if start_idx == -1:
                    raise ValueError(f"Target sequence not found in chain {chain_id}")

                logger.info(f"Found target sequence starting at position {start_idx}")
                target_residues = residues[start_idx:start_idx + len(target_sequence)]
            else:
                target_residues = residues

            ca_atoms = []
            for residue in target_residues:
                if 'CA' not in residue:
                    logger.warning(f"No CA atom in residue {residue.get_resname()} {residue.id[1]}")
                else:
                    ca_atoms.append(residue['CA'])

            if not ca_atoms:
                raise ValueError(f"No CA atoms found in chain {chain_id}")

            logger.info(f"Found {len(ca_atoms)} CA atoms")
            return ca_atoms

        except Exception as e:
            logger.error(f"Error getting CA atoms: {str(e)}")
            raise

    def analyze_structures(self,
                         ref_structure: PDB.Structure.Structure,
                         test_structure: PDB.Structure.Structure,
                         ref_chan_chain: str,
                         ref_bind_chain: str,
                         test_chan_chain: str,
                         test_bind_chain: str,
                         target_sequence: Optional[str] = None) -> AnalysisResult:
        """Analyze a pair of structures and return their RMSD values."""
        try:
            # Verify all chains exist
            ref_chains = [ref_chan_chain, ref_bind_chain]
            test_chains = [test_chan_chain, test_bind_chain]

            if not self.verify_chain_ids(ref_structure, ref_chains):
                raise ValueError("Missing required chains in reference structure")
            if not self.verify_chain_ids(test_structure, test_chains):
                raise ValueError("Missing required chains in test structure")

            # Verify target sequence is present in channel chains if provided
            if target_sequence:
                if not self.find_matching_chain(ref_structure, target_sequence, ref_chan_chain):
                    raise ValueError("Target sequence not found in reference channel chain")
                if not self.find_matching_chain(test_structure, target_sequence, test_chan_chain):
                    raise ValueError("Target sequence not found in test channel chain")

            # Get CA atoms for channel alignment
            ref_channel_ca = self.get_ca_atoms(ref_structure, ref_chan_chain, target_sequence)
            test_channel_ca = self.get_ca_atoms(test_structure, test_chan_chain, target_sequence)

            # Calculate and apply superposition
            logger.info("Performing superposition")
            self.sup.set_atoms(ref_channel_ca, test_channel_ca)
            for atom in test_structure.get_atoms():
                atom.set_coord(np.dot(atom.get_coord(), self.sup.rotran[0]) + self.sup.rotran[1])

            # Get binder atoms and calculate RMSD
            binder_ref_atoms = self.get_ca_atoms(ref_structure, ref_bind_chain)
            binder_test_atoms = self.get_ca_atoms(test_structure, test_bind_chain)

            # Calculate RMSDs
            min_length = min(len(binder_ref_atoms), len(binder_test_atoms))
            binder_ref_coords = np.array([atom.get_coord() for atom in binder_ref_atoms])[:min_length]
            binder_test_coords = np.array([atom.get_coord() for atom in binder_test_atoms])[:min_length]

            binder_rmsd = np.sqrt(np.mean(np.sum((binder_ref_coords - binder_test_coords) ** 2, axis=1)))
            channel_rmsd = self.sup.rms

            logger.info(f"Analysis complete - Binder RMSD: {binder_rmsd:.2f}, Channel RMSD: {channel_rmsd:.2f}")
            return AnalysisResult(
                binder_rmsd=binder_rmsd,
                channel_rmsd=channel_rmsd,
                binder_length_ref=len(binder_ref_atoms),
                binder_length_model=len(binder_test_atoms),
                channel_length=len(ref_channel_ca)
            )

        except Exception as e:
            logger.error(f"Error analyzing structures: {str(e)}")
            raise

def main():
    parser = argparse.ArgumentParser(description='Calculate RMSD between reference and test structures')
    parser.add_argument('-refpdb', required=True, help='Path to reference PDB/CIF file')
    parser.add_argument('-refchanchain', required=True, help='Chain ID of channel in reference structure')
    parser.add_argument('-refbindchain', required=True, help='Chain ID of binder in reference structure')
    parser.add_argument('-testpdb', required=True, help='Path to test PDB/CIF file')
    parser.add_argument('-testchanchain', required=True, help='Chain ID of channel in test structure')
    parser.add_argument('-testbindchain', required=True, help='Chain ID of binder in test structure')
    parser.add_argument('-targetseq', help='Target sequence to align (optional)')
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
            ref_chan_chain=args.refchanchain,
            ref_bind_chain=args.refbindchain,
            test_chan_chain=args.testchanchain,
            test_bind_chain=args.testbindchain,
            target_sequence=args.targetseq
        )

        # Print results
        print("\nAnalysis Results:")
        print(f"Binder RMSD: {result.binder_rmsd:.2f}")
        print(f"Channel RMSD: {result.channel_rmsd:.2f}")
        print(f"Binder length (reference): {result.binder_length_ref}")
        print(f"Binder length (test): {result.binder_length_model}")
        print(f"Channel length: {result.channel_length}")

    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

if __name__ == "__main__":
    main()
