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
class ChainMapping:
    channel_chain: str
    binder_chain: str

@dataclass
class AnalysisResult:
    binder_rmsd: float
    channel_rmsd: float
    binder_length_ref: int
    binder_length_model: int
    channel_length: int
    ref_mapping: ChainMapping
    test_mapping: ChainMapping

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

            chain_count = len(list(structure[0]))
            logger.info(f"Successfully loaded structure with {chain_count} chains")
            return structure
        except Exception as e:
            logger.error(f"Failed to load structure {filepath}: {str(e)}")
            raise

    def get_chain_sequences(self, structure: PDB.Structure.Structure) -> Dict[str, str]:
        """Get sequences for all chains in the structure."""
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
            logger.info(f"Chain {chain.id} sequence: {sequence}")
            logger.info(f"Chain {chain.id} sequence length: {len(sequence)}")
        return sequences

    def find_best_sequence_match(self, sequence: str, target: str) -> Optional[int]:
        """Find the best matching position of target sequence in main sequence."""
        def calculate_similarity(seq1: str, seq2: str) -> float:
            """Calculate sequence similarity score."""
            if len(seq1) != len(seq2):
                return 0.0
            matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
            return matches / len(seq1)

        target = target.upper()
        sequence = sequence.upper()

        # Try exact match first
        exact_match = sequence.find(target)
        if exact_match != -1:
            logger.info(f"Found exact sequence match at position {exact_match}")
            return exact_match

        # Try with flexible matching
        best_score = 0.0
        best_position = -1
        target_length = len(target)

        for i in range(len(sequence) - target_length + 1):
            window = sequence[i:i + target_length]
            score = calculate_similarity(window, target)
            logger.debug(f"Position {i}: {window} -> similarity score: {score:.2f}")
            if score > best_score and score > 0.8:  # 80% similarity threshold
                best_score = score
                best_position = i

        if best_position != -1:
            logger.info(f"Found best sequence match at position {best_position} with {best_score:.2f} similarity")
            logger.info(f"Matched sequence: {sequence[best_position:best_position + target_length]}")
            return best_position

        return None

    def identify_chains(self, structure: PDB.Structure.Structure,
                       target_sequence: str) -> ChainMapping:
        """Identify channel and binder chains based on sequence content."""
        sequences = self.get_chain_sequences(structure)

        # Find channel chain (contains target sequence)
        channel_chain = None
        best_match_position = None

        for chain_id, sequence in sequences.items():
            match_position = self.find_best_sequence_match(sequence, target_sequence)
            if match_position is not None:
                channel_chain = chain_id
                best_match_position = match_position
                logger.info(f"Found channel sequence in chain {chain_id} at position {match_position}")
                break

        if not channel_chain:
            raise ValueError("Could not find channel sequence in any chain")

        # The other chain must be the binder
        binder_chain = next(chain_id for chain_id in sequences.keys() if chain_id != channel_chain)
        logger.info(f"Using chain {binder_chain} as binder chain")

        return ChainMapping(channel_chain=channel_chain, binder_chain=binder_chain)

    def get_ca_atoms(self, structure: PDB.Structure.Structure, chain_id: str,
                     target_sequence: Optional[str] = None) -> List[PDB.Atom.Atom]:
        """Get CA atoms for specified chain, optionally filtering by target sequence."""
        logger.info(f"Getting CA atoms for chain {chain_id}")

        try:
            chain = structure[0][chain_id]
            residues = list(filter(lambda r: r.id[0] == ' ', chain.get_residues()))
            logger.info(f"Found {len(residues)} standard residues in chain")

            if target_sequence:
                sequence = self.get_chain_sequences(structure)[chain_id]
                match_position = self.find_best_sequence_match(sequence, target_sequence)

                if match_position is None:
                    raise ValueError(f"Target sequence not found in chain {chain_id}")

                logger.info(f"Using sequence match at position {match_position}")
                target_residues = residues[match_position:match_position + len(target_sequence)]
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
                         target_sequence: str) -> AnalysisResult:
        """Analyze a pair of structures and return their RMSD values."""
        try:
            # Identify chains in both structures
            ref_mapping = self.identify_chains(ref_structure, target_sequence)
            test_mapping = self.identify_chains(test_structure, target_sequence)

            # Get CA atoms for channel alignment
            ref_channel_ca = self.get_ca_atoms(ref_structure, ref_mapping.channel_chain, target_sequence)
            test_channel_ca = self.get_ca_atoms(test_structure, test_mapping.channel_chain, target_sequence)

            # Calculate and apply superposition
            logger.info("Performing superposition")
            self.sup.set_atoms(ref_channel_ca, test_channel_ca)
            for atom in test_structure.get_atoms():
                atom.set_coord(np.dot(atom.get_coord(), self.sup.rotran[0]) + self.sup.rotran[1])

            # Get binder atoms and calculate RMSD
            binder_ref_atoms = self.get_ca_atoms(ref_structure, ref_mapping.binder_chain)
            binder_test_atoms = self.get_ca_atoms(test_structure, test_mapping.binder_chain)

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
                channel_length=len(ref_channel_ca),
                ref_mapping=ref_mapping,
                test_mapping=test_mapping
            )

        except Exception as e:
            logger.error(f"Error analyzing structures: {str(e)}")
            raise

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

        # Print results
        print("\nAnalysis Results:")
        print(f"Reference structure mapping:")
        print(f"  Channel chain: {result.ref_mapping.channel_chain}")
        print(f"  Binder chain: {result.ref_mapping.binder_chain}")
        print(f"Test structure mapping:")
        print(f"  Channel chain: {result.test_mapping.channel_chain}")
        print(f"  Binder chain: {result.test_mapping.binder_chain}")
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
