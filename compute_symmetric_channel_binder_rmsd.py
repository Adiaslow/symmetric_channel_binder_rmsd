#!/usr/bin/env python3

from Bio import PDB
from Bio.PDB.Superimposer import Superimposer
import numpy as np
from pathlib import Path
import argparse
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Union

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
        if filepath.endswith('.pdb'):
            return self.parser.get_structure('structure', filepath)
        elif filepath.endswith('.cif'):
            return self.mmcif_parser.get_structure('structure', filepath)
        raise ValueError(f"Unsupported file format: {filepath}")

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
                    raise ValueError(f"Target sequence not found in chain {chain_id}")

                target_residues = residues[start_idx:start_idx + len(target_sequence)]
            else:
                target_residues = [res for res in residues if res.id[0] == ' ']

            ca_atoms = [residue['CA'] for residue in target_residues if 'CA' in residue]
            return ca_atoms

        except Exception as e:
            raise ValueError(f"Error getting CA atoms for chain {chain_id}: {str(e)}")

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
            # Get CA atoms for channel alignment
            ref_channel_ca = self.get_ca_atoms(ref_structure, ref_chan_chain, target_sequence)
            test_channel_ca = self.get_ca_atoms(test_structure, test_chan_chain, target_sequence)

            if not ref_channel_ca or not test_channel_ca:
                raise ValueError("Failed to get CA atoms for channel alignment")

            # Calculate and apply superposition
            self.sup.set_atoms(ref_channel_ca, test_channel_ca)
            for atom in test_structure.get_atoms():
                atom.set_coord(np.dot(atom.get_coord(), self.sup.rotran[0]) + self.sup.rotran[1])

            # Get binder atoms and calculate RMSD
            binder_ref_atoms = self.get_ca_atoms(ref_structure, ref_bind_chain)
            binder_test_atoms = self.get_ca_atoms(test_structure, test_bind_chain)

            if not binder_ref_atoms or not binder_test_atoms:
                raise ValueError("Failed to get CA atoms for binder comparison")

            # Calculate RMSDs
            min_length = min(len(binder_ref_atoms), len(binder_test_atoms))
            binder_ref_coords = np.array([atom.get_coord() for atom in binder_ref_atoms])[:min_length]
            binder_test_coords = np.array([atom.get_coord() for atom in binder_test_atoms])[:min_length]

            binder_rmsd = np.sqrt(np.mean(np.sum((binder_ref_coords - binder_test_coords) ** 2, axis=1)))
            channel_rmsd = self.sup.rms

            return AnalysisResult(
                binder_rmsd=binder_rmsd,
                channel_rmsd=channel_rmsd,
                binder_length_ref=len(binder_ref_atoms),
                binder_length_model=len(binder_test_atoms),
                channel_length=len(ref_channel_ca)
            )

        except Exception as e:
            raise ValueError(f"Error analyzing structures: {str(e)}")

def parse_args():
    parser = argparse.ArgumentParser(description='Calculate RMSD between reference and test structures')
    parser.add_argument('-refpdb', required=True, help='Path to reference PDB/CIF file')
    parser.add_argument('-refchanchain', required=True, help='Chain ID of channel in reference structure')
    parser.add_argument('-refbindchain', required=True, help='Chain ID of binder in reference structure')
    parser.add_argument('-testpdb', required=True, help='Path to test PDB/CIF file')
    parser.add_argument('-testchanchain', required=True, help='Chain ID of channel in test structure')
    parser.add_argument('-testbindchain', required=True, help='Chain ID of binder in test structure')
    parser.add_argument('-targetseq', help='Target sequence to align (optional)')
    return parser.parse_args()

def main():
    args = parse_args()
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
        print(f"Error: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()
