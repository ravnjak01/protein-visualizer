from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')  # <-- OBAVEZNO pre plt importa
import matplotlib.pyplot as plt
import io
import base64
def get_amino_acid_letter(three_letter_code):
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return aa_dict.get(three_letter_code, 'X')


def analyze_pdb(file_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", file_path)

    total_aa_counts = defaultdict(int)
    chain_data = {}

    for model in structure:
        for chain in model:
            aa_counts = defaultdict(int)
            residues = [res for res in chain if is_aa(res)]
            sequence = ""
            atom_count = 0

            for res in residues:
                try:
                    resname = res.get_resname()
                    aa_counts[resname] += 1
                    total_aa_counts[resname] += 1

                    sequence += get_amino_acid_letter(resname)
                    atom_count += len(res)
                except Exception as e:
                    print(f"Preskačem ostatak u lancu {chain.id}: {e}")
                    continue

            chain_data[chain.id] = {
                "residue_count": len(residues),
                "aa_counts": dict(aa_counts),
                "sequence": sequence,
                "atom_count": atom_count
            }

    return {
        "aa_counts": dict(total_aa_counts),
        "chains": chain_data
    }

def create_histogram(aa_counts):
    fig, ax = plt.subplots()
    aa_list = list(aa_counts.keys())
    counts = list(aa_counts.values())

    ax.bar(aa_list, counts)
    ax.set_xlabel("Amino Acid")
    ax.set_ylabel("Frequency")
    ax.set_title("Histogram of Amino Acids")

    plt.xticks(rotation=45, ha='right')

    # Spremanje slike u memoriju (kao bytes)
    buf = io.BytesIO()
    plt.tight_layout()
    plt.savefig(buf, format='png')
    buf.seek(0)

    # Pretvaranje slike u base64 string
    image_base64 = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    plt.close(fig)

    return image_base64
