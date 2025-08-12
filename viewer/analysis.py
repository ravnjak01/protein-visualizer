from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.vectors import calc_dihedral
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')  # <-- OBAVEZNO pre plt importa
import matplotlib.pyplot as plt
import io
import base64
import numpy as np
import math

def get_amino_acid_letter(three_letter_code):
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return aa_dict.get(three_letter_code, 'X')

def calculate_phi_psi_angles(residues):
    """Kalkuliše phi i psi uglove za svaki ostatak"""
    phi_psi_angles = []
    
    for i in range(1, len(residues) - 1):  # Ne možemo izračunati za prvi i poslednji ostatak
        try:
            # Potrebni atomi za phi ugao: C(i-1) - N(i) - CA(i) - C(i)
            prev_c = residues[i-1]['C']
            curr_n = residues[i]['N']
            curr_ca = residues[i]['CA']
            curr_c = residues[i]['C']
            
            # Potrebni atomi za psi ugao: N(i) - CA(i) - C(i) - N(i+1)
            next_n = residues[i+1]['N']
            
            # Kalkulacija phi i psi uglova
            phi = calc_dihedral(prev_c.get_vector(), curr_n.get_vector(), 
                              curr_ca.get_vector(), curr_c.get_vector())
            psi = calc_dihedral(curr_n.get_vector(), curr_ca.get_vector(), 
                              curr_c.get_vector(), next_n.get_vector())
            
            # Konverzija iz radijana u stepene
            phi_deg = math.degrees(phi)
            psi_deg = math.degrees(psi)
            
            phi_psi_angles.append((phi_deg, psi_deg, i))
            
        except KeyError as e:
            # Nema potrebne atome
            phi_psi_angles.append((None, None, i))
            
    return phi_psi_angles

def assign_secondary_structure_simple(phi_psi_angles):
    """
    Jednostavna doela sekundarne strukture na osnovu phi/psi uglova
    H = Alpha helix
    E = Beta sheet
    C = Coil/loop
    """
    ss_assignment = []
    
    for phi, psi, res_index in phi_psi_angles:
        if phi is None or psi is None:
            ss_assignment.append('C')
            continue
            
        # Alpha helix region: phi ~ -60, psi ~ -45
        if (-90 <= phi <= -30) and (-70 <= psi <= -20):
            ss_assignment.append('H')
        # Beta sheet region: phi ~ -120, psi ~ +120
        elif (-180 <= phi <= -90) and (90 <= psi <= 180):
            ss_assignment.append('E')
        # Dodatni beta sheet region
        elif (-180 <= phi <= -90) and (-180 <= psi <= -90):
            ss_assignment.append('E')
        else:
            ss_assignment.append('C')
            
    return ss_assignment

def hydrogen_bond_analysis(residues):
    """Jednostavna analiza vodoničnih veza za identifikaciju sekundarne strukture"""
    hbonds = []
    
    for i in range(len(residues)):
        for j in range(i + 3, len(residues)):  # Minimum 3 ostatka razmaka
            try:
                # Donor: N-H iz ostatka j
                # Akceptor: O=C iz ostatka i
                
                if 'N' in residues[j] and 'O' in residues[i]:
                    n_coord = residues[j]['N'].get_coord()
                    o_coord = residues[i]['O'].get_coord()
                    
                    # Kalkulacija rastojanja
                    distance = np.linalg.norm(n_coord - o_coord)
                    
                    # Hbond cutoff ~ 3.5 Å
                    if distance <= 3.5:
                        hbonds.append((i, j, distance))
                        
            except KeyError:
                continue
                
    return hbonds

def refine_secondary_structure(ss_simple, hbonds, residues):
    """Poboljšava dodeljivanje sekundarne strukture koristeći hydrogen bond patterns"""
    ss_refined = ss_simple.copy()
    
    # Alpha helix pattern: i to i+4 hydrogen bonds
    alpha_helix_regions = set()
    for donor, acceptor, dist in hbonds:
        if acceptor - donor == 4:  # Alpha helix pattern
            for k in range(donor, acceptor + 1):
                alpha_helix_regions.add(k)
    
    # Beta sheet pattern: longer range hydrogen bonds
    beta_sheet_regions = set()
    for donor, acceptor, dist in hbonds:
        if acceptor - donor > 4:  # Potential beta sheet
            beta_sheet_regions.add(donor)
            beta_sheet_regions.add(acceptor)
    
    # Update assignments
    for i, _ in enumerate(ss_refined):
        actual_res_idx = i + 1  # Jer počinjemo od drugog ostatka
        if actual_res_idx in alpha_helix_regions:
            ss_refined[i] = 'H'
        elif actual_res_idx in beta_sheet_regions and ss_refined[i] != 'H':
            ss_refined[i] = 'E'
    
    return ss_refined

def smooth_secondary_structure(ss_assignment, min_length=3):
    """Gladi dodeljivanje sekundarne strukture - uklanja kratke segmente"""
    smoothed = ss_assignment.copy()
    
    # Ukloni kratke helix i sheet segmente
    i = 0
    while i < len(smoothed):
        if smoothed[i] in ['H', 'E']:
            # Nađi kraj ovog segmenta
            j = i
            while j < len(smoothed) and smoothed[j] == smoothed[i]:
                j += 1
            
            # Ako je segment kraći od min_length, promeni u coil
            if j - i < min_length:
                for k in range(i, j):
                    smoothed[k] = 'C'
            i = j
        else:
            i += 1
    
    return smoothed

def calculate_secondary_structure_stats(ss_assignment):
    """Kalkuliše statistike sekundarne strukture"""
    total = len(ss_assignment)
    if total == 0:
        return {}
    
    helix_count = ss_assignment.count('H')
    sheet_count = ss_assignment.count('E')
    coil_count = ss_assignment.count('C')
    
    return {
        'helix_percentage': (helix_count / total) * 100,
        'sheet_percentage': (sheet_count / total) * 100,
        'coil_percentage': (coil_count / total) * 100,
        'helix_count': helix_count,
        'sheet_count': sheet_count,
        'coil_count': coil_count,
        'total': total
    }

def analyze_pdb_with_secondary_structure(file_path):
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
            
            # Kalkulacija sekundarne strukture
            phi_psi_angles = calculate_phi_psi_angles(residues)
            ss_simple = assign_secondary_structure_simple(phi_psi_angles)
            
            # Hydrogen bond analiza za poboljšanje
            hbonds = hydrogen_bond_analysis(residues)
            ss_refined = refine_secondary_structure(ss_simple, hbonds, residues)
            
            # Glaćenje
            ss_final = smooth_secondary_structure(ss_refined)
            
            # Statistike sekundarne strukture
            ss_stats = calculate_secondary_structure_stats(ss_final)

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
                "atom_count": atom_count,
                "secondary_structure": ''.join(ss_final),
                "ss_stats": ss_stats,
                "phi_psi_angles": phi_psi_angles,
                "hydrogen_bonds": hbonds
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

def create_secondary_structure_plot(chain_data):
    """Kreira vizualizaciju sekundarne strukture za sve lance"""
    fig, axes = plt.subplots(len(chain_data), 1, figsize=(12, 2 * len(chain_data)))
    if len(chain_data) == 1:
        axes = [axes]
    
    colors = {'H': 'red', 'E': 'blue', 'C': 'gray'}
    
    for idx, (chain_id, data) in enumerate(chain_data.items()):
        ax = axes[idx]
        ss = data['secondary_structure']
        
        # Kreiraj x i y koordinate za crtanje
        x_coords = list(range(len(ss)))
        y_coords = [0] * len(ss)  # Sve na istoj visini
        
        # Oboji prema tipu sekundarne strukture
        for i, ss_type in enumerate(ss):
            color = colors.get(ss_type, 'black')
            ax.scatter(x_coords[i], y_coords[i], c=color, s=50, alpha=0.8)
        
        ax.set_title(f'Secondary Structure - Chain {chain_id}')
        ax.set_xlabel('Residue Position')
        ax.set_ylabel('')
        ax.set_ylim(-0.5, 0.5)
        
        # Dodaj legendu
        legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                    markerfacecolor=colors[ss_type], markersize=10,
                                    label=f'{ss_type} ({ss_type_name})')
                         for ss_type, ss_type_name in [('H', 'Helix'), ('E', 'Sheet'), ('C', 'Coil')]]
        ax.legend(handles=legend_elements, loc='upper right')
        
        # Dodaj statistike kao tekst
        stats = data['ss_stats']
        stats_text = f"H: {stats['helix_percentage']:.1f}% | E: {stats['sheet_percentage']:.1f}% | C: {stats['coil_percentage']:.1f}%"
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    # Spremanje slike u memoriju
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=150)
    buf.seek(0)
    
    image_base64 = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    plt.close(fig)
    
    return image_base64

def create_ramachandran_plot(chain_data):
    """Kreira Ramachandran plot za sve lance"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    colors = {'H': 'red', 'E': 'blue', 'C': 'gray'}
    
    for chain_id, data in chain_data.items():
        phi_psi_angles = data['phi_psi_angles']
        ss = data['secondary_structure']
        
        for i, (phi, psi, _) in enumerate(phi_psi_angles):
            if phi is not None and psi is not None:
                ss_type = ss[i] if i < len(ss) else 'C'
                color = colors.get(ss_type, 'black')
                ax.scatter(phi, psi, c=color, s=20, alpha=0.6)
    
    ax.set_xlabel('Phi (degrees)')
    ax.set_ylabel('Psi (degrees)')
    ax.set_title('Ramachandran Plot')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.grid(True, alpha=0.3)
    
    # Dodaj legendu
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                markerfacecolor=colors[ss_type], markersize=8,
                                label=f'{ss_type} ({ss_type_name})')
                     for ss_type, ss_type_name in [('H', 'Helix'), ('E', 'Sheet'), ('C', 'Coil')]]
    ax.legend(handles=legend_elements)
    
    plt.tight_layout()
    
    # Spremanje slike u memoriju
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=150)
    buf.seek(0)
    
    image_base64 = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    plt.close(fig)
    
    return image_base64

# Primer korišćenja
if __name__ == "__main__":
    # Zameniti sa putanjom do vašeg PDB fajla
    pdb_file = "primer.pdb"
    
    try:
        results = analyze_pdb_with_secondary_structure(pdb_file)
        
        # Ispiši rezultate
        print("=== ANALIZA AMINOKISELINA ===")
        for aa, count in results["aa_counts"].items():
            print(f"{aa}: {count}")
        
        print("\n=== ANALIZA PO LANCIMA ===")
        for chain_id, data in results["chains"].items():
            print(f"\nLanac {chain_id}:")
            print(f"  Broj ostataka: {data['residue_count']}")
            print(f"  Broj atoma: {data['atom_count']}")
            print(f"  Sekvenca: {data['sequence'][:50]}{'...' if len(data['sequence']) > 50 else ''}")
            print(f"  Sekundarna struktura: {data['secondary_structure'][:50]}{'...' if len(data['secondary_structure']) > 50 else ''}")
            
            print(f"  Statistike sekundarne strukture:")
            stats = data['ss_stats']
            print(f"    Alpha helix: {stats['helix_percentage']:.1f}% ({stats['helix_count']} ostataka)")
            print(f"    Beta sheet: {stats['sheet_percentage']:.1f}% ({stats['sheet_count']} ostataka)")
            print(f"    Coil/loop: {stats['coil_percentage']:.1f}% ({stats['coil_count']} ostataka)")
        
        # Kreiranje grafika
        aa_histogram = create_histogram(results["aa_counts"])
        ss_plot = create_secondary_structure_plot(results["chains"])
        ramachandran = create_ramachandran_plot(results["chains"])
        
        print("\n=== GRAFICI KREIRANI ===")
        print("- Histogram aminokiselina")
        print("- Prikaz sekundarne strukture")
        print("- Ramachandran plot")
        
    except Exception as e:
        print(f"Greška pri analizi: {e}")

        # Dodajte ove funkcije u analysis.py

def create_text_based_secondary_structure(chain_data):
    """
    Kreira tekstualni prikaz sekundarne strukture koji je lakši za čitanje
    """
    results = {}
    
    for chain_id, data in chain_data.items():
        sequence = data.get('sequence', '')
        ss_structure = data.get('secondary_structure', '')
        
        if not sequence or not ss_structure:
            continue
            
        # Kreiraj formatiran prikaz
        formatted_display = format_ss_for_display(sequence, ss_structure)
        
        # Identifikuj segmente
        segments = identify_ss_segments(ss_structure)
        
        results[chain_id] = {
            'formatted_display': formatted_display,
            'segments': segments,
            'sequence_length': len(sequence)
        }
    
    return results

def format_ss_for_display(sequence, ss_structure, line_length=60):
    """
    Formatira sekundarnu strukturu u blokove za lakše čitanje
    """
    if len(sequence) == 0:
        return []
    
    # Dodaj prvi i poslednji ostatak kao Coil
    full_ss = 'C' + ss_structure + 'C'
    
    formatted_blocks = []
    
    for i in range(0, len(sequence), line_length):
        end_pos = min(i + line_length, len(sequence))
        
        # Sekvenca aminokiselina
        seq_block = sequence[i:end_pos]
        
        # Sekundarna struktura
        ss_block = full_ss[i:end_pos]
        
        # Pozicije (1-indexed)
        start_pos = i + 1
        end_pos_actual = i + len(seq_block)
        
        # Simboli za vizualni prikaz
        visual_block = ''
        for ss_char in ss_block:
            if ss_char == 'H':
                visual_block += '█'  # Pun blok za helix
            elif ss_char == 'E':
                visual_block += '▬'  # Linija za sheet
            else:
                visual_block += '·'  # Tačka za coil
        
        formatted_blocks.append({
            'start_pos': start_pos,
            'end_pos': end_pos_actual,
            'sequence': seq_block,
            'secondary_structure': ss_block,
            'visual': visual_block,
            'length': len(seq_block)
        })
    
    return formatted_blocks

def identify_ss_segments(ss_structure):
    """
    Identifikuje kontinuirane segmente sekundarne strukture
    """
    if not ss_structure:
        return []
    
    # Dodaj C na početak i kraj
    full_ss = 'C' + ss_structure + 'C'
    
    segments = []
    current_type = None
    current_start = None
    
    for i, ss_type in enumerate(full_ss):
        if ss_type != current_type:
            # Završi prethodna segment
            if current_type is not None and current_start is not None:
                segments.append({
                    'type': current_type,
                    'start': current_start + 1,  # 1-indexed
                    'end': i,  # 1-indexed
                    'length': i - current_start,
                    'type_name': get_ss_type_name(current_type)
                })
            
            # Počni novi segment
            current_type = ss_type
            current_start = i
    
    # Dodaj poslednji segment
    if current_type is not None and current_start is not None:
        segments.append({
            'type': current_type,
            'start': current_start + 1,
            'end': len(full_ss),
            'length': len(full_ss) - current_start,
            'type_name': get_ss_type_name(current_type)
        })
    
    # Filtriraj kratke segmente (< 3 ostatka) osim ako nisu sheet ili helix
    significant_segments = []
    for seg in segments:
        if seg['length'] >= 3 or seg['type'] in ['H', 'E']:
            significant_segments.append(seg)
    
    return significant_segments

def get_ss_type_name(ss_type):
    """Vraća puno ime tipa sekundarne strukture"""
    names = {
        'H': 'Alpha Helix',
        'E': 'Beta Sheet', 
        'C': 'Coil/Loop'
    }
    return names.get(ss_type, 'Unknown')

def create_interactive_ss_display(chain_data):
    """
    Kreira podatke za interaktivni prikaz sekundarne strukture
    """
    results = {}
    
    for chain_id, data in chain_data.items():
        sequence = data.get('sequence', '')
        ss_structure = data.get('secondary_structure', '')
        
        if not sequence or not ss_structure:
            continue
        
        # Kreiraj podatke za svaki ostatak
        residue_data = []
        full_ss = 'C' + ss_structure + 'C'
        
        for i, (aa, ss) in enumerate(zip(sequence, full_ss)):
            residue_data.append({
                'position': i + 1,
                'amino_acid': aa,
                'aa_full_name': get_aa_full_name(aa),
                'ss_type': ss,
                'ss_name': get_ss_type_name(ss),
                'css_class': get_ss_css_class(ss)
            })
        
        results[chain_id] = {
            'residues': residue_data,
            'total_residues': len(residue_data)
        }
    
    return results

def get_aa_full_name(aa_code):
    """Vraća puno ime aminokiseline"""
    aa_names = {
        'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic acid',
        'C': 'Cysteine', 'E': 'Glutamic acid', 'Q': 'Glutamine', 'G': 'Glycine',
        'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
        'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
        'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine'
    }
    return aa_names.get(aa_code, 'Unknown')

def get_ss_css_class(ss_type):
    """Vraća CSS klasu za bojenje"""
    classes = {
        'H': 'ss-helix',
        'E': 'ss-sheet',
        'C': 'ss-coil'
    }
    return classes.get(ss_type, 'ss-unknown')

def create_summary_table(chain_data):
    """
    Kreira tabelarni pregled sekundarne strukture
    """
    summary_data = []
    
    for chain_id, data in chain_data.items():
        ss_stats = data.get('ss_stats', {})
        segments = identify_ss_segments(data.get('secondary_structure', ''))
        
        # Brojanje segmenata po tipu
        helix_segments = [s for s in segments if s['type'] == 'H']
        sheet_segments = [s for s in segments if s['type'] == 'E']
        
        summary_data.append({
            'chain_id': chain_id,
            'total_residues': ss_stats.get('total', 0),
            'helix_percentage': ss_stats.get('helix_percentage', 0),
            'sheet_percentage': ss_stats.get('sheet_percentage', 0),
            'coil_percentage': ss_stats.get('coil_percentage', 0),
            'num_helix_segments': len(helix_segments),
            'num_sheet_segments': len(sheet_segments),
            'avg_helix_length': sum(h['length'] for h in helix_segments) / len(helix_segments) if helix_segments else 0,
            'avg_sheet_length': sum(s['length'] for s in sheet_segments) / len(sheet_segments) if sheet_segments else 0,
            'longest_helix': max((h['length'] for h in helix_segments), default=0),
            'longest_sheet': max((s['length'] for s in sheet_segments), default=0)
        })
    
    return summary_data

# Dodajte ove funkcije u views.py

def format_secondary_structure_display_enhanced(sequence, ss_structure):
    """
    Poboljšana verzija formatiranja sekundarne strukture
    """
    if not sequence or not ss_structure:
        return {
            'text_display': [],
            'segments': [],
            'interactive_data': [],
            'summary': {}
        }
    
    chain_data = {'temp': {
        'sequence': sequence,
        'secondary_structure': ss_structure,
        'ss_stats': calculate_secondary_structure_stats(list(ss_structure))
    }}
    
    # Različiti formati prikaza
    text_display = create_text_based_secondary_structure(chain_data)
    interactive_data = create_interactive_ss_display(chain_data)
    
    return {
        'text_display': text_display.get('temp', {}).get('formatted_display', []),
        'segments': text_display.get('temp', {}).get('segments', []),
        'interactive_data': interactive_data.get('temp', {}).get('residues', []),
        'summary': chain_data['temp']['ss_stats']
    }

def calculate_secondary_structure_stats(ss_assignment):
    """Kalkuliše statistike sekundarne strukture"""
    total = len(ss_assignment)
    if total == 0:
        return {}
    
    helix_count = ss_assignment.count('H')
    sheet_count = ss_assignment.count('E')
    coil_count = ss_assignment.count('C')
    
    return {
        'helix_percentage': (helix_count / total) * 100,
        'sheet_percentage': (sheet_count / total) * 100,
        'coil_percentage': (coil_count / total) * 100,
        'helix_count': helix_count,
        'sheet_count': sheet_count,
        'coil_count': coil_count,
        'total': total
    }