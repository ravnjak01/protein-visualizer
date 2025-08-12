from django.shortcuts import render
from .forms import ProteinFileForm
from django.shortcuts import render
from viewer import analysis

# Bezbedni import - radi sa postojećim kodom i omogućava dodavanje novih funkcija
try:
    # Pokušaj da importuješ nove funkcije
    from .analysis import analyze_pdb_with_secondary_structure, create_secondary_structure_plot, create_ramachandran_plot
    HAS_SECONDARY_STRUCTURE = True
    print("Sekundarna struktura je dostupna")
except ImportError as e:
    print(f"Sekundarna struktura nije dostupna: {e}")
    HAS_SECONDARY_STRUCTURE = False

# Ovi importi uvek rade jer su postojeći
from .analysis import analyze_pdb_with_secondary_structure,create_histogram

def upload_protein(request):
    context = {}
    
    if request.method == 'POST':
        form = ProteinFileForm(request.POST, request.FILES)
        if form.is_valid():
            protein = form.save()

            try:
                # Koristi proširenu analizu ako je dostupna, inače osnovnu
                if HAS_SECONDARY_STRUCTURE:
                    print("Koristim proširenu analizu...")
                    results = analyze_pdb_with_secondary_structure(protein.file.path)
                else:
                    print("Koristim osnovnu analizu...")
                    results = analyze_pdb(protein.file.path)
                
                print(f"Rezultati analize: {type(results)}")
                if results:
                    print(f"Ključevi u results: {results.keys()}")
                
                if results and isinstance(results, dict):
                    # Dodaj osnovne rezultate u context
                    context.update(results)
                    context['uploaded_file'] = protein
                  
                    if 'chains' in results:
                        chain_data = []
                        
                        for chain_id, data in results['chains'].items():
                            print(f"Obrađujem lanac {chain_id}")
                            print(f"Podaci za lanac: {list(data.keys())}")
                            
                            chain_info = {
                                "chain_id": chain_id,
                                "num_residues": data["residue_count"],
                                "atom_count": data["atom_count"],
                                "sequence": data["sequence"],
                                "aa_counts": data["aa_counts"],
                            }
                            
                            # Dodaj podatke o sekundarnoj strukturi ako su dostupni
                            if HAS_SECONDARY_STRUCTURE and "secondary_structure" in data:
                                ss_structure = data.get("secondary_structure", "")
                                sequence = data.get("sequence", "")
                                
                                print(f"SS struktura za lanac {chain_id}: '{ss_structure[:50]}...'")
                                print(f"Dužina sekvence: {len(sequence)}, Dužina SS: {len(ss_structure)}")
                                
                                if ss_structure and sequence:
                                    # Kreiraj formatiran prikaz
                                    formatted_ss = create_formatted_ss_display(sequence, ss_structure)
                                    
                                    chain_info.update({
                                        "secondary_structure": ss_structure,
                                        "ss_stats": data.get("ss_stats", {}),
                                        "num_hydrogen_bonds": len(data.get("hydrogen_bonds", [])),
                                        "formatted_ss": formatted_ss
                                    })
                                    
                                    print(f"Formatted SS kreiran za lanac {chain_id}")
                                else:
                                    print(f"Nema SS podataka za lanac {chain_id}")
                                    chain_info["formatted_ss"] = None
                            else:
                                print(f"SS nije dostupna za lanac {chain_id}")
                                chain_info["formatted_ss"] = None
                            
                            chain_data.append(chain_info)

                        context['chain_data'] = chain_data
                        context['has_secondary_structure'] = HAS_SECONDARY_STRUCTURE
                        
                        print(f"Ukupno lancova: {len(chain_data)}")
                        for i, chain in enumerate(chain_data):
                            print(f"Lanac {i}: {chain['chain_id']}, ima SS: {chain.get('formatted_ss') is not None}")
                        
                        # Kreiraj grafike
                        try:
                            # Osnovni histogram aminokiselina (uvek radi)
                            aa_histogram = create_histogram(results['aa_counts'])
                            context['histogram'] = aa_histogram
                            print("Histogram kreiran")
                            
                            # Dodatni grafici samo ako je moguće
                            if HAS_SECONDARY_STRUCTURE:
                                try:
                                    # Graf sekundarne strukture
                                    ss_plot = create_secondary_structure_plot(results['chains'])
                                    context['secondary_structure_plot'] = ss_plot
                                    print("SS plot kreiran")
                                    
                                    # Ramachandran plot
                                    ramachandran_plot = create_ramachandran_plot(results['chains'])
                                    context['ramachandran_plot'] = ramachandran_plot
                                    print("Ramachandran plot kreiran")
                                    
                                    # Summary statistike
                                    context['summary_stats'] = calculate_summary_statistics(results['chains'])
                                    print("Summary statistike kreirane")
                                    
                                except Exception as plot_error:
                                    print(f"Greška pri kreiranju SS grafika: {plot_error}")
                                    import traceback
                                    traceback.print_exc()
                                    context['plot_error'] = f"Greška pri kreiranju SS grafika: {str(plot_error)}"
                            
                        except Exception as plot_error:
                            print(f"Greška pri kreiranju grafika: {plot_error}")
                            import traceback
                            traceback.print_exc()
                            context['plot_error'] = f"Greška pri kreiranju grafika: {str(plot_error)}"
                        
                else:
                    context['error'] = "Greška pri analizi PDB fajla"
                    
            except Exception as e:
                print(f"Exception in upload_protein: {e}")
                import traceback
                traceback.print_exc()
                context['error'] = f"Greška pri obradi fajla: {str(e)}"
                
    else:
        form = ProteinFileForm()

    context['form'] = form
    return render(request, 'viewer/upload.html', context)

def create_formatted_ss_display(sequence, ss_structure):
    """
    Jednostavna funkcija za kreiranje formatiranog prikaza sekundarne strukture
    """
    try:
        print(f"Kreiram formatiran prikaz - seq len: {len(sequence)}, ss len: {len(ss_structure)}")
        
        if not sequence or not ss_structure:
            return {
                'text_display': [],
                'segments': [],
                'interactive_data': [],
                'summary': {}
            }
        
        # Dodaj C na početak i kraj za potpunu sekundarnu strukturu
        full_ss = 'C' + ss_structure + 'C'
        
        # Tekstualni prikaz po blokovima od 60 karaktera
        text_display = []
        line_length = 60
        
        for i in range(0, len(sequence), line_length):
            end_pos = min(i + line_length, len(sequence))
            
            seq_block = sequence[i:end_pos]
            ss_block = full_ss[i:end_pos]
            
            # Vizualni prikaz
            visual_block = ''
            for ss_char in ss_block:
                if ss_char == 'H':
                    visual_block += '█'  # Helix
                elif ss_char == 'E':
                    visual_block += '▬'  # Sheet
                else:
                    visual_block += '·'  # Coil
            
            text_display.append({
                'start_pos': i + 1,
                'end_pos': end_pos,
                'sequence': seq_block,
                'secondary_structure': ss_block,
                'visual': visual_block,
                'length': len(seq_block)
            })
        
        # Interaktivni prikaz - podatak po ostatku
        interactive_data = []
        ss_names = {'H': 'Helix', 'E': 'Sheet', 'C': 'Coil'}
        
        for i, (aa, ss) in enumerate(zip(sequence, full_ss)):
            interactive_data.append({
                'position': i + 1,
                'amino_acid': aa,
                'aa_full_name': get_aa_full_name(aa),
                'ss_type': ss,
                'ss_name': ss_names.get(ss, 'Unknown'),
                'css_class': get_ss_css_class(ss)
            })
        
        # Identifikuj segmente
        segments = []
        current_type = None
        current_start = None
        
        for i, ss_type in enumerate(full_ss):
            if ss_type != current_type:
                # Završi prethodnji segment
                if current_type is not None and current_start is not None:
                    length = i - current_start
                    if length >= 3:  # Samo značajne segmente
                        segments.append({
                            'type': current_type,
                            'start': current_start + 1,
                            'end': i,
                            'length': length,
                            'type_name': ss_names.get(current_type, 'Unknown')
                        })
                
                current_type = ss_type
                current_start = i
        
        # Dodaj poslednji segment
        if current_type is not None and current_start is not None:
            length = len(full_ss) - current_start
            if length >= 3:
                segments.append({
                    'type': current_type,
                    'start': current_start + 1,
                    'end': len(full_ss),
                    'length': length,
                    'type_name': ss_names.get(current_type, 'Unknown')
                })
        
        # Statistike
        total = len(full_ss)
        helix_count = full_ss.count('H')
        sheet_count = full_ss.count('E')
        coil_count = full_ss.count('C')
        
        summary = {
            'helix_percentage': (helix_count / total) * 100 if total > 0 else 0,
            'sheet_percentage': (sheet_count / total) * 100 if total > 0 else 0,
            'coil_percentage': (coil_count / total) * 100 if total > 0 else 0,
            'helix_count': helix_count,
            'sheet_count': sheet_count,
            'coil_count': coil_count,
            'total': total
        }
        
        result = {
            'text_display': text_display,
            'segments': segments,
            'interactive_data': interactive_data,
            'summary': summary
        }
        
        print(f"Kreiran formatiran prikaz: {len(text_display)} blokova, {len(segments)} segmenata, {len(interactive_data)} ostataka")
        return result
        
    except Exception as e:
        print(f"Greška u create_formatted_ss_display: {e}")
        import traceback
        traceback.print_exc()
        return {
            'text_display': [],
            'segments': [],
            'interactive_data': [],
            'summary': {}
        }

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

def calculate_summary_statistics(chains_data):
    """Kalkuliše ukupne statistike za sve lance"""
    if not HAS_SECONDARY_STRUCTURE:
        return {}
        
    total_residues = 0
    total_helix = 0
    total_sheet = 0
    total_coil = 0
    total_hbonds = 0
    
    for chain_id, data in chains_data.items():
        ss_stats = data.get('ss_stats', {})
        total_residues += ss_stats.get('total', 0)
        total_helix += ss_stats.get('helix_count', 0)
        total_sheet += ss_stats.get('sheet_count', 0)
        total_coil += ss_stats.get('coil_count', 0)
        total_hbonds += len(data.get('hydrogen_bonds', []))
    
    if total_residues == 0:
        return {}
        
    return {
        'total_residues': total_residues,
        'total_chains': len(chains_data),
        'helix_percentage': (total_helix / total_residues) * 100 if total_residues > 0 else 0,
        'sheet_percentage': (total_sheet / total_residues) * 100 if total_residues > 0 else 0,
        'coil_percentage': (total_coil / total_residues) * 100 if total_residues > 0 else 0,
        'total_hydrogen_bonds': total_hbonds,
        'avg_hbonds_per_residue': total_hbonds / total_residues if total_residues > 0 else 0
    }