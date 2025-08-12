from django.shortcuts import render
from .forms import ProteinFileForm
from .analysis import analyze_pdb
from django.shortcuts  import  render
from .analysis import create_histogram
from viewer import analysis
def upload_protein(request):
   
    context = {
   
    }
    
    if request.method == 'POST':
        form = ProteinFileForm(request.POST, request.FILES)
        if form.is_valid():
            protein = form.save()

            try:
                # Analiza fajla
                results = analyze_pdb(protein.file.path)
                
                
                
                if results and isinstance(results, dict):
                    # Dodaj rezultate u context
                    context.update(results)
                    context['uploaded_file'] = protein
                  
                    if 'chains' in results:
                       chain_data = []
                    for chain_id, data in results['chains'].items():
                      chain_data.append({
                          "chain_id": chain_id,
                    "num_residues": data["residue_count"],
                    "atom_count": data["atom_count"],
                    "sequence": data["sequence"],
                    "aa_counts": data["aa_counts"]
                      })

                    context['chain_data'] = chain_data  
                    image_data = create_histogram(results['aa_counts'])
                    context['histogram'] = image_data
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

