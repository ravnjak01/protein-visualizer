# Protein Structure Analysis from PDB Files

##  Overview
This project provides tools for **comprehensive analysis of protein structures** stored in **PDB (Protein Data Bank)** format.  
It focuses on **structural, chain-specific, and amino acid composition analysis**, along with **3D interactive visualization** of proteins.


---

## Features
- **Load and parse PDB files** using [Biopython](https://biopython.org/).
- **Detailed per-chain analysis**:
  - Chain ID
  - Number of residues
  - Amino acid composition
- **Histogram of amino acids** for each chain.
- **3D interactive visualization** of protein structures .
- **Secondary structure extraction** (helix, sheet, coil) directly from the PDB file.
- **Summary table** showing structural statistics.

---

 ## Project Structure
```
PROTEIN-VISUALIZER
├─── .gitignore
├─── db.sqlite3
├─── manage.py
├─── protein_viewer_db.sqlite3
├─── requirements.txt
├─── media/
├─── protein_viewer/
│   ├─── __pycache__/
│   ├─── __init__.py
│   ├─── asgi.py
│   ├─── settings.py
│   ├─── urls.py
│   └─── wsgi.py
├─── venv/
├─── viewer/
│   ├─── __pycache__/
│   ├─── migrations/
│   ├─── templates/
│   │   └─── viewer/
│   │       └─── upload.html
│   ├─── __init__.py
│   ├─── admin.py
│   ├─── analysis.py
│   ├─── apps.py
│   ├─── forms.py
│   ├─── models.py
│   ├─── tests.py
│   └─── urls.py
```
## Technologies used
-- **Django**
- **Python**
- **SQLite**
- **HTML,CSS,JS,NGL.JS** 
- **Biopython**
- **Matplotlib**
- **venv**
- **Git**
- **Numpy**
- **Pandas**
- **Pillow**

## Limitations
Works only with valid .pdb files.
Secondary structure extraction relies on existing PDB annotations.
No prediction of unknown structures — only analysis of provided PDB files.



