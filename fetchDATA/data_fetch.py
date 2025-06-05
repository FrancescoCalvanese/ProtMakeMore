import os
import requests
from Bio.PDB import PDBParser, DSSP



def parse_identifier(identifier):
    """Extract UniProt ID and residue range."""
    identifier = identifier.strip(">").strip()
    base, rng = identifier.split('/')
    uniprot_id = base.split('_')[0]
    start, end = map(int, rng.split('-'))
    return uniprot_id, start, end

def fetch_uniprot_sequence(uniprot_id):
    """Fetch the full sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Could not fetch sequence for UniProt ID {uniprot_id}")
        return "token"
    fasta = response.text
    sequence = ''.join(fasta.split('\n')[1:])
    return sequence

def remove_dbref_lines_from_file(file_path):
    """
    Removes lines starting with 'DBREF' from a PDB file in-place.

    Parameters:
    - file_path: str, path to the PDB file to clean.
    """
    with open(file_path, 'r') as fin:
        lines = fin.readlines()

    with open(file_path, 'w') as fout:
        for line in lines:
            if not line.startswith("DBREF"):
                fout.write(line)

def download_alphafold_pdb(uniprot_id, out_dir="pdbs"):
    if uniprot_id == "token":
        return "token"
    """Download the AlphaFold PDB model."""
    os.makedirs(out_dir, exist_ok=True)
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        out_path = os.path.join(out_dir, f"{uniprot_id}_AF.pdb")
        with open(out_path, 'w') as f:
            f.write(response.text)
            print(f"Available AlphaFold model for {uniprot_id} to {out_path}")
        return out_path
    else:
        print(f"Not Available AlphaFold model for {uniprot_id}")
        return "token"


def extract_secondary_structure(pdb_path, start, end):
    if pdb_path == "token":
        return "???", "???"
    remove_dbref_lines_from_file(pdb_path)
    """Use DSSP to extract secondary structure assignments."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("AF", pdb_path)
    model = structure[0]
    dssp = DSSP(model, pdb_path)

    secstr = ""
    sequence = ""
    
    for key in dssp.keys():
        chain_id, res_id = key
        resnum = res_id[1]
        if start <= resnum <= end:
            aa, ss = dssp[key][1], dssp[key][2]
            sequence += aa
            secstr += ss if ss != ' ' else 'C'  # Convert blanks to coil

    os.remove(pdb_path)
    return sequence, secstr

def data_fetch(identifier, tk = 1):
  if tk == 1:
    uniprot_id, start, end = parse_identifier(identifier)
    #print(f"UniProt ID: {uniprot_id}, Range: {start}-{end}")

    full_seq = fetch_uniprot_sequence(uniprot_id)
    sub_seq = full_seq[start-1:end]
    #print(f"\nSubsequence from UniProt [{start}-{end}]:\n{sub_seq}")

    pdb_path = download_alphafold_pdb(uniprot_id)
    #print(f"PDB downloaded to: {pdb_path}")

    ss_seq, ss_struct = extract_secondary_structure(pdb_path, start, end)
    #print(f"\nSecondary structure from PDB [{start}-{end}]:\n{ss_struct}")

    return sub_seq, ss_struct
  else:
    uniprot_id, start, end = parse_identifier(identifier)
    #print(f"UniProt ID: {uniprot_id}, Range: {start}-{end}")

    full_seq = fetch_uniprot_sequence(uniprot_id)
    sub_seq = full_seq[start-1:end]
    #print(f"\nSubsequence from UniProt [{start}-{end}]:\n{sub_seq}")
    return sub_seq

def extract_gt_lines(file_path):
    """
    Reads a file and returns a list of lines that start with '>'.

    Parameters:
    - file_path: str, path to the input file

    Returns:
    - list_gt: list of strings, each starting with '>'
    """
    list_gt = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                list_gt.append(line.strip())  # remove trailing newline
    return list_gt


def simplify_secondary_structure(dssp_str):
    mapping = {
        'H': 'H', 'G': 'H', 'I': 'H',  # helix types
        'E': 'E', 'B': 'E',    
        '?' : '?',        # strand types
        'T': 'C', 'S': 'C', '-': 'C', 'P': 'C'  # coil/irregular
    }
    return ''.join(mapping.get(c, c) for c in dssp_str)