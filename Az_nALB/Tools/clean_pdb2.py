
from __future__ import print_function
''' This script cleans PDBs for Rosetta by removing extraneous information, converting residues names and renumbering.

It outputs both the cleaned PDB and a fasta file of the cleaned sequence.

Required parameters are the name of the PDB you want to clean, and the chain ids of the chains you want.

The PDB name may be specified with or without the .pdb file handle and may be provided as a gziped file.
If the PDB isn't found locally, the given 4 letter code will be fetched from the internet.

Chain id: only the specified chains will be extracted. You may specify more than one: "AB" gets you chain A and B,
and "C" gets you just chain C. Special notations are "nochain" to remove chain identiry from the output, and "ignorechain"
to get all the chains.

(Script written by Phil Bradley, Rhiju Das, Michael Tyka, TJ Brunette, and James Thompson from the Baker Lab. Edits done by Steven Combs, Sam Deluca, Jordan Willis and Rocco Moretti from the Meiler Lab.)
'''

# Function of this script: "clean" raw pdb file by following tasks so that rosetta modeling becomes easier

## starts residue number at 1
## translates certain residues to their cannonical amino acid equivalents
## removes unknown residues
## removes residues with 0 occupancy
## generates a fasta file
## and leaves the 1st model among many NMR models

#!/usr/bin/env python
'''
This script cleans PDBs for Rosetta by:
- Renumbering residues starting from 1
- Converting non-standard residues to canonical equivalents
- Removing unknown residues and residues with 0 occupancy
- Generating FASTA files
- Keeping only the first NMR model

Batch processing mode added to handle multiple PDB files at once.

Single: python clean_pdb.py input.pdb A
To use in batch mode: python clean_pdb.py --batch --input-dir ./raw_pdbs --output-dir ./clean A

'''
import sys
import os
import glob
from os import popen, system, unlink, makedirs
from os.path import exists, basename, join, abspath, splitext
from optparse import OptionParser

# Local package imports
from amino_acids import longer_names
from amino_acids import modres

# Remote host for downloading PDBs
remote_host = ''

# Global variables (reset for each PDB in batch mode)
shit_stat_insres = False
shit_stat_altpos = False
shit_stat_modres = False
shit_stat_misdns = False  
fastaseq = {}
pdbfile = ""
files_to_unlink = []

def download_pdb(pdb_id, dest_dir):
    """Download a PDB from RCSB."""
    url = f'http://www.rcsb.org/pdb/files/{pdb_id.upper()}.pdb.gz'
    dest = f'{abspath(dest_dir)}/{pdb_id}.pdb.gz'
    wget_cmd = f'wget --quiet {url} -O {dest}'
    print(wget_cmd)
    if remote_host:
        wget_cmd = f'ssh {remote_host} {wget_cmd}'
    popen(wget_cmd).readlines()
    return dest if exists(dest) else None

def check_and_print_pdb(count, residue_buffer, residue_letter):
    """Check backbone atoms and write residue lines."""
    global pdbfile, fastaseq
    hasCA = hasN = hasC = False
    for line in residue_buffer:
        atomname = line[12:16]
        occupancy = float(line[55:60])
        if atomname == " CA " and occupancy > 0.0: hasCA = True
        if atomname == " N  " and occupancy > 0.0: hasN = True
        if atomname == " C  " and occupancy > 0.0: hasC = True

    if hasCA and hasN and hasC:
        for line in residue_buffer:
            newnum = f'{count:4d} '
            pdbfile += line[:22] + newnum + line[27:]
        chain = residue_buffer[0][21]
        fastaseq[chain] = fastaseq.get(chain, '') + residue_letter
        return True
    return False

def get_pdb_filename(name):
    """Find a PDB file locally."""
    for ext in ['', '.pdb', '.pdb.gz', '.pdb1.gz']:
        if exists(name + ext):
            return name + ext
    return None

def open_pdb(name):
    """Open a PDB file (local or downloaded)."""
    filename = get_pdb_filename(name)
    if not filename:
        print(f"File for {name} doesn't exist, downloading...")
        filename = download_pdb(name[0:4].upper(), '.')
        if filename:
            files_to_unlink.append(filename)
        else:
            print(f"Failed to download {name}", file=sys.stderr)
            return None, None
    
    stem = basename(filename)
    for ext in ['.gz', '.pdb1', '.pdb']:
        stem = stem.replace(ext, '')
    
    if filename.endswith('.gz'):
        lines = popen(f'zcat {filename}', 'r').readlines()
    else:
        with open(filename, 'r') as f:
            lines = f.readlines()
    return lines, stem

def process_pdb(pdb_input, chainid, options, output_dir='.'):
    """Process a single PDB file."""
    global pdbfile, fastaseq, shit_stat_altpos, shit_stat_insres, shit_stat_modres, shit_stat_misdns
    
    # Reset global variables
    pdbfile = ""
    fastaseq = {}
    shit_stat_altpos = shit_stat_insres = shit_stat_modres = shit_stat_misdns = False
    
    lines, stem = open_pdb(pdb_input)
    if not lines:
        return False

    oldresnum = '   '
    count = 1
    residue_buffer = []
    residue_letter = ''

    for line in lines:
        if line.startswith('ENDMDL'): break
        if len(line) > 21 and (line[21] in chainid or options.allchains):
            if line[:4] != "ATOM" and line[:6] != 'HETATM': continue

            line_edit = line
            resn = line[17:20]

            # Handle modified residues
            if resn in modres:
                orig_resn = resn
                resn = modres[resn]
                line_edit = 'ATOM  ' + line[6:17] + resn + line[20:]

                if orig_resn == "MSE":
                    if line_edit[12:14] == 'SE':
                        line_edit = line_edit[:12] + ' S' + line_edit[14:]
                    if len(line_edit) > 75 and line_edit[76:78] == 'SE':
                        line_edit = line_edit[:76] + ' S' + line_edit[78:]
                else:
                    shit_stat_modres = True

            if resn not in longer_names:
                continue

            resnum = line_edit[22:27]

            if not resnum == oldresnum:
                if residue_buffer:
                    if not check_and_print_pdb(count, residue_buffer, residue_letter):
                        shit_stat_misdns = True
                    else:
                        count += 1
                residue_buffer = []
                residue_letter = longer_names[resn]

            oldresnum = resnum

            if line[26] != ' ':
                shit_stat_insres = True

            altpos = line[16]
            if altpos != ' ':
                shit_stat_altpos = True
                if altpos == 'A':
                    line_edit = line_edit[:16] + ' ' + line_edit[17:]
                else:
                    continue

            if options.removechain:
                line_edit = line_edit[:21] + ' ' + line_edit[22:]

            if options.keepzeroocc:
                line_edit = line_edit[:55] + " 1.00" + line_edit[60:]

            residue_buffer.append(line_edit)

    if residue_buffer:
        if not check_and_print_pdb(count, residue_buffer, residue_letter):
            shit_stat_misdns = True

    # Print summary
    flags = [
        "ALT" if shit_stat_altpos else "---",
        "INS" if shit_stat_insres else "---",
        "MOD" if shit_stat_modres else "---",
        "DNS" if shit_stat_misdns else "---",
        "OK" if fastaseq else "BAD"
    ]
    print(f"{stem} {chainid} {len(''.join(fastaseq.values())):5d} {' '.join(flags)}")

    # Save outputs
    if not options.nopdbout and pdbfile:
        outfile = join(output_dir, f"{stem}_{chainid}.pdb")
        with open(outfile, 'w') as f:
            f.write(pdbfile + "TER\n")

    if fastaseq:
        fasta_file = join(output_dir, f"{stem}_{chainid}.fasta")
        with open(fasta_file, 'w') as f:
            if options.allchains:
                f.write(f">{stem}_{chainid}\n{''.join(fastaseq.values())}\n")
            else:
                for chain in fastaseq:
                    f.write(f">{stem}_{chain}\n{fastaseq[chain]}\n")
    
    return True

def batch_process(input_dir, output_dir, chainid, options):
    """Process all PDB files in a directory."""
    if not exists(input_dir):
        print(f"Input directory not found: {input_dir}", file=sys.stderr)
        return False

    if not exists(output_dir):
        makedirs(output_dir)

    processed = 0
    for pdb_file in glob.glob(join(input_dir, "*.pdb*")):
        if pdb_file.endswith('.pdb') or pdb_file.endswith('.pdb.gz'):
            if process_pdb(pdb_file, chainid, options, output_dir):
                processed += 1

    print(f"\nProcessed {processed} PDB files. Output saved to: {output_dir}")
    return True

def main():
    parser = OptionParser(usage="%prog [options] <pdb> <chain id>", description=__doc__)
    parser.add_option("--nopdbout", action="store_true", help="Don't output a PDB.")
    parser.add_option("--allchains", action="store_true", help="Use all chains from the input PDB.")
    parser.add_option("--removechain", action="store_true", help="Remove chain information from output PDB.")
    parser.add_option("--keepzeroocc", action="store_true", help="Keep zero occupancy atoms in output.")
    parser.add_option("--batch", action="store_true", help="Process all PDBs in a directory.")
    parser.add_option("--input-dir", help="Directory containing PDB files (required for --batch).")
    parser.add_option("--output-dir", default="cleaned_pdbs", help="Output directory (default: ./cleaned_pdbs).")

    options, args = parser.parse_args()

    if options.batch:
        if not options.input_dir:
            parser.error("--batch requires --input-dir")
        chainid = args[0] if args else "A"
        batch_process(options.input_dir, options.output_dir, chainid, options)
    else:
        if len(args) != 2:
            parser.error("Must specify both the pdb and the chain id")
        process_pdb(args[0], args[1], options, options.output_dir)

    # Cleanup downloaded files
    for f in files_to_unlink:
        try:
            unlink(f)
        except:
            pass

if __name__ == "__main__":
    main()