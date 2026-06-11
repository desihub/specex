import os
import sys
import re
import subprocess

# Ensure we use the current workspace code
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))

from specex.io import read_lamp_lines

def demonstrate_quirk():
    lamp_file = 'py/specex/data/specex_linelist_desi.txt'
    
    # 1. Simple count of lines not starting with '#'
    with open(lamp_file, 'r') as f:
        clean_lines = [l for l in f if l.strip() and not l.strip().startswith('#')]
    
    # 2. Use our ported parser (which now matches C++ behavior)
    ported_lines = read_lamp_lines(lamp_file)
    
    print("--- Specex Lamp Line Parsing Comparison ---")
    print(f"File: {lamp_file}")
    print(f"Lines NOT starting with '#':      {len(clean_lines)}")
    print(f"Lines loaded by Specex (Ported):  {len(ported_lines)}")
    
    diff = len(ported_lines) - len(clean_lines)
    print(f"\nQuirk: Specex loads {diff} extra lines that appear commented out.")
    
    print("\nExample of 'commented' lines Specex actually uses:")
    # These are lines where the name starts with # but the rest of the line is valid data
    for l in ported_lines:
        if l['name'].startswith('#'):
            print(f"  Ion: {l['name']:<10} Wave: {l['wave']:<10} Score: {l['score']}")

    print("\nReason: The C++ code uses 'is >> ion >> wave >> score'.")
    print("If 'ion' is '#ArI', it is simply treated as a name, not a comment.")

if __name__ == "__main__":
    demonstrate_quirk()
