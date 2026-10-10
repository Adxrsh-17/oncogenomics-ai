import pandas as pd
from pathlib import Path

MAF_FILE = Path(__file__).resolve().parent.parent / "results" / "tcga_luad_merged.maf"

def main():
    if not MAF_FILE.exists():
        print(f"MAF file not found: {MAF_FILE}")
        return
        
    df = pd.read_csv(MAF_FILE, sep="\t", comment="#", low_memory=False)
    
    total_variants = len(df)
    sift_col = "SIFT" if "SIFT" in df.columns else None
    polyphen_col = "PolyPhen" if "PolyPhen" in df.columns else None
    
    print(f"Total variants: {total_variants}")
    
    if sift_col:
        sift_valid = df[sift_col].notna().sum()
        print(f"Variants with SIFT: {sift_valid} ({sift_valid/total_variants:.1%})")
    else:
        print("SIFT column not found.")
        
    if polyphen_col:
        polyphen_valid = df[polyphen_col].notna().sum()
        print(f"Variants with PolyPhen: {polyphen_valid} ({polyphen_valid/total_variants:.1%})")
    else:
        print("PolyPhen column not found.")

if __name__ == "__main__":
    main()
