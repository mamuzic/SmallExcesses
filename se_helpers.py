import pandas as pd
import os, sys
from pathlib import Path


def load_scan(dirname, fname):
    os.chdir(dirname)

    # Read the CSV
    # The file seems to have an extra comma at the start of the header,
    # so we use index_col=0 to ignore the first unnamed column.
    df = pd.read_csv(fname, index_col=0)

    print("✅ Data loaded successfully!")
    if False:
        print(df.head())

    if False:
        # Optional: check column names
        print("\nColumns:")
        print(df.columns.tolist())

    return df

def extract_se(dirname,fname):
    os.chdir(dirname)
    df = load_scan(dirname, fname)

    n_models, n_info = df.shape
    print("Total number of models:", n_models)
    print("Number of info columns:", n_info)

    columns = df.columns.tolist()
    if False: print("Columns:", columns)

    analyses = []
    for col in columns:
        if col.endswith('__ExpCLs'):
            name = col.split('__')[0]
            if f"{name}__ObsCLs" in columns:
                analyses.append(name)

    print("Analyses with ExpCLs or ObsCLs columns:")
    print(analyses)

    results = {}

    for ana in analyses:
        exp_col = f"{ana}__ExpCLs"
        obs_col = f"{ana}__ObsCLs"

        if exp_col in df.columns and obs_col in df.columns:
            # Find models satisfying ExpCLs >= 0.05 and ObsCLs < 0.05
            mask = (df[exp_col] <= 0.05) & (df[obs_col] > 0.05)
            #subset = df.loc[mask, ["Model_number", exp_col, obs_col]]
            subset = df.loc[mask]
        
            results[ana] = subset
            print(f"\n🔹 {ana}: found {len(subset)} models with ExpCLs <= 0.05 and ObsCLs > 0.05")
        else:
            print(f"\n⚠️ Skipping {ana} (missing columns)")

    if False: print(results["2L0J"].head())

    excess_all = pd.concat(results.values(), ignore_index=True)
    excess_all = excess_all.drop_duplicates(subset=["Model_number"])
    if False: print(excess_all.head())

    n_models, n_info = excess_all.shape
    print("Number of SE models", n_models)

    base, ext = os.path.splitext(fname)
    se_fname = dirname+"/"+base + "_SE" + ext
    excess_all.to_csv(se_fname, index=False)
    if os.path.exists(se_fname):
        print("✅ Created Small Excesses skim csv file", se_fname)
    else:
        print("❌ Error producing the csv file.", se_fname)

    return excess_all

def copy_model_slha(slha_dir, modelid):
    # Make new folder name with _SE suffix
    slha_dir = Path(slha_dir)
    model = str(modelid) + ".slha"
    se_slha_dir = slha_dir.parent / (slha_dir.name + "_SE")
    se_slha_dir.mkdir(parents=True, exist_ok=True)

    for f in slha_dir.rglob(model):
        print("✅ Found:", f)

        if False:
            with open(f, "r") as f:
                print(f.read())

        # Destination file path (same filename inside the new folder)
        dst = se_slha_dir / f.name

        # Copy file contents (no shutil)
        dst.write_bytes(f.read_bytes())

        print(f"✅ Copied {f.name} to {dst}")

def copy_all_model_slha(slha_dir, excess_all):
    # Loop on all models
    for modelid in excess_all.Model_number:
        print(modelid)
        copy_model_slha(slha_dir, modelid)

def main():
    dirname = "/Users/mamuzic/MyWork/LJUBLJANA/WP2_pMSSMML_SMASH/SmallExcesses"
    bino_csv = "../EWKpMSSM_HepDATA/Bino-DM.csv"
    bino_slha_dir = "../EWKpMSSM_HepDATA/Bino-DM"

    # load_scan(dirname, bino_csv) # Test
    excess_all = extract_se(dirname, bino_csv)
    copy_all_model_slha(dirname+"/"+bino_slha_dir, excess_all)


# This ensures code runs only when executed directly,
# not when imported as a module.
if __name__ == "__main__":
    main()