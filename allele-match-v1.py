import pandas as pd
# Currently use, last use 27 Nov 2025
# Updated 10/2025
# =================================================================================
# FUNCTION 1: Prepare the initial cDNA set 
# =================================================================================
def prepare_cDNA_set(sample_file_path, sheet_name):
    """
    Reads and processes a sample Excel file to generate a set of all cDNA changes
    and a separate set for homozygous changes.
    """
    # --- Configuration for the sample file ---
    COLUMN_A_GENE = 'Gene'
    COLUMN_N_REGION = 'Region'
    COLUMN_K_NUCLEOTIDE = 'Nucleotide'
    COLUMN_TO_JOIN = 'Nucleotide'
    COLUMN_M_ZYGOSITY = 'Zygosity'  
    GENE_TO_FIND = 'ABO'
    REGION_TO_FIND = 'Exonic'
    
    print(f"--- [Part 1] Reading data from '{sample_file_path}' ---")
    try:
        df = pd.read_excel(sample_file_path, sheet_name=sheet_name)
    except FileNotFoundError:
        print(f"--- ERROR: The file '{sample_file_path}' was not found. ---")
        return None, None 

    # sort ABO GENE and EXON
    print(f"--- [Part 1] Sorting and Filtering for Gene='{GENE_TO_FIND}' and Region='{REGION_TO_FIND}' ---")
    df_sorted = df.sort_values(by=[COLUMN_A_GENE, COLUMN_N_REGION, COLUMN_K_NUCLEOTIDE])
    filtered_df = df_sorted[
        (df_sorted[COLUMN_A_GENE] == GENE_TO_FIND) 
        
        # & # ==== uncomment to select only exonic region - will lose some allele !!!!
        # (df_sorted[COLUMN_N_REGION] == REGION_TO_FIND) # ==== uncomment to select only exonic region
        
    ]

    if filtered_df.empty:
        print("--- [Part 1] No matching rows found in sample file. ---")
        print("--- [Converting to homozygous c.261del] ---")
        my_cDNA_set_ini = {"c.261del"}
        homo_set = {"c.261del"}

    # Create the set of ALL variants
    else:
        initial_join = "','".join(filtered_df[COLUMN_TO_JOIN].dropna().astype(str))
        my_cDNA_set_ini = set(initial_join.split("','"))
    
    # --- NEW: Create the set of HOMOZYGOUS variants ---
        homo_df = filtered_df[filtered_df[COLUMN_M_ZYGOSITY] == 'hom']
        homo_set = set(homo_df[COLUMN_TO_JOIN].dropna().astype(str))

    # --- END OF NEW SECTION ---
    print('--- [Part 1] converting reference allele from O01.01 (Genomic) to A1.01 (ISBT) ---')
    
    print('\n--- Genomic cDNA set\n---', my_cDNA_set_ini)
    if homo_set:
        print('--- Homozygous set = ', homo_set)
    else:
        print('--- Homozygous set = none')
    print('\n--- [ c.260_262= indicates non c.261del allele. ]')
    print('--- [ Allele without c.260_262= will be converted to c.261del ]')

    if ('c.260_262=' in my_cDNA_set_ini):
        if 'c.260_262=' in homo_set:
            my_cDNA_set_ini.remove('c.260_262=')
            homo_set.remove('c.260_262=')
        else:
            my_cDNA_set_ini.remove('c.260_262=')
            my_cDNA_set_ini.add('c.261del')
    else:
        my_cDNA_set_ini.add('c.261del')
        homo_set.add('c.261del')
        # pass

    print('\n--- ISBT cDNA set\n---', my_cDNA_set_ini)
    if homo_set:
        print('--- Homozygous set = ', homo_set)
    else:
        print('--- Homozygous set = none')

    print(f"\n--- [Part 1] Successfully created initial cDNA set with {len(my_cDNA_set_ini)} items. ---")
    print(f"--- [Part 1] Found {len(homo_set)} homozygous variants. ---")
    
    return my_cDNA_set_ini, homo_set

# =================================================================================
# FUNCTION 2: Find the matching allele pair from the allele database
# =================================================================================
def find_allele_pair(allele_db_path, sheet_name, my_cDNA_set, homo_set_final, sample_name):
    """
    Reads an allele database, and finds a pair of alleles that match the provided cDNA set.
    """
    # --- Configuration for the allele database file ---
    COLUMN_D_NAME = 'Allele'
    COLUMN_L_NAME = 'DNA Change'

    print(f"--- [Part 2] Reading data from '{allele_db_path}' ---")
    try:
        df = pd.read_excel(allele_db_path, sheet_name=sheet_name)
    except FileNotFoundError:
        print(f"--- ERROR: The file '{allele_db_path}' was not found. ---")
        return

    print("--- [Part 2] Creating dictionary from allele database ---")
    df[COLUMN_L_NAME] = df[COLUMN_L_NAME].fillna('').str.split('\n')
    df_exploded = df.explode(COLUMN_L_NAME)
    df_exploded = df_exploded[df_exploded[COLUMN_L_NAME].str.strip() != '']
    data_dict = df_exploded.groupby(COLUMN_D_NAME)[COLUMN_L_NAME].apply(set).to_dict()

    print(f"\n--- [Part 2] Comparing your set for sample: {sample_name} ---")
    print(f"--- Total positions to match: {len(my_cDNA_set)} ---\n")
    print('cDNA set = ', my_cDNA_set)
    # ========
    # ========
    # Identify variants represented in the allele database
    database_variants = set().union(*data_dict.values())

    # Homozygous variants present in the allele database
    homo_in_database = homo_set_final & database_variants

    # Homozygous variants absent from the allele database
    homo_not_in_database = homo_set_final - database_variants

    print("Homozygous variants in database:", homo_in_database)
    print("Homozygous variants not in database:", homo_not_in_database)
    # ========
    # ========
    # if homo_set_final:
    #     print('Homozygous set = ', homo_set_final)
    # else:
    #     print('Homozygous set = none')
    # ========
    # ========
    # Find all alleles that are a subset of the main set
    found_match = False
    exact_match = False

    for key, value_set in data_dict.items():

        # Check if the set from the Excel file is exactly equal to your set
        if my_cDNA_set == value_set:
            found_match = True
            exact_match = True
            print(f"\n\n✅✅✅ \nExact Match Found for allele: '{key}'\n")
            print(f"The DNA change was: {value_set} \n")
            print(f"The your change was: {value_set} \n\n")
            key_exact = key
            value_exact = value_set
            
    for key, value_set in data_dict.items():
        # Find the intersection between your set and the set from the Excel file
        intersection = my_cDNA_set.intersection(value_set)
        
        # If there are any common items, we have a match
        if intersection:
            found_match = True
            print(f"---------\n✅'Possible allele: {key} \n {value_set}'\n position count = {len(value_set)} ")
            print(f"Remain = {len(my_cDNA_set - intersection)}:{list(my_cDNA_set - intersection)}\n")
    
    if not found_match:
        print("--- No matching cDNA positions were found in the Excel file. ---")

## added 28 AUG ##
    # print('-- Recording alleles name --')
    keys = list(data_dict.keys())
    allele_issubset = {}

    # print('-- Counting alleles name --')
    count_allele = len(data_dict)

    # print('-- Finding all possible alleles --')
    # print(f'Your cDNA set: {my_cDNA_set}')

    for i in range(count_allele):
        if data_dict[keys[i]].issubset(my_cDNA_set): 
            allele_issubset[f'{keys[i]}']=data_dict[keys[i]]
            # print(f'{keys[i]}: {data_dict[keys[i]]}')  # ====== This will print all possible alleles ======= #
    

    subset_keys = list(allele_issubset.keys())
    count_issubset = len(subset_keys)

    found_pair = False
    print('\n======== \nFinding all possible pairs \n========')
    for i in range(count_issubset):
        for j in range(i + 1, count_issubset):
            allele1_set = allele_issubset[subset_keys[i]]
            allele2_set = allele_issubset[subset_keys[j]]

            intersect = allele1_set.intersection(allele2_set)
            union = allele1_set.union(allele2_set)

            if intersect == homo_in_database and union.issubset(my_cDNA_set) and union != my_cDNA_set:
                found_pair = True

                other_cDNA = my_cDNA_set - union
                original_other_cDNA = other_cDNA.copy()

                other_cDNA = {
                    i for i in other_cDNA
                    if any(i in variants for variants in data_dict.values())
                }

                discarded_cDNA = original_other_cDNA - other_cDNA

                if other_cDNA == set():
                    print('\n✅✅✅ Found Perfect Allele Pair!  ✅✅✅')
                    print(f'Allele 1: {subset_keys[i]} -> {allele1_set}')
                    print(f'Allele 2: {subset_keys[j]} -> {allele2_set}')
                    print(f'\nTotal cDNA: {union}')
                    print(f"\nRemaining cDNA changes in allele table: {other_cDNA}")
                    print(f"\nDiscarded cDNA changes not found in allele table: {discarded_cDNA}")

                    if intersect is True: 
                        print(f'Homozygous: {intersect}\n')
                    else:
                        print('Homozygous: none\n')

                if other_cDNA != set():
                    print('\n✅ Partial Matched Allele Pair! ✅')
                    print(f'Allele 1: {subset_keys[i]} -> {allele1_set}')
                    print(f'Allele 2: {subset_keys[j]} -> {allele2_set}')
                    print(f'\nTotal cDNA: {union}')
                    print(f"\nRemaining cDNA changes in allele table: {other_cDNA}")
                    print(f"\nDiscarded cDNA changes not found in allele table: {discarded_cDNA}")
                    if intersect is True: 
                        print(f'Homozygous: {intersect}\n')
                    else:
                        print('Homozygous: none\n')
               
            if intersect == homo_in_database and union == my_cDNA_set:
                print('\n✅✅✅ Found Perfect Allele Pair!  ✅✅✅')
                print(f'Allele 1: {subset_keys[i]} -> {allele1_set}')
                print(f'Allele 2: {subset_keys[j]} -> {allele2_set}')
                print(f'\nTotal cDNA: {union}')

    if exact_match and homo_set_final == my_cDNA_set:
        print("\n✅✅✅ Found Perfect Allele Pair! ✅✅✅")
        print(f"Allele 1: {key_exact}\ncDNA change: {value_exact}")
        print(f"Allele 2: {key_exact}\ncDNA change: {value_exact}")    
    elif exact_match and not found_pair:
        print("\n✅✅✅ Found Perfect Allele Pair! ✅✅✅")
        print("Allele 1: ABO*A1.01 \ncDNA change: none")
        print(f"Allele 2: {key_exact}\ncDNA change: {value_exact}")
    elif not exact_match and not found_pair: 
        print('--- Cannot find any allele pair that meets the criteria. ---')


# =================================================================================
# MAIN EXECUTION BLOCK
# =================================================================================
if __name__ == "__main__":

    # --- Database Configuration ---
    ALLELE_DATABASE_FILE = 'ABO_allele_table_v8.xlsx'

    print(f'\n\n****\nISBT allele table: {ALLELE_DATABASE_FILE}')
    print('PLEASE CHECK CURRENT VERSION BEFORE PROCEED\n****')

    # --- Overall Configuration ---
    SAMPLE_FILE = input('\nENTER ANNOTATED SAMPLE NAME (e.g., sample.xlsx):')
    SAMPLE_SHEET = 'sheet1'
    ALLELE_SHEET = 'Sheet1'
    VARIANTS_TO_REMOVE = {'c.260_262='}
    # This line is no longer needed as the homo set is generated from the file
    # HOMO_SET_CRITERIA = set()

    ## ============= HEADER OF PROGRAM ================================= ##
    print(f'\n=== ANALYZING SAMPLE: {SAMPLE_FILE} ===') 
    print(f'=== ALLELE_DATABASE: {ALLELE_DATABASE_FILE} === \n')

    # --- Step 1: Run the first part to get the initial and homo sets ---
    # CHANGED: Unpack the returned tuple into TWO variables
    initial_cDNA_set, homo_set_from_sample = prepare_cDNA_set(
        sample_file_path=SAMPLE_FILE, 
        sheet_name=SAMPLE_SHEET
    )
    
    # Proceed only if the first part was successful
    if initial_cDNA_set is not None:
        final_cDNA_set = initial_cDNA_set
        # --- Step 2: Modify the set as needed ---

        # --- Step 3: Run the second part using the final set ---
        find_allele_pair(
            allele_db_path=ALLELE_DATABASE_FILE,
            sheet_name=ALLELE_SHEET,
            my_cDNA_set=final_cDNA_set,
            homo_set_final=homo_set_from_sample,
            sample_name=SAMPLE_FILE
        )
       
    else:
        print("\n--- Script stopped because the initial cDNA set could not be created. ---")