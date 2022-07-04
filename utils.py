import pandas as pd

def read_df_file(f_path):
    if f_path.endswith(".xlsx"):
        data = pd.read_excel(f_path)
    elif f_path.endswith(".csv"):
        data = pd.read_csv(f_path)
    else:
        assert False, "ERROR Parsing GENEpop format file"
    return data