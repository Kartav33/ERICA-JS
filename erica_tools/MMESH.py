import pandas as pd

def readMMESH_fromFile(fullfilename):
    """
    Parameters
    ----------
    fullfilename : string
        The fully qualified (i.e., absolute) filename pointing to the MMESH 
        output .csv file

    Returns
    -------
    mmesh_results : pandas DataFrame
        The results of the MMESH run
    """
    
    #   Read the csv, parsing MultiIndex columns correctly
    mmesh_results = pd.read_csv(fullfilename, comment='#', header=[0,1], index_col=[0])
    
    #   Reset the index to be datetimes, rather than integers
    mmesh_results.index = pd.to_datetime(mmesh_results.index)
    
    return mmesh_results


def MMESH_slicing(df_MMESH, start, end, diag: str):
    diag_available = ['B_mag','B_mag_a', 'B_mag_loc', 'B_mag_scale', 
                      'n_tot','n_tot_a', 'n_tot_loc', 'n_tot_scale',
                      'P_dyn','P_dyn_a', 'P_dyn_loc', 'P_dyn_scale',
                      'u_mag','u_mag_a', 'u_mag_loc', 'u_mag_scale']
    
    if diag not in diag_available:
        return df_MMESH.loc[start:end]
        #raise ValueError(f"{diag} not found")
    
    i = diag_available.index(diag)
    return df_MMESH.loc[start:end, df_MMESH.columns[i]]