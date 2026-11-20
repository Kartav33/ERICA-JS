import os
from pathlib import Path

# Research workspace on the server
WORKSPACE = Path(os.environ.get("ERICA_WORKSPACE", "/obs/kkesri"))

# Data hierarchy
DATA_ROOT = WORKSPACE / "data_and_outputs"
INPUT_DATA = DATA_ROOT / "spacecraft_and_model"
OUTPUTS = DATA_ROOT / "datasets_ERICA-JS"
EPHEM = DATA_ROOT / "ephem"

MMESH = INPUT_DATA / "MMESH_data"
TFCAT = INPUT_DATA / "TFCat_catalog_2.1"

# Folder organization

'''
/obs/kkesri
    -- ERICA-JS
        -- tfcat_tools
        -- erica_tools
        -- 01.ipynb, 02.ipynb ....
        -- paths.py
    
    -- data_and_outputs
        -- spacecraft_and_model
            -- MMESH
            -- TFCat_catalog
        -- outputs_ERICA-JS
        -- plots
        -- ephem
    
    -- pyKOM
'''