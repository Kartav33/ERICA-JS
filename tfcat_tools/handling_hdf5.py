import numpy
from h5py import File
from tfcat_tools import *
from tfcat import TFCat
from spacepy import pycdf

def data_to_hdf5(file_name: str,
                time: numpy.ndarray,
                freqs: numpy.ndarray,
                data: numpy.ndarray,):
        """
        Saves the data to disk as a pre-processed HDF5 file.
        """
        output_file = File(file_name, 'w')
        output_file.attrs.create('observer', 'Juno')

        # Has to be done differently as this is an Astropy quantity
        output_file.create_dataset('Time', data=time)
        output_file['Time'].attrs.create('units', 'unix')

        output_file.create_dataset('Frequency', data=freqs)
        output_file['Frequency'].attrs.create('units',  'kHz')

        output_file.create_dataset('Data', data=data)
        output_file['Data'].attrs.create('units',  'W/m2/Hz')


        output_file.close()

def load_hdf5_data(file_name: str):
    """
    Reads a datafile in HDF5 format using the config file set earlier.
    """
    file: File = File(file_name)
    time = file['Time']
    freqs = file['Frequency']
    flux = file['Flux']
    return time, freqs, flux


# Function to read and concatenate variables from HDF5 files
def concatenate_hdf5_files(file_list, variable_names, path_files = './'):
    concatenated_data = {var: [] for var in variable_names}

    for file in file_list:
        with File(path_files+file, 'r') as f:
            for var in variable_names:
                concatenated_data[var].append(f[var][:])

    # Concatenate lists into numpy arrays
    concatenated_data = {var:numpy.concatenate(concatenated_data[var], axis=0) for var in variable_names}

    return concatenated_data

# Function to recursively flatten the nested structure
def flatten(nested_list):
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten(item))
        else:
            flat_list.append(item)
    return flat_list

def sort_cat(cat: TFCat):
    crs = cat.crs
    properties = cat.properties
    fields= cat.data.fields
    schema = cat.data["$schema"]
    event = FeatureCollection(cat.data, properties=properties, fields=fields, crs=crs, schema=schema).as_table()
    event.sort(['tmin'])
    return TFCat(FeatureCollection().from_table(table=event, properties=properties, fields=fields, crs=crs, schema=schema))

def load_intensity_from_hdf5_file(file_name: str):
    """
    Reads a datafile in HDF5 format using the config file set earlier.
    """
    file: File = File(file_name)
    normalized_intensity_bKOM = numpy.array(file['normalized_intensity_bKOM'])

    return normalized_intensity_bKOM  


def load_background_data(cdf_file_path):
    with pycdf.CDF(cdf_file_path) as cdf_file:
        background = numpy.array(cdf_file["Background"])
        sigma = numpy.array(cdf_file["Sigma"])
        frequency = numpy.array(cdf_file["Frequency"])
        
        #for attr in cdf_file["Background"].attrs:        
        #    print(f"{attr}: {cdf_file['Background'].attrs[attr]}")
    
    return background, sigma, frequency
