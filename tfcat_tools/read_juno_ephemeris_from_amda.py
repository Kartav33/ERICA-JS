import numpy
import datetime

@numpy.vectorize
def amdadate_to_datetime(date):
    return datetime.datetime.strptime(date,"%Y-%m-%dT%H:%M:%S.%f")

def juno_ephemeris_from_amda(file):
    data_file = numpy.loadtxt(file,dtype="str")
    date = amdadate_to_datetime(data_file[:,0][:])
    x_coord = data_file[:,1].astype(float)
    y_coord = data_file[:,2].astype(float)
    z_coord = data_file[:,3].astype(float)

    return (date,x_coord,y_coord,z_coord)