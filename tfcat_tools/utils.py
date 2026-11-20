import numpy
from datetime import *
from astropy.time import Time
from webgeocalc import Calculation

def nozeros(x):
    return numpy.where(x != 0, x, 1)

def xyz_to_rtp(xyz):
    rtp = numpy.zeros(shape=xyz.shape)
    rtp[..., 0] = numpy.sqrt(xyz[..., 0]**2+xyz[..., 1]**2+xyz[..., 2]**2)
    rtp[..., 1] = numpy.arcsin(xyz[..., 2]/nozeros(rtp[..., 0]))*180/numpy.pi
    rtp[..., 2] = numpy.arctan2(xyz[..., 1], xyz[..., 0])*180/numpy.pi
    return rtp

def set_ephemeris(time, ephemeris, local_time_ephemeris):
    tmin, tmax = time.min(), time.max()
    new_ephemeris_time = ephemeris[0][(ephemeris[0] >= tmin) & (ephemeris[0] <= tmax)]
    new_local_time_ephemeris_time = local_time_ephemeris[0][(local_time_ephemeris[0] >= tmin) & (local_time_ephemeris[0] <= tmax)]
    _ephemeris_xyz = numpy.zeros(shape=(time.shape[0],3))
    _local_time_ephemeris_xyz = numpy.zeros(shape=(time.shape[0],3))
    for i in range(3) :
        _ephemeris = ephemeris[i+1][(ephemeris[0] >= tmin) & (ephemeris[0] <= tmax)]
        _ephemeris_xyz[:, i] = numpy.interp(time, new_ephemeris_time, _ephemeris)
        
        _ephemeris = local_time_ephemeris[i+1][(local_time_ephemeris[0] >= tmin) & (local_time_ephemeris[0] <= tmax)]
        _local_time_ephemeris_xyz[:, i]  = numpy.interp(time, new_local_time_ephemeris_time, _ephemeris)
        
    _ephemeris = numpy.zeros(shape=(time.shape[0],4))
    _ephemeris[:, 0:3] = xyz_to_rtp(_ephemeris_xyz)
    del _ephemeris_xyz
    _local_time_ephemeris = xyz_to_rtp(_local_time_ephemeris_xyz)
    _ephemeris[:, 3] = _local_time_ephemeris[:, 2]
    return _ephemeris

@numpy.vectorize
def webgeocalc_date_to_datetime(date_YMD,date_hours):
    return datetime.datetime.strptime(date_YMD+"T"+date_hours,"%Y-%m-%dT%H:%M:%S.%f")

@numpy.vectorize
def unix_date_to_datetime(date_unix):
    return datetime.datetime.fromtimestamp(date_unix)


def juno_ephemeris_from_webgeocalc(file, planetary_radius=False):
    if planetary_radius == False:
        print("### User needs to give a planetary radius to have ephemeris in planetary units ###")
    data_file = numpy.loadtxt(file,dtype="str")
    date = webgeocalc_date_to_datetime(data_file[:,0],data_file[:,1])
    x_coord = data_file[:,5].astype(float)
    y_coord = data_file[:,6].astype(float)
    z_coord = data_file[:,7].astype(float)

    return (date,x_coord[:]/planetary_radius,y_coord[:]/planetary_radius,z_coord[:]/planetary_radius)

#thresholding data
def thresholding_data(data,background):

    occurrence = numpy.zeros(data.shape)
    for ifreq in range(len(background)):
        mask = data[:, ifreq] > background[ifreq]
        occurrence[mask, ifreq] = 1
    return occurrence

@numpy.vectorize
def from_UTC_to_astropy_Time(time_UTC):
    return Time(time_UTC[:-4]).unix

def get_ephem_from_wgc(observer, time, api = None):
    """
    Get ephemeris data from WGC servers through API
    :param observer:
    :type observer: dict
    :param time:
    :type time: dict
    :return:
    :rtype: dict
    """

    wgc_kernels = [
                    "pds/wgc/mk/ground_stations_v0018.tm",
                    "pds/wgc/mk/solar_system_v0056.tm",
                    "pds/wgc/mk/latest_lsk_v0004.tm",
                    ]

    time_interval = [time['time_start'].strftime("%Y-%m-%d %H:%M:%S"), time['time_end'].strftime("%Y-%m-%d %H:%M:%S")]                    

    if api == None:
        api = ''
    elif api.lower() == 'obspm':
        api = 'https://voparis-webgeocalc2.obspm.fr/geocalc/api'
        #api = "https://wgc.jpl.nasa.gov:8443/webgeocalc/api"
        #api = 'https://wgc.obspm.fr/webgeocalc/api'
    calc = Calculation(
        api = api,
        kernel_paths =  wgc_kernels,
        intervals = time_interval,
        time_step = time['time_step'],
        time_step_units = time['time_units'],
        calculation_type = 'STATE_VECTOR',
        target = observer['NAME'],
        observer = observer['PARENT'],
        reference_frame = observer['referenceFrame'],
        aberration_correction = 'NONE',
        state_representation = 'PLANETOGRAPHIC',
    )

    calc.payload
    calc.run()

    longitude = numpy.array(calc.results['LONGITUDE'])  # must be in LongW for ExPRES
    latitude = numpy.array(calc.results['LATITUDE'])
    distance = numpy.array(calc.results['ALTITUDE'])
    time_ephem = numpy.array(calc.results['DATE'])

    result_data = {
                    'time': time_ephem,
                    'longitude': -longitude, # So that the longitude is counting West (SIII for Jupiter)
                    'latitude': latitude,
                    'distance': distance,

                }
    return result_data
    