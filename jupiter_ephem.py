from astropy.coordinates import Longitude
from astropy.time import Time
from astropy.units import Unit
import os
from webgeocalc import StateVector, GFCoordinateSearch

__ALL__ = ['get_cml', 'get_meridian_transit', 'get_jovian_moon_phase']

#os.environ["HTTP_PROXY"] = "http://localhost:3128"
os.environ["WGC_URL"] = 'https://voparis-webgeocalc2.obspm.fr/geocalc/api'

_CORR = {
    True: "LT",
    False: "XLT",
    None: "NONE"
}


def _get_time_longitude(epochs: Time, observer: str, target: str, to_observer=True) -> (Time, Longitude):
    """
    Compute the time_at_target and the longitude using WGC StateVector API.
    :param epochs: Time array at observer's location
    :param observer: Observer's name (EARTH, IO, GANYMEDE, EUROPA, CALLISTO, AMALTHEA)
    :param target: Target's name (JUPITER)
    :param to_observer: True when correcting for light time to observer. Set to False to correct for the reverse path,
    and to None for no correction.
    :return: Time at target and Longitude of Sub-Observer's point (in IAU_JUPITER system III longitude)
    """

    kernels = 1
    if observer.lower() == 'juno':
        kernels = 'JUNO archive'
    calc = StateVector(
        api = 'https://voparis-webgeocalc2.obspm.fr/geocalc/api',
        kernels=kernels,
        times=list(epochs.isot),
        target=target.upper(),
        observer=observer.upper(),
        reference_frame='IAU_JUPITER',
        aberration_correction=_CORR[to_observer],
        state_representation="LATITUDINAL"
    )
    print(calc.api)
    result = calc.run()
    # NB: Longitude is in IAU_JUPITER but centered on observer's location.
    # Hence, Longitude of sub-observer point must be flipped by 180°
    return Time([item.replace('UTC','').strip() for item in result['TIME_AT_TARGET']]), \
           Longitude([item + 180 for item in result['LONGITUDE']], unit=Unit('deg'))


def _get_meridian_time(start_epochs: Time, end_epochs: Time, time_step: int, time_step_units: str, observer: str, target: str, earth_longitude: Longitude) -> Time:
    """
    Compute the Meridian Transit times of target for an observer at a given Longitude

    :param start_epochs: start times of observation
    :param end_epochs: end times of observation
    :param target: Observation target
    :param earth_longitude: Observer's Longitude (Earth coordinates)
    :return: Meridian transit epochs of target at observer's location.
    """

    # We need arrays of epochs, so we fix the Scalar input case
    if start_epochs.ndim < 1:
        start_epochs = Time([start_epochs])
    if end_epochs.ndim < 1:
        end_epochs = Time([end_epochs])

    kernels = 1
    if observer.lower() == 'juno':
        kernels = 'JUNO archive'
    calc = GFCoordinateSearch(
        api='http://voparis-webgeocalc2.obspm.fr:8080/geocalc/api',
        kernels=kernels,
        intervals=[{'startTime': start, 'endTime': end} for start, end in zip(start_epochs.isot, end_epochs.isot)],
        time_step=time_step,
        time_step_units=time_step_units,
        target=target,
        observer=observer,
        reference_frame="EARTH_FIXED",
        aberration_correction="LT",
        output_duration_units="SECONDS",
        condition={
            "coordinateSystem": "PLANETOGRAPHIC",
            "coordinate": "LONGITUDE",
            "relationalCondition": "=",
            "referenceValue": earth_longitude.value
        }
    )
    print(calc.api)
    return calc.run()


def _cml(longitude: Longitude) -> list:
    """
    Reverses the sign, since CML is a west-ward longitude system.
    NB: CML = Central Meridian Longitude = Longitude of the Sub-Observer point
    """
    return (360 - longitude.value) % 360


def get_jovian_moon_phase(epochs: Time, moon_name: str, observer: str) -> list:
    """
    Compute the Moon phase.

    :param epochs: Observation Epoch (at observer, i.e., Earth)
    :param moon_name: Name of the selected Jovian moon
    :return: Phase in degrees
    """
    tjup, ljup = _get_time_longitude(epochs, observer, 'JUPITER', True)
    _, lmoon = _get_time_longitude(tjup, moon_name, 'JUPITER', True)
    # NB: the phase reference (phase = 0°) is the opposition
    return (lmoon - ljup).value + 180


def get_cml(epochs: Time, observer: str) -> list:
    """
    Compute the CML of Earth for a list of Epochs

    :param epochs: list of times
    :return: list of CML
    """
    _, ljup = _get_time_longitude(epochs, observer, 'JUPITER')
    return _cml(ljup)


def get_jupiter_meridian_transit(start_epochs: Time, end_epochs: Time) -> list:

    target = "JUPITER"
    nancay_longitude = Longitude(2.19503, unit="deg")
    result = _get_meridian_time(start_epochs, end_epochs, target, nancay_longitude)

    return [" ".join(item.split(' ')[:-1]) for item in result['DATE']] if len(result['DATE']) > 0 else None
