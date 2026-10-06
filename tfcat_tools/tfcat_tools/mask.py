import numpy
from shapely.geometry import MultiPoint, MultiPolygon, Polygon, Point
from shapely.ops import unary_union
from astropy.time import Time
import os

from tfcat.codec import load
from tfcat import TFCat

from spacepy import pycdf
from pathlib import Path
from typing import Union
from tqdm import tqdm

class Mask:
    def __init__(self,
                cat: TFCat,
                time: numpy.ndarray,
                freqs: numpy.ndarray,
                data: numpy.ndarray):
        
        
        sorted_indices = numpy.argsort(time)
        sorted_time = time[sorted_indices]
        sorted_data = data[sorted_indices]
        del time, data
        time, data = numpy.array(sorted_time), numpy.array(sorted_data)
        self.cat = cat
        self.time = time
        self.freqs = freqs
        print(self.freqs.shape)
        self.data = data
        self._compute_data_points()
        if self.cat.properties['instrument_host_name'] == 'Juno' and self.cat.properties['instrument_name'] == 'Waves' :
            self.mask = self._compute_mask(list(self.cat.data.fields['feature_type']['values'].keys()))
        else:
            self.mask = self._compute_mask(list(self.cat.data.fields['feature_type']['values'].keys()))
    @classmethod
    def from_file(cls,
                cat_file_name: Union[Path, str],
                *data_file_name: str):
        cat = TFCat.from_file(cat_file_name)
        time_list = []
        data_list = []
        check_list = []
        for file_name in data_file_name :
            if not os.path.basename(file_name) in check_list :
                data_list.append(pycdf.CDF(file_name)['Data'][:,:])
                time_list.append(pycdf.CDF(file_name)['Epoch'][:])
                check_list.append(os.path.basename(file_name))
        time = Time(numpy.concatenate(time_list), format='datetime').unix
        data = numpy.concatenate(data_list)
        freqs = pycdf.CDF(data_file_name[0])['Frequency'][:]
        return cls(cat, time, freqs, data)

    @property
    def time(self):
        return self._time
    
    @time.setter
    def time(self, time):
        if self.cat.properties['instrument_host_name'] == 'Juno' and self.cat.properties['instrument_name'] == 'Waves' :
            self._time = time[::15]
        else:
            self._time = time

    @property
    def data(self):
        return self._data
    
    @data.setter
    def data(self, data):
        if self.cat.properties['instrument_host_name'] == 'Juno' and self.cat.properties['instrument_name'] == 'Waves' :
            self._data = numpy.zeros(shape=(len(self.time), data.shape[-1]))
            for i in range(data.shape[-1]):
                self._data[:, i] = numpy.interp(numpy.arange(0, data.shape[0], 15), numpy.arange(0, data.shape[0], 1), data[:,i])
        else:
            self._data = data

    def _compute_mask(self, 
            feature_type: list):
        polygon_dict = self._build_polygons(feature_type)
        self._mask = {}
        for key in feature_type :
            _mask = numpy.zeros((self.data.size,))
            try:
                if not polygon_dict[key].is_empty :
                    select_points = polygon_dict[key].intersection(self._data_points)
                    _mask[[int(geom.z) for geom in select_points.geoms]] = 1
                    _mask = _mask.reshape(self.time.shape[0], self.freqs.shape[0])
            except AttributeError:
                _mask = _mask.reshape(self.time.shape[0], self.freqs.shape[0])
            self._mask[key] = (_mask > 0).reshape(self.time.shape[0], self.freqs.shape[0])
        return self._mask


    def _compute_data_points(self):
        times_coords = numpy.zeros_like(self.data)
        freqs_coords = numpy.zeros_like(self.data)
        times_coords[:,:], freqs_coords[:,:] = numpy.copy(self.time[:,None]), numpy.copy(self.freqs[None, :])
        coords = [(t, f) for t, f in zip(times_coords.flatten(), freqs_coords.flatten())]
        _data = self.data.flatten()
        index = numpy.arange(len(_data), dtype=int)
        self._data_points = MultiPoint([Point(x, y, z) for (x, y), z in zip(coords, index)]) 

    def _build_polygons(self,
                    feature_type: list):
        polygon_dict = {}
        for key in feature_type :
            features_polygon = []
            for i in range(len(self.cat)):
                label = self.cat.data.features[i]['properties']['feature_type']
                time_points = numpy.array(self.cat.data.features[i]['geometry']['coordinates'][0])[:,0]
                if any(time_points >= self.time.min()) and any(time_points <= self.time.max()):
                    if label in key:
                        features_polygon.append(Polygon(self.cat.feature(i).geometry.coordinates[0]))
            polygon_dict[key] = unary_union(MultiPolygon(features_polygon))
        return polygon_dict

    def catalog_bounds(self,
                        *feature_type: str,
                        invert: bool = False):
        
        for key in feature_type:
            if not "mask" in locals() :
                mask = numpy.copy(self.mask[key])
            else : 
                mask += self.mask[key]
        if invert == True :
            return ~mask
        else : 
            return mask

    def apply(self, 
            *feature_type: str,
            invert : bool = False):
        for key in feature_type :
            if not "mask" in locals() :
                mask = numpy.copy(self.mask[key])
            else : 
                mask += self.mask[key]
        if invert == True :
            return self.data*~mask
        else : 
            return self.data*mask




