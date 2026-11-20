from astropy.table import vstack
import numpy
from tfcat import TFCat 
from tfcat_tools.featureaddons import FeatureCollection

def merge(*cats: TFCat):
    crs = cats[0].crs
    properties = cats[0].properties
    fields= cats[0].data.fields
    schema = cats[0].data["$schema"]
    for cat in cats :
        if not crs == cat.crs or not properties == cat.properties or not fields == cat.data.fields:
           raise TypeError("catalogs to merged should have the same crs, properties and fields")
        event = FeatureCollection(cat.data, properties=properties, fields=fields, crs=crs, schema=schema).as_table()
        if not "merged_event" in locals() :
            merged_event = event
        else :
            merged_event = vstack([merged_event, event])
        
    return TFCat(FeatureCollection().from_table(table=merged_event, properties=properties, fields=fields, crs=crs, schema=schema))

def remove_feature(cat: TFCat, 
            time: float,
            frequency: float, 
            feature_type: str):
    crs = cat.crs
    properties = cat.properties
    fields= cat.data.fields
    schema = cat.data["$schema"]
    event = FeatureCollection(cat.data, properties=properties, fields=fields, crs=crs, schema=schema).as_table()
    removed_event = event[numpy.where(~((event['tmin'] <= time) & (event['tmax'] >= time)
                          & (event['fmin'] <= frequency) & (event['fmax'] >= frequency)
                          & (event['feature_type'] == feature_type)))]
    return TFCat(FeatureCollection().from_table(table=removed_event, properties=properties, fields=fields, crs=crs, schema=schema))

def keep_features_from_year(cat: TFCat,
                            time: [float, float]):
    crs = cat.crs
    properties = cat.properties
    fields= cat.data.fields
    schema = cat.data["$schema"]
    event = FeatureCollection(cat.data, properties=properties, fields=fields, crs=crs, schema=schema).as_table()
    keeped_event = event[numpy.where((event['tmin'] >= time[0]) & (event['tmax'] < time[-1]))]
    return TFCat(FeatureCollection().from_table(table=keeped_event, properties=properties, fields=fields, crs=crs, schema=schema))

def keep_features_from_frequency(cat: TFCat,
                                 frequency: [float, float]):
    crs = cat.crs
    properties = cat.properties
    fields= cat.data.fields
    schema = cat.data["$schema"]
    event = FeatureCollection(cat.data, properties=properties, fields=fields, crs=crs, schema=schema).as_table()
    keeped_event = event[numpy.where((event['fmin'] >= frequency[0]) & (event['fmax'] < frequency[-1]))]
    return TFCat(FeatureCollection().from_table(table=keeped_event, properties=properties, fields=fields, crs=crs, schema=schema))


def keep_features_from_type(cat: TFCat,
                            feature_type: str):
    crs = cat.crs
    properties = cat.properties
    fields= cat.data.fields
    schema = cat.data["$schema"]
    event = FeatureCollection(cat.data, properties=properties, fields=fields, crs=crs, schema=schema).as_table()
    keeped_event = event[numpy.where((event['feature_type'] == feature_type))]
    return TFCat(FeatureCollection().from_table(table=keeped_event, properties=properties, fields=fields, crs=crs, schema=schema))




def rename_feature(cat: TFCat, 
            time: float,
            frequency: float, 
            feature_type: str,
            new_feature_type: str):
    crs = cat.crs
    properties = cat.properties
    fields= cat.data.fields
    schema = cat.data["$schema"]
    event = FeatureCollection(cat.data, properties=properties, fields=fields, crs=crs, schema=schema).as_table()
    event['feature_type'][numpy.where((event['tmin'] <= time) & (event['tmax'] >= time)
                          & (event['fmin'] <= frequency) & (event['fmax'] >= frequency)
                          & (event['feature_type'] == feature_type))] = new_feature_type
    return TFCat(FeatureCollection().from_table(table=event, properties=properties, fields=fields, crs=crs, schema=schema))

#cut
#divide