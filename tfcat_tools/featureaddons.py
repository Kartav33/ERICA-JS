"""
This class is adapted from: https://gitlab.obspm.fr/maser/catalogues/tfcat/-/blob/master/tfcat/feature.py
"""



from tfcat import TFCat
from tfcat.crs import DefaultCRS, CRS
from tfcat.validate import JSONSCHEMA_URI
from tfcat.feature import Feature, FeatureCollection, _FeatureColumn
from astropy.table import Table, Column

class FeatureCollection(FeatureCollection):
    """Class for a TFcat feature collection.

    This class is derived from the tfcat FeatureCollection class, adding new capabilities to
    manage data format and storage.

    :param features: List of features to constitute the FeatureCollection.
    :type features: iterable
    :param properties: Feature collection global properties
    :type properties: Union[dict, None]
    :param fields: Fields defining the feature properties
    :type fields: Union[list[Field], None]
    :param crs: Coordinate reference system
    :type crs: Union[dict, CRS, None]
    :param schema: TFCat JSON Schema URI (default to the current supported JSON Schema)
    :type schema: str
    :return: FeatureCollection object
    :rtype: FeatureCollection
    """

    @classmethod
    def from_table(cls, table, properties=None, fields=None, crs=None, schema=JSONSCHEMA_URI):
        feature  = [Feature(id=i, geometry=table[i]['feature'], properties={key:table[i][key] for key in fields}) for i in range(len(table))]
        return cls(feature, properties=properties, fields=fields, crs=crs, schema=schema)



    def as_table(self, features=None, properties=None, fields=None, crs=None, schema=JSONSCHEMA_URI):

        cols = [
                _FeatureColumn(
                    name='feature',
                    data=[item.geometry for item in self.features],
                    description='Feature',
                    crs=self.crs
                ),
                Column(
                    name='tmin',
                    data=[self.crs.time_converter(item.bbox[0]).unix for item in self.features],
                    description='Feature Start Time'
                ),
                Column(
                    name='tmax',
                    data=[self.crs.time_converter(item.bbox[2]).unix  for item in self.features],
                    description='Feature End Time'
                ),
                Column(
                    name='fmin',
                    data=[item.bbox[1] for item in self.features],
                    unit=self.crs.properties['spectral_coords']['unit'],
                    description='Feature lower spectral bound'
                ),
                Column(
                    name='fmax',
                    data=[item.bbox[3] for item in self.features],
                    unit=self.crs.properties['spectral_coords']['unit'],
                    description='Feature upper spectral bound'
                )
        ]

        for key in self.fields:
            cur_col = Column(
                name=key,
                dtype=self.fields[key]['datatype'],
                data=[item.properties[key] for item in self.features],
                unit=self.fields[key].get('unit', None),
                description=self.fields[key].get('info', None)
            )
            if 'values' in self.fields[key]:
                cur_col.meta = {'options': self.fields[key].get('values', None)}
            cols.append(cur_col)

        return Table(cols)

