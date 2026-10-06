import numpy
from astropy.table import Table

from tfcat import TFCat
from tfcat_tools.featureaddons import Feature, FeatureCollection


class Monitor(TFCat):
    """
        Object to monitor (indexing, sorting, errors, missspelling)
    """

    def check(self, date_start, date_end, feature_type) :
        index = numpy.array(list(range(len(self.data.features))))
    def clean_errors(self): 


    def check_errors(self):
        table_event = self.data.as_table()
        

    def _spelling_error(self, table_event) :
        index = numpy.array(list(range(len(self.data.features))))
        for key in self.data.fields :
            if 'values' in self.data.fields[key]:
                values = list(self.data.fields[key]['values'].keys())
        return index[not table_event[key] in values]

    def correct(self, )
    


