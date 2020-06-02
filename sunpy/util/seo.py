"""
This module implements functions that can be used to match rows across various Data Catalogues.
All these functions will be coalesced into the Search Events Object.
"""

import pandas as pd
import numpy as np
from astropy.io.votable import parse
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
import warnings
from sunpy.util import SunpyUserWarning

class Sunspotter:
    """
    For the Sunspotter Dataset.
    """
    def __init__(self, *, rankings=None, classifications=None,
                 properties=None, timefits=None, delimiter=';'):
        """
        Parameters
        ----------
        rankings: `str`
            file path to the `rankings` CSV file.
        classifications: `str`
            file path to the `classifications` CSV file.
        properties: `str`
            file path to the `lookup_properties` CSV file.
        timefits: `str`
            file path to the `lookup_timefits` CSV file.
        delimiter : `str`
            Delimiter while reading the CSV files.
            Default delimiter is `;`
        """
        self.rankings = rankings
        self.classifications = classifications
        self.properties = properties
        self.timefits = timefits

        self.delimiter = delimiter

        self.get_dataframes()

    def get_dataframes(self):
        """
        Get dataframes from the sunspotter CSV files.
        """
        self.rankings = pd.read_csv(self.rankings, delimiter=self.delimiter)
        self.classifications = pd.read_csv(self.classifications, delimiter=self.delimiter)
        self.properties = pd.read_csv(self.properties, delimiter=self.delimiter)
        self.timefits = pd.read_csv(self.timefits, delimiter=self.delimiter)

    def get_timefits_id(self, obsdate):
        """
        Identifies the id for the observations for a given observation date.
        This id can be used to identify the particular observation across 
        different sunspotter datasets.

        Parameters
        ----------
        obsdate: `str` of format `yyyy-mm-dd hh-mm-ss`
            The date and time of a particular observation.

        Returns
        -------
        idx: `numpy.int64`
            The id for the observations for a given observation date.
        """
        return self.timesfits[self.timesfits.obs_date == obsdate].get(key='#id').iloc[0]

    def get_properties(self, idx):
        """
        Gets the list of properties of all ARs corresponding to a particular id.

        Parameters
        ----------
        idx: `numpy.int64`
            The id corresponding to a particular observation.

        Returns
        -------
        properties_list: `pd.DataFrame`
            DataFrame corresponding to all the observations for a given id.
        """
        return self.properties[self.properties.id_filename == idx]


class HFC:
    """
    For the HELIO Feature Catalogue.
    """
    def __init__(self, vot_file=None):
        """
        Parameters
        ----------
        vot_file: `str`
            file path to the vot file corresponding to the given observation.
        
        #TODO: Make getting the VOT table completely online using urllib.
        """
        self.vot_table = vot_file

        self.vot_table_to_pandas()

    def vot_table_to_pandas(self):
        """
        Reads the VOT table and returns the corresponding pandas DataFrame.

        Returns
        -------
        vot_table: `pd.Dataframe`
            Pandas Dataframe corresponding to the VOT table from HFC.
        """
        self.vot_table = parse(self.vot_table)
        table = self.vot_table.get_first_table().to_table(use_names_over_ids=True)
        return table.to_pandas()

