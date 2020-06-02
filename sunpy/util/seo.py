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


class TableMatcher:
    """
    Table Matcher Object for finding corresponding rows across two distinct datasets.
    """
    def __init__(self, df_1, df_2, match_type=None):
        """
        Parameters
        ----------
        df_1: `pd.DataFrame`
            First DataFrame to match the rows from.
        df_2: `pd.DataFrame`
            Second DataFrame to match the rows from.
        match_type: `str`
            Matching algorithm to match the rows of df_1 with rows of df_2.
        """
        self.df_1 = df_1
        self.df_2 = df_2

        self.match_type = match_type
        if self.match_type is None:
            self.match_type = 'cosine'
        
        elif self.match_type is not None and self.match_type != 'euclidean':
            raise SunpyUserWarning('Incorrect matching algorithm specified.')

    def prepare_tables(self, feature_1, feature_2):
        """
        Prepates tables for matching.

        Parameters
        ----------
        feature_1: `list`
            List of columns from df_1 to match the rows with.
        feature_2: `list`
            List of columns from df_2 to match the rows with.

        Returns
        -------
        df_1: `pd.DataFrame`
            df_1 with only those columns that will be used for matching. 
        df_2: `pd.DataFrame`
            df_2 with only those columns that will be used for matching. 

        #TODO: Add pairwise type checks to see if corresponding features
        are in the right alignment.
        """
        if len(feature_1) != len(feature_2):
            raise SunpyUserWarning("The number of columns to match the rows on must be the same.")
        
        try:
            df_1 = self.df_1[feature_1]
        except KeyError:
            raise SunpyUserWarning("The features specified for table 1 do not "
                                    "correspond to any columns in table 1.")
        try:
            df_2 = self.df_2[feature_2]
        except KeyError:
            raise SunpyUserWarning("The features specified for table 2 do not "
                                    "correspond to any columns in table 2.")

        return df_1, df_2

    def match_cosine(self, df_1, df_2):
        """
        Finds Cosine similarity between the rows of the two dataframes.

        Parameters
        ----------
        df_1: `pd.DataFrame`
            First DataFrame to match the rows from.
        df_2: `pd.DataFrame`
            Second DataFrame to match the rows from.
        
        Returns
        -------
        result: `numpy.ndarray`
            Array of size `(n,)` where n is the number of rows in df_1.
            Contains indices of rows from df_2 that best correspond to rows from df_1.
        match_score: `numpy.ndarray`
            Array of size `(n,)` where n is the number of rows in df_1.
            Contains match score for  corresponding best matches.
        """
        
        cosine = cosine_similarity(X=df_1, Y=df_2)
        result = np.argmax(cosine, axis=1)
        match_score = np.max(cosine, axis=1)

        return result, match_score

    def match_euclidean(self, df_1, df_2):
        """
        Finds euclidean distance between the rows of the two dataframes.

        Parameters
        ----------
        df_1: `pd.DataFrame`
            First DataFrame to match the rows from.
        df_2: `pd.DataFrame`
            Second DataFrame to match the rows from.

        Returns
        -------
        result: `numpy.ndarray`
            Array of size `(n,)` where n is the number of rows in df_1.
            Contains indices of rows from df_2 that best correspond to rows from df_1.
        match_score: `numpy.ndarray`
            Array of size `(n,)` where n is the number of rows in df_1.
            Contains match score for  corresponding best matches.
        """
        euclidean = euclidean_distances(X=df_1, Y=df_2)
        result = np.argmin(euclidean, axis=1)
        match_score = np.min(euclidean, axis=1)

        return result, match_score

    def verify(self, match_score, threshold):
        """
        Verify matching quality. If any match score is less than the threshold,
        raises Sunpy User Warnings.

        Parameters
        ----------
        match_score: `numpy.ndarray`
            Array of size `(n,)` where n is the number of rows in df_1.
            Contains match score for  corresponding best matches.
        threshold: `float`
            Minimum score for considering a proper match.
        """
        for index, score_value in enumerate(match_score):
            if score_value < threshold:
                warnings.warn(SunpyUserWarning(f"\nMatch at Index {index} is likely to be incorrect\n"))

    def match(self, feature_1, feature_2, threshold):
        """
        Finds best match between the rows of the two dataframes.
        Raises warning id matching is dubious.

        Parameters
        ----------
        feature_1: `list`
            List of columns from df_1 to match the rows with.
        feature_2: `list`
            List of columns from df_2 to match the rows with.
        threshold: `float`
            Minimum score for considering a proper match.

        Returns
        -------
        result: `numpy.ndarray`
            Array of size `(n,)` where n is the number of rows in df_1.
            Contains indices of rows from df_2 that best correspond to rows from df_1.
        """
        df_1, df_2 = self.prepare_tables(feature_1, feature_2)

        if self.match_type == 'cosine':
            result, match_score = self.match_cosine(df_1, df_2)
        elif self.match_type == 'euclidean':
            result, match_score = self.match_euclidean(df_1, df_2)

        self.verify(match_score, threshold)

        return result
