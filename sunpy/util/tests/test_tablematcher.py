import numpy as np
import pandas as pd
import pytest

from sunpy.util import SunpyUserWarning
from sunpy.util.tablematcher import TableMatcher


@pytest.fixture
def df_1():
	data = {
		'feat_a': [415.31357, 	788.04966, 	615.53375 ,	194.83354 ,	776.57878 ,	831.94859],
		'feat_b': [266.63121,	343.50743, 	357.68856 ,	386.77800 ,	662.85391 ,	739.92834],
		'feat_c': [-191.75614,	546.61535, 	204.86999 ,	-628.51599 ,	523.89210 ,	633.57692],
		'feat_d': [-484.68983, 	-332.40190, 	-304.30979, 	-246.68507 ,	300.20745 ,	452.88802],
		'feat_e': [28.0535, 	38.6332, 	18.4522, 	42.1355, 	40.4294, 	55.406],
		'feat_f': [121.139, 	23.7468, 	104.284, 	62.3888, 	41.531, 	39.5828]
		}
	return pd.DataFrame(data)


@pytest.fixture
def df_2():
    data = {
        'feat_aa': [693., 190., 593., 146., 755., 417., 478., 196., 106., 764., 150.,
                    389., 691., 816., 855., 616., 641., 679.],
        'feat_bb': [586., 205., 647., 546., 332., 268., 811., 388., 469., 662., 705.,
                    746., 435., 325., 374., 359., 721., 622.],
        'feat_cc': [357.563 , -637.184 ,  160.185 , -724.916 ,  481.672 , -189.093 ,
                    -67.6126, -625.524 , -805.133 ,  498.883 , -717.188 , -242.919 ,
                    354.999 ,  601.707 ,  678.369 ,  205.842 ,  254.527 ,  330.072],
        'feat_dd': [147.131 , -607.368 ,  268.166 ,   68.1751, -354.912 , -482.907 ,
                    593.771 , -244.81  ,  -84.6278,  297.731 ,  383.462 ,  464.929 ,
                    -151.913 , -369.151 , -271.767 , -301.903 ,  415.062 ,  219.31 ],
        'feat_xx': [567., 184., 632., 450., 281., 191., 795., 335., 452., 637., 686.,
                    727., 413., 312., 354., 303., 704., 606.],
        'feat_yy': [329.343994, -644.799988,  162.688004, -656.703979,  392.832001,
                    -319.424011,  -75.391998, -668.607971, -809.471985,  577.343994,
                    -716.223999, -242.048004,  343.231995,  599.16803 ,  676.544006,
                    255.936005,  255.936005,  327.359985],
        'feat_zz': [109.120003, -650.752014,  238.080002, -123.008003, -458.303986,
                    -636.864014,  561.471985, -351.167999, -119.040001,  248.      ,
                    345.216003,  426.559998, -196.416   , -396.799988, -313.471985,
                    -414.656006,  380.928009,  186.496002]
        }
    return pd.DataFrame(data)


@pytest.fixture
def feature_1():
	return ['feat_a', 'feat_b', 'feat_c', 'feat_d']


@pytest.fixture
def feature_2():
	return ['feat_aa', 'feat_bb', 'feat_cc', 'feat_dd']

@pytest.fixture
def best_match():
	return [5, 13, 15,  7,  9,  9]

def test_match_default_algorithm(df_1, df_2):
	matcher = TableMatcher(df_1, df_2)
	assert matcher.match_type == 'cosine'

def test_match_incorrect_algorithm(df_1, df_2):
	with pytest.raises(SunpyUserWarning):
		TableMatcher(df_1, df_2, match_type='notAnAlgorithm')

def test_match_cosine_no_threshold(df_1, df_2, feature_1, feature_2, best_match):
	matcher = TableMatcher(df_1, df_2)
	result = matcher.match(feature_1, feature_2)

	assert np.array_equal(result, best_match)

def test_match_cosine(df_1, df_2, feature_1, feature_2, best_match):
	matcher = TableMatcher(df_1, df_2)

	with pytest.raises(SunpyUserWarning):
		matcher.match(feature_1, feature_2, threshold=0.999)

def test_match_euclidean_no_threshold(df_1, df_2, feature_1, feature_2, best_match):
	matcher = TableMatcher(df_1, df_2, match_type='euclidean')
	result = matcher.match(feature_1, feature_2)

	assert np.array_equal(result, best_match)

def test_match_euclidean(df_1, df_2, feature_1, feature_2, best_match):
	matcher = TableMatcher(df_1, df_2, match_type='euclidean')

	with pytest.raises(SunpyUserWarning):
		matcher.match(feature_1, feature_2, threshold=3)
