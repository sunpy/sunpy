import numpy as np
import pandas as pd
import urllib2


def cdaw(url="http://cdaw.gsfc.nasa.gov/CME_list/UNIVERSAL/text_ver/univ_all.txt"):
    """
    Download the CDAW CME catalog text file and return it as a pandas
    DataFrame.

    Parameters
    ----------
    url : string
        The URL of the full CDAW CME catalog as a text file.

    :return: `~pandas.DataFrame`
        The CDAW CME catalog returned as pandas dataframe object.

    References
    ----------
    * `SOHO LASCO CME Catalog <http://cdaw.gsfc.nasa.gov/CME_list/>`_
    * `Catalog description <http://cdaw.gsfc.nasa.gov/CME_list/catalog_description.htm>`_
    """
    f = urllib2.urlopen(url)

    # Read in the header lines first.
    for i in range(4):
        f.readline()

    # Initialise the data list
    cols = ['date', 'time', 'cpa', 'width', 'lin_speed',
            'quad_speed_init', 'quad_speed_final','quad_speed_20',
            'accel', 'mass', 'kin_energy', 'mpa', 'remarks']
    source = {k: list() for k in cols}

    # Parse the text.
    for line in f:
        line = line.replace('Halo','360').strip()
        columns = line.split()
        columns = [col.strip('*') if ("--" not in col) and ("**" not in col) else np.NaN for col in columns]
        if len(columns) < 13:
            columns.append(' ')
        for i, (k,val) in enumerate(zip(cols,columns)):
            if i == 12:
                source[k].append(' '.join(columns[i:]))
            else:
                source[k].append(val)
    f.close()

    # Return a pandas dataframe
    return pd.DataFrame(source, columns=cols)


