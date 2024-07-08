import asdf
import asdf_astropy

__all__ = ["write"]
def write(fname, data, header, **kwargs):
    # Create the ASDF file
        """
    Take data and header pairs and save it to asdf file 
    inspired from https://docs.sunpy.org/en/stable/generated/gallery/saving_and_loading_data/genericmap_in_asdf.html


    Parameters
    ----------
    fname : `str`
        File name, with extension.
    data : `numpy.ndarray`
        n-dimensional data array.
    header : `dict`
        A header dictionary.

    """

    meta = dict(header)
    with asdf.AsdfFile() as af:
        af.tree={"meta":meta,"data":data}

        af.write_to(fname)
