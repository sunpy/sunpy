import asdf
import asdf_astropy

__all__ = ["write"]
def write(fname, data, header, **kwargs):
    # Create the ASDF file
    "This function will write a asdf file "

    meta = dict(header)
    with asdf.AsdfFile() as af:
        af.tree={"meta":meta,"data":data}

        af.write_to(fname)
