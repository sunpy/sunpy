import glymur
from scipy.ndimage import zoom


def resize_jp2(jp2file: str, scale: float) -> None:
    """
    Resize a JPEG2000 file by a given scale factor.

    It uses zoom from scipy.ndimage to resize the image data, so factor is the inverse.

    It will copy the XML metadata from the original file to the new one.
    It will not overwrite the original file.
    """
    original_jp2 = glymur.Jp2k(jp2file)
    rescaled_data = zoom(original_jp2[:], scale)
    new_file = glymur.Jp2k(f'{jp2file}_rescale.jp2', rescaled_data)
    for box in original_jp2.box:
        if isinstance(box, glymur.jp2box.XMLBox):
            new_file.append(box)
