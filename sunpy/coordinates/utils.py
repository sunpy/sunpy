from sunpy.coordinates import Helioprojective, sun

__all__ = ['solar_angular_radius', 'coordinate_is_on_solar_disk']


def _verify_coordinate_helioprojective(coordinates):
    """
    Raises an error if the coordinate is not in the
    `~sunpy.coordinates.frames.Helioprojective` frame.

    Parameters
    ----------
    coordinates : `~astropy.coordinates.SkyCoord`, `~astropy.coordinates.BaseCoordinateFrame`
    """
    frame = coordinates.frame if hasattr(coordinates, 'frame') else coordinates
    if not isinstance(frame, Helioprojective):
        raise ValueError(f"The input coordinate(s) is of type {type(frame).__name__}, "
                         "but must be in the Helioprojective frame.")


def solar_angular_radius(coordinates):
    """
    Calculates the solar angular radius as seen by the observer.

    The tangent vector from the observer to the edge of the Sun forms a
    right-angle triangle with the radius of the Sun as the far side and the
    Sun-observer distance as the hypotenuse. Thus, the sine of the angular
    radius of the Sun is ratio of these two distances.

    Parameters
    ----------
    coordinates : `~astropy.coordinates.SkyCoord`, `~sunpy.coordinates.frames.Helioprojective`
        The input coordinate. The coordinate frame must be
        `~sunpy.coordinates.Helioprojective`.

    Returns
    -------
    angle : `~astropy.units.Quantity`
        The solar angular radius.
    """
    _verify_coordinate_helioprojective(coordinates)
    return sun._angular_radius(coordinates.rsun, coordinates.observer.radius)


@u.quantity_input
def coordinate_is_on_solar_disk(coordinates):
    """
    Checks if the helioprojective Cartesian coordinates are on the solar disk.

    The check is performed by comparing the coordinate's angular distance
    to the angular size of the solar radius. The solar disk is assumed to be
    a circle i.e., solar oblateness and other effects that cause the solar disk to
    be non-circular are not taken in to account.

    Parameters
    ----------
    coordinates : `~astropy.coordinates.SkyCoord`, `~sunpy.coordinates.frames.Helioprojective`
        The input coordinate. The coordinate frame must be
        `~sunpy.coordinates.Helioprojective`.

    Returns
    -------
    `~bool`
        Returns `True` if the coordinate is on disk, `False` otherwise.
    """
    _verify_coordinate_helioprojective(coordinates)
    # Calculate the radial angle from the center of the Sun (do not assume small angles)
    # and compare it to the angular radius of the Sun
    return np.arccos(np.cos(coordinates.Tx) * np.cos(coordinates.Ty)) <= solar_angular_radius(coordinates)
