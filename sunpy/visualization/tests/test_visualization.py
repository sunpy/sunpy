from matplotlib.figure import Figure

import astropy.units as u
from astropy.wcs import WCS

from sunpy.tests.helpers import figure_test
from sunpy.visualization import axis_labels_from_physical_type, show_hpr_impact_angle


def test_axis_labels_from_physical_type():
    assert axis_labels_from_physical_type('custom:pos.helioprojective.lon', None) == \
        'Helioprojective Longitude (Solar-X)'
    assert axis_labels_from_physical_type('custom:pos.helioprojective.lat', None) == \
        'Helioprojective Latitude (Solar-Y)'
    assert axis_labels_from_physical_type('custom:pos.heliographic.stonyhurst.lon', None) == \
        'Heliographic Longitude'
    assert axis_labels_from_physical_type('custom:pos.heliographic.carrington.lon', None) == \
        'Carrington Longitude'
    assert axis_labels_from_physical_type('custom:pos.helioprojectiveradial.lon', None) == \
        'Helioprojective Position Angle'
    assert axis_labels_from_physical_type('custom:pos.helioprojectiveradial.lat', None) == \
        'Helioprojective Declination'
    # Unknown type should fall back to the physical type string itself
    assert axis_labels_from_physical_type('custom:pos.unknown', None) == 'custom:pos.unknown'
    # None physical type returns empty string
    assert axis_labels_from_physical_type(None, None) == ''
    # Unit is appended when provided
    assert axis_labels_from_physical_type('custom:pos.helioprojective.lon', 'deg') == \
        'Helioprojective Longitude (Solar-X) [deg]'


@figure_test
def test_show_hpr_impact_angle():
    wcs = WCS({
        'ctype1': 'HRLN-CAR',
        'ctype2': 'HRLT-CAR',
        'cdelt1': 1,
        'cdelt2': 1,
        'crpix1': 180.5,
        'crpix2': 90.5,
        'crval1': 180,
        'crval2': 0,
    })

    fig = Figure(figsize=(9, 4))
    ax1 = fig.add_subplot(121, projection=wcs)
    ax1.grid()

    ax1.set_xlim(-0.5, 360 - 0.5)
    ax1.set_ylim(-0.5, 180 - 0.5)
    ax1.coords[0].set_ticks(spacing=90*u.deg)
    ax1.coords[1].set_ticks(spacing=30*u.deg)
    ax1.set_aspect(2)

    ax1.set_title("Shows declination by default")

    ax2 = fig.add_subplot(122, projection=wcs)
    ax2.grid()

    show_hpr_impact_angle(ax2.coords[1])

    ax2.set_xlim(-0.5, 360 - 0.5)
    ax2.set_ylim(-0.5, 180 - 0.5)
    ax2.coords[0].set_ticks(spacing=90*u.deg)
    ax2.coords[1].set_ticks(spacing=30*u.deg)
    ax2.set_aspect(2)

    ax2.set_title("Changed to show impact angle")

    return fig
