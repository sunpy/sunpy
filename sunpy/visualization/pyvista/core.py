import functools

from astropy.constants import R_sun
import numpy as np
import pyvista as pv

from sunpy.coordinates import HeliocentricInertial
from sunpy.map.maputils import all_corner_coords_from_map


__all__ = ['PyVistaPlotter']


class PyVistaPlotter:
    """
    A plotter for 3D data.

    Since pyvista is not coordinate aware, everything is converted to
    the `~sunpy.coordinates.HeliocentricInertial` frame, and in distance units
    of :math:`R_{sun} = 1`. The z-axis is aligned with the solar rotation axis.

    Attributes
    ----------
    plotter : pyvista.Plotter
    coordinate_frame : astropy.coordinates.BaseFrame
        Coordinate frame of the plot. The x, y, z axes of the pyvista plotter
        will be the x, y, z axes in this coordinate system.
    """
    def __init__(self, coordinate_frame=None):
        if coordinate_frame is None:
            coordinate_frame = HeliocentricInertial()
        self._coordinate_frame = coordinate_frame
        self.plotter = pv.Plotter()

    @property
    def coordinate_frame(self):
        """
        Coordinate frame of the plot.
        """
        return self._coordinate_frame

    @functools.wraps(pv.Plotter.show)
    def show(self, *args, **kwargs):
        """
        Show the plot.
        """
        self.plotter.show(*args, **kwargs)

    def _coords_to_xyz(self, coords):
        coords = coords.transform_to(self.coordinate_frame)
        coords.representation_type = 'cartesian'
        return np.column_stack((coords.x.to_value(R_sun),
                                coords.y.to_value(R_sun),
                                coords.z.to_value(R_sun)))

    def _pyvista_mesh(self, m):
        """
        Create a mesh from a map.

        Parameters
        ----------
        m : sunpy.map.Map
        """
        corner_coords = all_corner_coords_from_map(m)
        nodes = self._coords_to_xyz(corner_coords.ravel())
        grid = pv.StructuredGrid()
        grid.points = nodes
        grid.dimensions = [m.data.shape[0] + 1,
                           m.data.shape[1] + 1,
                           1]
        data = m.data.T.reshape(-1)
        grid['data'] = m.plot_settings['norm'](data)
        return grid

    def plot_map(self, m, **kwargs):
        """
        Plot a map.

        Parameters
        ----------
        m : sunpy.map.Map
            Map to be plotted.
        **kwargs :
            Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        cmap = kwargs.pop('cmap', m.cmap)
        mesh = self._pyvista_mesh(m)
        self.plotter.add_mesh(mesh, cmap=cmap, **kwargs)

    def plot_line(self, coords, **kwargs):
        """
        Plot a line from a set of coordinates.

        Parameters
        ----------
        coords : astropy.coordinates.SkyCoord
        **kwargs :
            Keyword arguments are handed to `pyvista.Plotter.add_mesh`.

        Notes
        -----
        This plots a `pyvista.Spline` object.
        """
        points = self._coords_to_xyz(coords)
        spline = pv.Spline(points)
        self.plotter.add_mesh(spline, **kwargs)

    def plot_solar_axis(self, length=2.5, arrow_kwargs={}, **kwargs):
        """
        Parameters
        ----------
        length : float
            Length of the arrow in multiples of solar radii.
        arrow_kwargs : dict
            Keyword arguments to be handed to `pyvista.Arrow`.
            ``start``, ``direction``, and ``scale`` cannot be manually
            specified, as they are automatically set.
        **kwargs :
            Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        defaults = {'shaft_radius': 0.01,
                    'tip_length': 0.05,
                    'tip_radius': 0.02}
        defaults.update(arrow_kwargs)
        arrow = pv.Arrow(start=(0, 0, -length / 2),
                         direction=(0, 0, length),
                         scale='auto',
                         **defaults)
        self.plotter.add_mesh(arrow, **kwargs)
