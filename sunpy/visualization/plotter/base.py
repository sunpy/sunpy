import abc


class BasePlotter(abc.ABC):
    """
    Base plotter GenericMap plotter objects
    """

    @property
    @abc.abstractmethod
    def cmap(self):
        ...

    @abc.abstractmethod
    def draw_grid(selfd, *args, **kwargs):
        ...

    @abc.abstractmethod
    def draw_limb(selfd, *args, **kwargs):
        ...

    @abc.abstractmethod
    def draw_quadrangle(selfd, *args, **kwargs):
        ...

    @abc.abstractmethod
    def draw_contours(selfd, *args, **kwargs):
        ...

    @abc.abstractmethod
    def peek(selfd, *args, **kwargs):
        ...

    @abc.abstractmethod
    def plot(selfd, *args, **kwargs):
        ...
