-------
Spectra
-------

Spectrograms
------------
SunPy currently supports reading dynamic spectra from Callisto instruments.
The main class that is used for this is
:py:class:`sunpy.spectra.sources.callisto.CallistoSpectrogram`. SunPy also
comes with an example image that shows a radio burst observerd by BIR that
can be found in sunpy.CALLISTO_IMAGE

    >>> from matplotlib import pyplot as plt
    >>> import sunpy
    >>> from sunpy.spectra.sources.callisto import CallistoSpectrogram
    >>> image = CallistoSpectrogram.read(sunpy.CALLISTO_IMAGE)

You can now view the image by using the :py:meth:`show` method.

    >>> image.show()

.. image:: ../images/spectra_ex1.png

You can then perform automatic constant background subtraction by using the
:py:meth:`subtract_bg` method. The resulting image will be clipped at 0
using the min_ parameter of show in order to avoid negative values.

    >>> nobg = image.subtract_bg()
    >>> nobg.show(min_=0)

.. image:: ../images/spectra_ex2.png

If you want to see the background determined by the automatic subtraction,
you can use the :py:meth:`auto_const_bg` method and visualize the resulting
data using :py:func:`pyplot.plot`.

    >>> bg = image.auto_const_bg()
    >>> plt.plot(image.freq_axis, bg)
    >>> plt.show() # This might not be necessary if you are using pylab.

.. image:: ../images/spectra_ex3.png

Now let us say we want to isolate the interesting bit (which starts around
10:38) from the boring background; there is a method called
:py:meth:`in_interval` that allows us to take the part of an image that is
within a specified interval. Leaving out the second interval border defaults
to the end.

    >>> interesting = nobg.in_interval("10:38")
    >>> interesting.show(min_=0)

.. image:: ../images/spectra_ex4.png

To get rid of the noise, we could also clip low intensities.

    >>> interesting.show(min_=20)

.. image:: ../images/spectra_ex5.png

