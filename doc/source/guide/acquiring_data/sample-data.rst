-----------------------
Downloading Sample Data
-----------------------

SunPy provides a number of sample data files which are referenced in the 
documentation and examples. These files are available to download onto your 
local machine so that you can try out the code in the documentation. To 
download the sample data simply run the following command::

    import sunpy.data
    sunpy.data.download_sample_data()

This will download the data to your sample-data directory which can be 
customized by editing the sunpyrc file (see :doc:`../customization`). 
After the files have downloaded SunPy provides custom 
shortcuts to these files in the `sunpy.data.sample` module which you 
can access by simply import this module.
