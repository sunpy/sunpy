=====================
SSW/SunPy Cheat Sheat
=====================

`SolarSoft (SSW) <http://sohowww.nascom.nasa.gov/solarsoft/>`_ is a  
popular IDL software library for solar data analysis, and in fact, many parts 
of SunPy are inspired by data structions and functions in SSW.

In order to help make it easier for users with SSW experience to work with
SunPy, we have put together a list of tips for new Python/SunPy users, as well
as a table with the names of some functions with comparable IDL/SSW 
counterparts.

NOTE: This article is currently a work-in-progress. 

Caveats
-------
Some things to be careful about when switching from IDL/SSW to Python/SunPy.

* Copy by value vs. reference (Jack?)
* Top-left vs. bottom-left origins.
* etc.

Mapping of terms
----------------
Below is a table which describes some of the functions and classes available
in IDL/SSW and the names and locations of comparable versions in Python/SunPy.

**Functions**

=========  ========   ============= 
 IDL/SSW    SunPy      Module 
=========  ========   ============= 
congrid    resample   sunpy.image 
=========  ========   =============

