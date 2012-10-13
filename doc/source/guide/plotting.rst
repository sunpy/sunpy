-----------------
Plotting in SunPy
-----------------

SunPy makes use of `matplotlib <http://matplotlib.org/>` for all of its plotting needs 
as such tries to follow the matplotlib plotting philosophy. 
It is therefore useful to go over how matplotlib works as background.

1. Matplotlib Tutorial
----------------------
The tutorial provided here is a summary of one that can be found in the `matplotlib
usage documentation <http://matplotlib.org/faq/usage_faq.html/>`.

Matplotlib provides two main pathways for plotting. One is meant for interactive use
(e.g. command-line) and the other for non-interactive use (e.g. modules). It is important
to recognize though that the interactive-use pathway (referred to as pyplot) just
provides shortcuts doing many of the more advanced non-interactive functions in the 
background. It is therefore possible to switch between the two as necessary and
convenient and it is possible to use pyplot in non-interactive way. In this manner pyplot
is just a shortcut making it quicker to set up plot axes and figures. 
In order to get access to the full interactive capabilities of pyplot it is
necessary to turn this feature on. This will be discussed below.
Pylab is another matplotlib usage scenario but it is essentially just pyplot with the
interactive capabilities turned on and numpy and matplotlib imported into the main 
namespace.

2. Pyplot
---------
Here is an example of pyplot usage ::

    import matplotlib.pyplot as plt

    plt.plot(range(10), range(10))
    plt.title("A simple Plot")
    plt.show()

The show() command open a plot on the screen and blocks execution (meaning you can't 
do anything with the prompt or your script freezes) until the plot window is closed. For 
reasons that are not very clear, the creators of matplotlib designed the show() command
so that it would only work once. If you were to call show() on the plt object again 
after the above code is executed nothing happens. Apparently, this confusing behavior 
is something that the matplotlib devs get complaints about often and so this may change
in the future (or may already have changed depending on your choice of backend). 
A discussion about this can be found `here 
<http://stackoverflow.com/questions/5524858/matplotlib-show-doesnt-work-twice>`.
Don't be confused by another command called draw(). This is only used while in interactive
mode. 

To turn on interactivity for pyplot use the command ::
    
    plt.ion(). 
    
In interactive mode, the plot will appear at the first plot() command and most 
commands will update the plot as you call them. Here is an example ::

    plt.plot(range(10), range(10))
    plt.title("Simple Plot")
    
In this example, you'll see that the title appears right on the plot when you call it.
Note that in this case the show command is useless as the plot shows up right when you
create it. Also note that some commands will not automatically update the plot and
you have to use the draw() command. The following command ::

    plt.ioff()
    
turns off interactivity.

3. Advanced Pyplot
------------------
If you need more fine-grained control over plots the recommended path is to use pyplot
and access the figures and axes objects. Here is an example ::

    import matplotlib.pyplot as plt
    import numpy as np
    x = np.arange(0, 10, 0.2)
    y = np.sin(x)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    plt.show()

Figure is the top-level container for all plot elements and axes is the top-level container
for a particular plot. So the above example, creates a figure then creates an axes
and then uses the pyplot plot() method to populate the plot in ax. You generally don't need
to mess with the figure object but with this method you now have your hands on the ax
object so you can do this like change the labels on the x and y axes or add a legend, etc.
do whatever you want to it. In the previous section, pyplot took care of creating these
objects for you so you don't have to worry about creating them yourself.

4. SunPy Plotting Standards
---------------------------

To be consistent with matplotlib, SunPy has developed a standard plotting policy which 
supports both simple and advanced matplotlib usage. 

5. show()
---------

For quick and easy access to a plots
all sunpy base objects (e.g. map, spectra, lightcurve) define their own show() command.
For example you can do the following ::

    import sunpy
    map = sunpy.make_map(sunpy.EIT_195_IMAGE)
    map.show()

This creates a plot window with all axes defined, a plot title, and the image of the map
data all defined by the contents of the map. As this is a show() command, the plot window
blocks and must be closed before doing anything else. This is meant as a quick way to 
visualize the contents of an object you've created.

6. plot()
---------

For more advanced plotting the base sunpy objects also provide a plot() command. This
command is similar to the pyplot plot() command in that it will create a figure and axes
object for you if you haven't already. It returns a figure object and does not create a
plot window. With the figure object in your hands you can reach in and grab the axes
and therefore manipulate the plot as you see fit. Here is an example of this at work ::

    import sunpy
    import matplotlib.pyplot as plt
    map = sunpy.make_map(sunpy.EIT_195_IMAGE)
    fig = map.plot()
    plt.show()

This output of this example is equivalent to one in the previous section. If we want
to make changes to the plot that is possible by using the gca() command on the figure
object. This returns the axes object.    

    ax = fig.gca()
    ax.plot([-1000,1000], [0,0])
    plt.show()

The above a plot of line across the map. Using the fig.gca() command to get access to the
axes object most anything can be done to the plot and the plot can be displayed as usual
using the show() command. Here is another example ::

    from matplotlib import patches
    fig = map.plot()
    ax = fig.gca()
    rect = patches.Rectangle([-350, -650], 500, 500, color = 'white', fill=False)
    ax.add_artist(rect)
    plt.show()
    
Finally, here is a more complex example ::

    from matplotlib import patches
    import sunpy
    import matplotlib.pyplot as plt
    map = sunpy.make_map(sunpy.AIA_171_IMAGE)
    smap = map.submap([-100-250, -100+250], [-400-250, -400+250])
    rect = patches.Rectangle([-200, -400], 500, 500, color = 'white')
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    map.plot(ax1)
    ax1.add_artist(rect)
    ax2 = fig.add_subplot(112)
    smap.plot(ax2)
    plt.show()
    
The above example creates two side by side plots one with the overall view of the Sun
with a small area marked with a white box. That smaller view is then shown in the plot
next to it.