"""
==================================
scrapping files with bs4 and urlib
==================================
"""
#############################imports################################
# bs4 is the initials to BeautifulSoup4, a powerful web scraper tool
import bs4
# urlib it the default python library to get content from the web
import urllib.parse
import urllib.request

##############################Code######################################

# this is a example url for the scrapping
base_url = "https://umbra.nascom.nasa.gov/pub/lasco/lastimage/200722/c2/"
# open the url with urlib
index = urllib.request.urlopen(base_url)
# indicate to the bs4 to read the open url
# there will be a warning saying about the parser choice, but it can be ignored
tree = bs4.BeautifulSoup(index.read())
# let's define a simple function to scape the files we want and send them with a generator


def retrive_files():
    # find all the occurencies of the tag 'a' in source url
    for a in tree.find_all("a"):
        # get the href porperty and assign to a different variable
        i = a.get("href")
        # and return as generator if its ends with .fts, in that case,it is the file we want
        if i.endswith("fts"):
            # for performace of the software, in the case of a lot of files, or huge ones, to send them one by one
            yield i
    # in a shorter version without generator
    # files = [a for a in [a.get("href") for a in tree.find_all("a")] if a.endswith("fts")]


# here the files are comming according to the time, we scraped them
files = retrive_files()
for filename in files:
    # here you can see the print of all files, one by one, to check if you get all
    print(filename)
# let's save each file in temp to use it later
for filename in files:
    urllib.request.urlretrieve(urllib.parse.urljoin(base_url, filename), f"/tmp/{filename}")

# for more information, read the docs
# Beautiful Soup 4
# https://www.crummy.com/software/BeautifulSoup/bs4/doc/

# urlib
# https://docs.python.org/3/library/urllib.html

#in case you need to interact with the Javascript of the page, you will need the selenium library
# Selenium
# https://www.selenium.dev/documentation/en/