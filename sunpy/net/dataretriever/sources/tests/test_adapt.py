from sunpy.net import attrs
from sunpy.net.dataretriever.sources.adapt import carrington_time
from sunpy.net.dataretriever.attrs.adapt import *
import os

def test_client(CR=2193, frames=1, ask=False):
    res = test_search(CR, frames)
    print("Found {} files".format(len(res)))
    print(res)
    # ask the user to continue or not
    if len(res) > 0:
        if ask:
            print("Do you want to continue?")
            inp = input("y/n: ")
            if inp == 'y':
                out = test_fetch(res)
            else:
                print("Aborting")
                return None
        else:
            out = test_fetch(res)
    print(out)
    return out

def test_search(CR=2193, frames=1):
    from sunpy.net.fido_factory import Fido
    # Get the date range
    get_date, get_date_end = carrington_time(CR, frames)
    LngType = '0' # 0 is carrington, 1 is central meridian
    res = Fido.search(attrs.Instrument('adapt'), attrs.Time(get_date, get_date_end), ADAPTLngType(LngType))
    print(res)
    return res

def test_fetch(res, path="data"):
    from sunpy.net.fido_factory import Fido

    directory = os.path.join(os.getcwd(), path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    ret =  Fido.fetch(res, path=directory)
    test_download(ret)
    return ret


def test_download(out):
    print(out)
    pass


# print(__name__)
if __name__ == "__main__":
    # import os
    # os.execv('/usr/bin/python3', ['python3'] + sys.argv)
    # /opt/homebrew/anaconda3/envs/sunpy-dev/bin/python
    pass
    test_client()