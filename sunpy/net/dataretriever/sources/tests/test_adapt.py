import pytest
from sunpy.net import attrs
from sunpy.net.dataretriever.sources import adapt as adapt_source
from sunpy.net.dataretriever.attrs.adapt import ADAPTLonType
import os
from sunpy.net.fido_factory import Fido

global result
result = None

@pytest.mark.remote_data
def test_search(CR=2193, frames=1, test=True):

    # Get the date range
    get_date, get_date_end = adapt_source.carrington_time(CR, frames)
    LngType = '0' # 0 is carrington, 1 is central meridian
    print (f"\n\n\n\n\n\n\nSearching for images from {get_date} to {get_date_end}\n")
    res = Fido.search(attrs.Instrument('adapt'), attrs.Time(get_date, get_date_end), ADAPTLonType(LngType))
    print(res)
    global result
    result = res
    length_test(result, test=test)
    return result

@pytest.mark.remote_data
def test_fetch(search_result):
    print(f"Fetching search_result:")
    ret =  Fido.fetch(search_result, path=make_dir())
    length_test(ret, test=True)
    print("Success!\n\n")


@pytest.fixture
def search_result():
    global result
    return result
    # return test_search()

# Helper Functions ####################################################
def length_test(ret, test=False):
    """Checks to see if the variable exists and has a length

    Args:
        ret (_type_): _description_
        test (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
    """
    if test:
        return None
    try:
        assert len(ret) > 0
    except TypeError:
        return None
    return ret

#pytest fixtures
def make_dir(folder="data"):
    directory = os.path.join(os.getcwd(), folder)
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory


def test_client(CR=2193, frames=1, ask=False):
    res = test_search(CR, frames)
    # print("Found {} files".format(len(res)))
    # print(f"Client res is {res}")
    # ask the user to continue or not
    assert res
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
    return out


if __name__ == "__main__":
    test_client()


    # from sunpy.tests.helpers import remote_data
# from pytest.mark import remote_data


