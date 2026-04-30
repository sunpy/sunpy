import json
import pathlib

from astroquery.utils.tap.core import TapPlus

soar = TapPlus(url="http://soar.esac.esa.int/soar-sl-tap/tap")


__all__ = ["update_soar_attr_data"]


def get_cdf_descriptors():
    # Get data descriptors for CDF files
    print("Updating CDF descriptors...")
    job = soar.launch_job("select * from soar.cdf_dataset")
    res = job.get_results()
    descriptors = {}
    for row in res:
        desc = row["logical_source"].split("_")[-1]
        descriptors[desc] = row["logical_source_description"]
    return descriptors


def get_fits_descriptors():
    # Get data descriptors for FITS files
    print("Updating FITS descriptors...")
    soar = TapPlus(url="http://soar.esac.esa.int/soar-sl-tap/tap")
    job = soar.launch_job("select * from soar.fits_dataset")
    res = job.get_results()
    descriptors = {}
    for row in res:
        desc = row["logical_source"].split("_")[-1]
        # Currently no way to get a description from the FITS table
        descriptors[desc] = ""
    return descriptors


def get_all_descriptors():
    desc = get_cdf_descriptors()
    desc.update(get_fits_descriptors())
    return desc


def get_all_instruments():
    # Get the unique instrument names
    SOAR = TapPlus(url="http://soar.esac.esa.int/soar-sl-tap/tap")
    job = SOAR.launch_job("select * from soar.instrument")
    res = job.get_results()

    instruments = [
        "EPD",
        "EUI",
        "MAG",
        "METIS",
        "PHI",
        "RPW",
        "SOLOHI",
        "SPICE",
        "STIX",
        "SWA",
    ]
    instr_desc = {}
    for r in res:
        if r["name"] not in instruments:
            pass
        else:
            instr_desc[r["name"]] = r["long_name"]
    return instr_desc


def get_all_sensors():
    """
    Get a dict of unique sensor names and descriptions from SOAR. Parse to json
    with something similar to:

        sensors = get_all_sensors()
        json.dumps(get_all_sensors(), indent=2)

    The contents of this output should be placed in
    ./sunpy_soar/data/sensor_attrs.json

    Returns
    -------
    sensor_names: Dict[str, str]
    """
    print("Retrieving sensor descriptors...")
    job = soar.launch_job("SELECT * FROM soar.sensor")
    res = job.get_results()
    sensor_names = {}
    for row in res:
        sensor_names[row["name"]] = str(row.get("long_name", ""))
    return sensor_names


def get_all_soops():
    # Get the unique soop names
    print("Updating SOOP descriptors...")
    SOAR = TapPlus(url="http://soar.esac.esa.int/soar-sl-tap/tap")
    job = SOAR.launch_job("select * from soar.soop")
    res = job.get_results()

    soop_names = {}
    for row in res:
        soop_names[row["soop_name"]] = ""

    return soop_names


def update_soar_attr_data():
    import sunpy.net.soar.data

    base_path = pathlib.Path(sunpy.net.soar.data.__file__).parent
    attr_file = base_path / "attrs.json"

    descriptors = get_all_descriptors()
    with attr_file.open("w") as attrs_file:
        json.dump(dict(sorted(descriptors.items())), attrs_file, indent=2)

    instr_file = base_path / "instrument_attrs.json"
    instr_descriptors = get_all_instruments()
    with instr_file.open("w") as instrs_file:
        json.dump(dict(sorted(instr_descriptors.items())), instrs_file, indent=2)

    sensor_file = base_path/ "sensor_attrs.json"
    sensor_descriptors = get_all_sensors()
    with sensor_file.open("w") as sensors_file:
        json.dump(dict(sorted(sensor_descriptors.items())), sensors_file, indent=2)

    soop_file = base_path/ "soop_attrs.json"
    soop_descriptors = get_all_soops()
    with soop_file.open("w") as soops_file:
        json.dump(dict(sorted(soop_descriptors.items())), soops_file, indent=2)
