def _test_private_date_setters(map_object):
    map_object._set_date("2012-03-04 05:06:07")
    assert map_object.date.isot == "2012-03-04T05:06:07.000"

    map_object._set_reference_date("2012-03-04 05:06:07")
    assert map_object.reference_date.isot == "2012-03-04T05:06:07.000"
