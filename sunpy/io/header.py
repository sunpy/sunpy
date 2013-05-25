# -*- coding: utf-8 -*-
from __future__ import absolute_import
from collections import OrderedDict

class FileHeader(OrderedDict):
    def __init__(self, *args, **kwargs):
        OrderedDict.__init__(self, *args, **kwargs)