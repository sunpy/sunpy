from lxml import etree
from zeep import Plugin

from sunpy import log

__all__ = ['SunPyLoggingZeepPlugin']


class SunPyLoggingZeepPlugin(Plugin):
    def ingress(self, envelope, http_headers, operation):
        log.debug("VSO Response:\n " + etree.tostring(envelope, pretty_print=True).decode("utf-8"))
        return envelope, http_headers

    def egress(self, envelope, http_headers, operation, binding_options):
        log.debug("VSO Request:\n " + etree.tostring(envelope, pretty_print=True).decode("utf-8"))
        return envelope, http_headers
