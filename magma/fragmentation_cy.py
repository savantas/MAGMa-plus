# 2014.10.20 12:27:50 CEST

def __bootstrap__():
    global __file__
    global __bootstrap__
    global __loader__
    import sys
    import pkg_resources
    import imp
    __file__ = pkg_resources.resource_filename(__name__, 'fragmentation_cy.so')
    __loader__ = None
    del __bootstrap__
    del __loader__
    imp.load_dynamic(__name__, __file__)


__bootstrap__()
