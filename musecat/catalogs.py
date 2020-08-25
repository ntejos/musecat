from astropy.table import Table
"""Module for Catalog classes"""

class Cat(object):
    """Class to handle a generic MUSE catalog"""
    def __init__(self, tab):
        self.tab = tab  # astropy table object

        #check if it has the minimum column names
        self._check_minimum_req(self.tab)

    def _check_minimum_req(self, tab):
        """Checks for minimum requirements before initializing"""
        colnames = tab.colnames
        min_req = ['RA_d', 'Dec_d', 'z']
        for req in min_req:
            if req not in colnames:
                raise ValueError('`{}` columname must be in catalog!'.format(req))



class CatLSDCat(Cat):
    """Class to handle main outputs of LSDCat"""

    def __init__(self, tab):
        super(CatLSDCat, self).__init__(tab)

    @classmethod
    def from_file(cls, filename):
        """Init from a filename"""
        tab = Table.read(filename)
        return CatLSDCat(tab)


class CatPyMUSECat(Cat):
    """Class to handle main outputs of XXX"""

    def __init__(self, tab):
        super(CatPyMUSECat,self).__init__(tab)

    @classmethod
    def from_file(cls, filename):
        """Init from a filename"""
        tab = Table.read(filename)
        return CatPyMUSECat(tab)



