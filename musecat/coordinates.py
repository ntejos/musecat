"""
Coordinate management
musecat 0.1.0
arctomo
C.Moya
"""

# basic astropy
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u

#third party
from musecat import catalogs as cats

def ds9_to_coord(ds9):
    pass


def read_coord(ra, dec):
    """

    :param ra: string of a RA
    :param dec: string of a Dec
    :return point: SkyCoord

    Create a SkyCoord from a tuple of strings RA,Dec

    """
    if ":" in ra and ":" in dec:
        units = [u.hourangle, u.deg]
    elif " " in ra and " " in dec:
        units = [u.hourangle, u.deg]
    else:
        try:
            ra, dec = float(ra), float(dec)
            units = [u.deg, u.deg]
        except ValueError:
            print(
                "Unable to parse coordinates\n available formats are: hh:mm:ss.s dd:mm:ss.s\n "
                "hh mm ss.s dd mm ss.s\n ddd.dd dd.dd")
            return
    point = SkyCoord(ra, dec, unit=units)
    return point


class Matcher:
    def __init__(self,cat1,cat2):
        self.catpy = None
        self.catLSD = None

        #assing variables
        if type(cat1)== cats.CatPyMUSECat:
            self.catpy = cat1

        elif type(cat1) == cats.CatLSDCat:
            self.catLSD = cat1

        else:
            print("wrong input types")
            return


        if type(cat2) == cats.CatPyMUSECat:
            self.catpy = cat2

        elif type(cat1) == cats.CatLSDCat:
            self.catLSD = cat2

        else:
            print("wrong input types")
            return

        coords1 = []
        coords2 = []

        for i in self.catpy.tab:
            ra,dec  = str(i['RA_d']),str(i['Dec_d'])
            s = read_coord(ra,dec)
            coords1.append(s)

        for i in self.catLSD.tab:
            ra, dec = str(i['RA_d']), str(i['Dec_d'])
            s = read_coord(ra, dec)
            coords2.append(s)

        self.pycoords = SkyCoord(coords1)
        self.LSDcoords = SkyCoord(coords2)

    def match(self):
        """
        Returns the indices and distance in arcsec  of LSD sources on the PyMUSE catalog
        :return:
        """
        #max_sep = 1.0 * u.arcsec
        # idx, d2d, d3d = cc1.match_to_catalog_3d(cc2)
        idx, d2d, d3d = self.pycoords.match_to_catalog_sky(self.LSDcoords)
        #sep_constraint = d2d < max_sep


        return idx,d2d.arcsec


