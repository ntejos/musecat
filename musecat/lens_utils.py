from astropy.cosmology import Planck15 as cosmo
import numpy as np


class NumericalLensModel:

    """
    class to handle a numerical lens-model
    """

    def __init__(self, z_lens, z_source, alpha_x_array, alpha_y_array, wcs, cosmo=cosmo, mag_array=None):

        """
        Parameters
        ----------
        z_lens : float
            redshift of the modeled lens
        z_source : float
            redshift used to normalize the deflection matrices
        alpha_x_array : numpy array or list of list (2D matrix)
            deflections matrix in x axis direction, normalized with z_source redshift. Deflections in units of arcsec
        alpha_y_array : numpy array or list of list (2D matrix)
            deflections matrix in y axis direction, normalized with z_source redshift. Deflections in units of arcsec
        wcs : PyMuse WCS object (at the moment) (it could be changed to be an astropy wcs instead)
            wcs object to make conversions between pixel positions and global coordinates and vice-versa
        cosmo : astropy.cosmology.core.Cosmology objetc
            cosmology to make distance calculations based on redshift

        """
        self.z_lens = z_lens
        self.z_source = z_source
        self.alpha_x_array = alpha_x_array  # in arcsec
        self.alpha_y_array = alpha_y_array
        self.wcs = wcs  # PYMUSE wcs
        self.cosmo = cosmo
        self.mag_array = mag_array

    def dls(self, z_source):

        """
        Angular diameter distance between the lens redshift and a source at redshift 'z_source'

        Parameters
        ----------
        z_source : float

        Returns
        -------
        float
            distance in Mpc
        """

        return self.cosmo.angular_diameter_distance_z1z2(self.z_lens, z_source)

    def ds(self, z_source):

        """
        Angular diameter distance of a source at redshift 'z_source'

        Parameters
        ----------
        z_source : float

        Returns
        -------
        float
            distance in Mpc

        """

        return self.cosmo.angular_diameter_distance(z_source)

    def _alpha_x(self, ra, dec):

        """
        returns the value of the deflection matrix 'alpha_x_array' in the world coordinates 'ra', 'dec'.
        if pixel coordinates corresponding to 'ra', 'dec' coordinates are not integers, decimals are truncated

        Parameters
        ----------
        ra: float
            ra coordinate in degrees
        dec: float
            dec coordinate in degrees
        Returns
        -------
        float
            value of 'alpha_x_array' in the x, y positions corresponding to 'ra', 'dec'

        """

        y, x = self.wcs.sky2pix([dec, ra])[0] #ra, dec --> x, y
        x, y = int(x), int(y) # truncate decimals

        return self.alpha_x_array[y][x]

    def _alpha_y(self, ra, dec):

        """
        returns the value of the deflection matrix 'alpha_y_array' in the world coordinates 'ra', 'dec'.
        if pixel coordinates corresponding to 'ra', 'dec' coordinates are not integers, decimals are truncated

        Parameters
        ----------
        ra: float
            ra coordinate in degrees
        dec: float
            dec coordinate in degrees
        Returns
        -------
        float
            value of 'alpha_x_array' in the x, y positions corresponding to 'ra', 'dec'

        """

        y, x = self.wcs.sky2pix([dec, ra])[0] #ra, dec --> x, y
        x, y = int(x), int(y) # truncate decimals

        return self.alpha_y_array[y][x]

    def alpha_x_norm(self, ra, dec, z_source):

        """
        scales the deflections by redshit in x axis.
        Assumes that self.alpha_x_array and self.alpha_y_array are normalized to self.z_source redshift such that

        alpha_matrix = original_alpha*dls/ds(self.x_source).

        The normalizations is done by removing the old normalization and making a new in the following way

        new_alpha = alpha_matrix/(dls/ds(self.x_source))*(dls/ds(z_source))

        where dls and ds are angular diameter distances between the lens and the source, and the observer and the lens
        respectively (calculated using self.cosmo)


        Parameters
        ----------
        ra: float
            ra coordinate in degrees
        dec: float
            dec coordinate in degrees
        z_source: float
            redshift of the source

        Returns
        -------
        float
            value of 'alpha_x_array' deflection in pixels in the x, y positions corresponding to 'ra', 'dec' position
            scaled by redshift.

        """


        dls1 = self.dls(self.z_source)
        ds1 = self.ds(self.z_source)

        dls2 = self.dls(z_source)
        ds2 = self.ds(z_source)

        arcsec_deflection = self._alpha_x(ra, dec) / (dls1 / ds1) * (dls2 / ds2)
        pix_scale = self.wcs.get_step()[0] * 3600

        pix_deflection = arcsec_deflection/pix_scale

        return pix_deflection

    def alpha_y_norm(self, ra, dec, z_source):
        """
        scales the deflections by redshit in y axis.
        Assumes that self.alpha_x_array and self.alpha_y_array are normalized to self.z_source redshift such that

        alpha_matrix = original_alpha*dls/ds(self.x_source).

        The normalizations is done by removing the old normalization and making a new in the following way

        new_alpha = alpha_matrix/(dls/ds(self.x_source))*(dls/ds(z_source))

        where dls and ds are angular diameter distances between the lens and the source, and the observer and the lens
        respectively (calculated using self.cosmo)


        Parameters
        ----------
        ra: float
            ra coordinate in degrees
        dec: float
            dec coordinate in degrees
        z_source: float
            redshift of the source

        Returns
        -------
        float
            value of 'alpha_y_array' deflection in pixels in the x, y positions corresponding to 'ra', 'dec' position
            scaled by redshift.

        """

        dls1 = self.dls(self.z_source)
        ds1 = self.ds(self.z_source)

        dls2 = self.dls(z_source)
        ds2 = self.ds(z_source)

        arcsec_deflection = self._alpha_y(ra, dec) / (dls1 / ds1) * (dls2 / ds2)
        pix_scale = self.wcs.get_step()[0] * 3600

        pix_deflection = arcsec_deflection / pix_scale

        return pix_deflection

    def source_pos(self, ra, dec, z_source):

        """
        Applies the lens equation to calculate the source position in the source plane, given the coordinates in the
        image plane using the deflections matrices.
        The lens equation is

        beta_x = theta_x - alpha_x
        beta_y = theta_y - alpha_y

        here, theta_x and theta_y are the coordinates in the image plane, alpha_x, alpha_y the deflections scaled by
        redshift, and beta_x, beta_y are the coordinates of the source in the source plane+

        Parameters
        ----------
        ra: float
            ra coordinate in degrees in the image plane
        dec: float
            dec coordinate in degrees in the image plane
        z_source: float
            redshift of the source

        Returns
        -------
        float
            ra, dec coordinates of the source in the source plane
        """

        alpha_x, alpha_y = self.alpha_x_norm(ra, dec, z_source), self.alpha_y_norm(ra, dec, z_source)
        theta_y, theta_x = self.wcs.sky2pix((dec, ra))[0]
        beta_x, beta_y = theta_x - alpha_x, theta_y - alpha_y

        return self.wcs.pix2sky((beta_y, beta_x))[0][::-1]