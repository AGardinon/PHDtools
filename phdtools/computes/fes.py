# -------------------------------------------------- #
# Computes - fes module
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import numpy as np
import phdtools.plots as phdplot
from typing import Union, List
from phdtools.computes import misc

# -------------------------------------------------- #
# --- FES

class BaseFES:
    """Base FES class
    """

    # class variable for units conversions
    kb_units = dict(
        kJ = 0.00831446261815324,
        kcal = 0.00198720425864083,
        kb =  1.0
    )


    def __init__(self, 
                 temperature: Union[int, float],
                 units: str = 'kb'):
        
        self._temp = temperature
        if units in any(self.kb_units.keys()):
            self._unit = units
        else:
            raise NameError("Invalid unit type.\n"
                            f"Types available: {self.kb_units.keys()}")
        pass

    @property
    def temp(self) -> Union[int, float]:
        """Temperature in K.

        :return: temperature used for the FES calculation.
        :rtype: Union[int, float]
        """
        return self._temp
    
    @temp.setter
    def temp(self, 
             value: Union[int, float]) -> None:
        """Set Temperature.

        :param value: temperature in K.
        :type value: Union[int, float]
        :raises ValueError: if temperature is negative.
        """
        if value >= 0:
            self._temp = value
        else:
            raise ValueError("Temperature must be positive")
        pass

    @property
    def unit(self) -> str:
        """Energy unit.

        :return: energy unit.
        :rtype: str
        """
        return self._unit

    @unit.setter
    def unit(self, 
             value: str) -> None:
        """Set the energy unit.

        :param value: energy unit.
        :type value: str
        """
        if value in any(self.kb_units.keys()):
            self._unit = value
        else:
            raise NameError("Invalid unit type.\n"
                            f"Types available: {self.kb_units.keys()}")
        pass

    @property
    def kbT(self) -> float:
        """Computes the kbT constant.

        :return: kbT value in the chosen energy units
        :rtype: float
        """
        if self._unit == 'kb':
            return self._unit
        else:
            return self._unit * self._temp
        
# -------------------------------------------------- #

def histo_to_fes(histo: np.ndarray, 
                 kbt: float, 
                 zero_level: Union[str, float] ='min',
                 fill_empty=True) -> np.ndarray:
    """_summary_

    :param histo: _description_
    :type histo: np.ndarray
    :param kbt: _description_
    :type kbt: float
    :param zero_level: _description_, defaults to 'min'
    :type zero_level: Union[str, float], optional
    :param fill_empty: _description_, defaults to True
    :type fill_empty: bool, optional
    :raises ValueError: _description_
    :return: _description_
    :rtype: np.ndarray
    """
    zeta = -1 * kbt * np.log(histo)
    # scaling
    if isinstance(zero_level, str):
        if zero_level == 'min':
            zeta = zeta - np.min(zeta)
        elif zero_level == 'max':
            zeta = zeta - np.max(zeta)
    elif isinstance(zero_level, (float, int)):
        zeta = zeta - zero_level
    else:
        raise ValueError("Unrecognized `zero_level` value.\n"
                         "Value acceppted: 'max', 'min', numerical.")
    # fill empty
    if fill_empty:
        _max = np.max(zeta[zeta != np.inf])
        zeta[zeta == np.inf] = _max
    else:
        pass

    return zeta


def mesh_grid_2d(X: np.ndarray,
                 Y: np.ndarray, 
                 bins: int) -> (np.ndarray, np.ndarray):
    """_summary_

    :param X: _description_
    :type X: _type_
    :param np: _description_
    :type np: _type_
    :return: _description_
    :rtype: _type_
    """
    # boundaries - first dimension
    x_min, x_max = X.min(), X.max()
    xx = np.arange(x_min, x_max, 
                   ((x_max - x_min) / bins))
    # boundaries - second dimension
    y_min, y_max = Y.min(), Y.max()
    yy = np.arange(y_min, y_max, 
                   ((y_max - y_min) / bins))
    # mesh
    XX, YY = np.meshgrid(xx, yy)
    
    return XX, YY

# ---

class FES(BaseFES):
    """Computes the FES for a give set of data.

    :param BaseFES: base class
    :type BaseFES: class
    """

    def __init__(self, 
                 temperature: Union[int, float], 
                 units: str = 'kb'):
        super().__init__(temperature, units)
        self.fes_dict = None
        pass

    @misc.my_timer
    def fit(self,
            X: np.ndarray,
            Y: np.ndarray,
            bins: int,
            range: (float, float) =None,
            zero_level: Union[str, float] ='min',
            weights: np.ndarray =None,
            fill_empty=True) -> dict:
        """_summary_

        :param X: _description_
        :type X: np.ndarray
        :param Y: _description_
        :type Y: np.ndarray
        :param bins: _description_
        :type bins: int
        :param range: _description_, defaults to None
        :type range: float, float, optional
        :param zero_level: _description_, defaults to 'min'
        :type zero_level: Union[str, float], optional
        :param weights: _description_, defaults to None
        :type weights: np.ndarray, optional
        :param fill_empty: _description_, defaults to True
        :type fill_empty: bool, optional
        :return: _description_
        :rtype: dict
        """
        
        if not Y:
            print("Computing 1D ...")
            self.fes_dict = self._fit1D
            return self.fes_dict
        if Y:
            print("Computing 2D ...")
            self.fes_dict = self._fit2D
            return self.fes_dict
        else:
            raise ValueError("There are some problems.")


    def _fit1D(self, 
               X: np.ndarray,
               bins: int,
               range: (float, float) =None,
               zero_level: Union[str, float] ='min',
               weights: np.ndarray =None,
               fill_empty=True) -> dict:
        """_summary_

        :param X: _description_
        :type X: np.ndarray
        :param bins: _description_
        :type bins: int
        :param range: _description_, defaults to None
        :type range: float, float, optional
        :param zero_level: _description_, defaults to 'min'
        :type zero_level: Union[str, float], optional
        :param weights: _description_, defaults to None
        :type weights: np.ndarray, optional
        :param fill_empty: _description_, defaults to True
        :type fill_empty: bool, optional
        :return: _description_
        :rtype: dict
        """
        
        # compute histo
        hist, edges = np.histogram(a=X, bins=bins, 
                                   range=range, weights=weights,
                                   density=True)
        # compute fes (inverted histo)
        zeta = histo_to_fes(histo=hist, kbt=self.kbT,
                            zero_level=zero_level, 
                            fill_empty=fill_empty)
        # output
        fes_dict = dict(
            dim = 1,
            fes = zeta.T,
            grid = edges
        )
        return fes_dict
    

    def _fit2D(self,
               X: np.ndarray,
               Y: np.ndarray,
               bins: int,
               range: (float, float) =None,
               zero_level: Union[str, float] ='min',
               weights: np.ndarray =None,
               fill_empty=True) -> dict:
        """_summary_

        :param X: _description_
        :type X: np.ndarray
        :param Y: _description_
        :type Y: np.ndarray
        :param bins: _description_
        :type bins: int
        :param range: _description_, defaults to None
        :type range: float, float, optional
        :param zero_level: _description_, defaults to 'min'
        :type zero_level: Union[str, float], optional
        :param weights: _description_, defaults to None
        :type weights: np.ndarray, optional
        :param fill_empty: _description_, defaults to True
        :type fill_empty: bool, optional
        :return: _description_
        :rtype: dict
        """

        # compute histo
        hist, xedges, yedges = np.histogram2d(x=X, y=Y, bins=bins,
                                              range=range, weights=weights, 
                                              density=True)
        # compute fes
        zeta = histo_to_fes(histo=hist, kbt=self.kbT,
                            zero_level=zero_level,
                            fill_empty=fill_empty)
        XX, YY = mesh_grid_2d(X=xedges, Y=yedges, bins=bins)
        # output
        fes_dict = dict(
            dim = 2,
            fes = zeta.T,
            grid = (XX, YY)
        )
        return fes_dict

# -------------------------------------------------- #
# --- Plot

# import