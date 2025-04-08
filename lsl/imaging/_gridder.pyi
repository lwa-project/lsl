import numpy as np

from typing import Tuple

def WProjection(u: np.ndarray, v: np.ndarray, w: np.ndarray, data: np.ndarray, wgt: np.ndarray, uvSize: float=80.0, uvRes: float=0.5, wRes: float=0.1) -> Tuple[np.ndarray,np.ndarray,np.ndarray]: ...
