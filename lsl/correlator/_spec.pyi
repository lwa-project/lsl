import numpy as np

from typing import Callable, Optional

def FPSD(signals: np.ndarray, LFFT: int=64, overlap: int=1, clip_level: int=0, window: Optional[int]=None) -> np.ndarray: ...

def PFBPSD(signals: np.ndarray, LFFT: int=64, overlap: int=1, clip_level: int=0, window: Optional[int]=None) -> np.ndarray: ...
