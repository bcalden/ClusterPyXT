import numpy as np

class Image:
    """
    The `Image` class in CPXT encodes attributes such as the filename, header,
    and the actual underlying data. This allows for working with the image
    data and attributes in a more 'pythonic' way. 
    Attributes
    ----------
    filename : str
        The filename should be the full path as to where it is stored on disk
    filepath : Path
    header: FITSHeader
    data: np.array
    
    Methods
    -------
    <method> : <return type>
    """
    
    def __init__(self):
        pass