import numpy as np

def time_correlation(
        time: np.ndarray, quant: np.ndarray, 
        window_shift : int = 10, time_split : float = 0.5) -> np.ndarray:
    """
    Compute the time correlation of a quantity.

    Parameters
    ----------
    time : np.ndarray
        The time array.
    quant : np.ndarray
        The quantity array.
    window_shift : int, optional
        The shift of the window. The default is 10.
    time_split : float, optional
        The time split. The default is 0.5.

    Returns
    -------
    np.ndarray
        The time correlation array.

    """

    raise NotImplementedError

    if time.shape != quant.shape:
        raise ValueError('time and quant must have the same shape')
    
    time_window_lenght = int(time_split * time.shape[0])
    time_window = np.arange(time_window_lenght)

    time_correlation = np.zeros(time_window_lenght)
    counts = np.zeros(time_window_lenght)

    init_time_index = 0
    while (init_time_index + time_window_lenght) < time.shape[0]):
        time_correlation += np.sum(
            quant[init_time_index + time_window] * quant[init_time_index: init_time_index + time_window_lenght], axis=0)
        counts += time_window_lenght
        init_time_index += window_shift