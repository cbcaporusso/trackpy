import numpy as np

# def time_correlation(
#         time: np.ndarray, quant: np.ndarray, 
#         window_shift : int = 10, time_split : float = 0.5) -> np.ndarray:
#     """
#     Compute the time correlation of a quantity.

#     Parameters
#     ----------
#     time : np.ndarray
#         The time array.
#     quant : np.ndarray
#         The quantity array.
#     window_shift : int, optional
#         The shift of the window. The default is 10.
#     time_split : float, optional
#         The time split. The default is 0.5.

#     Returns
#     -------
#     np.ndarray
#         The time correlation array.

#     """

#     raise NotImplementedError

#     if time.shape != quant.shape:
#         raise ValueError('time and quant must have the same shape')
    
#     time_window_lenght = int(time_split * time.shape[0])
#     time_window = np.arange(time_window_lenght)

#     time_correlation = np.zeros(time_window_lenght)
#     counts = np.zeros(time_window_lenght)

#     init_time_index = 0
#     while (init_time_index + time_window_lenght) < time.shape[0]):
#         time_correlation += np.sum(
#             quant[init_time_index + time_window] * quant[init_time_index: init_time_index + time_window_lenght], axis=0)
#         counts += time_window_lenght
#         init_time_index += window_shift

# import numpy as np

def cg_field_to_cells(pos: np.ndarray, cell_size: float, quantity):
    """
    Compute the coarse-grained field from a field defined on a set of particles.

    Parameters
    ----------
    pos : np.ndarray
        The position array.
    cell_size : float
        The cell size.
    quantity : np.ndarray
        The quantity array.

    Returns
    -------
    np.ndarray
        The coarse-grained field.

    """
      
    box_size = np.amax(pos,axis=0)
    
    cells = pos/cell_size
    cells = cells.astype(int)
    num_part = len(cells[:,0])
    
    num_cells = np.amax(cells,axis=0)+1

    coarsed_quantity = np.zeros(num_cells)
    nn = np.zeros(num_cells)

    for i in range(num_part):
        coarsed_quantity[tuple(cells[i])] += quantity[i] 
        nn[tuple(cells[i])] += 1 


    xpos = []
    ypos = []
    q_array = []

    trim=2

    for xx in range(trim,num_cells[0]-trim):
        for yy in range(trim,num_cells[1]-trim):
            xpos.append( cell_size * (.5 + xx))
            ypos.append( cell_size * (.5 + yy))
            
            q_array.append(coarsed_quantity[xx,yy]/nn[xx,yy]) if nn[xx,yy]>0 else q_array.append(0)

    pos = np.column_stack((np.array(xpos),np.array(ypos)))

    return pos, np.array(q_array)

