def generate_average_from_paraview_slice():
    """
    Pass in cell data and point data files for slice
    Process files into objects
    Sort point data object by pos.y because they aren't always in the right order.
    There are now n cells and n+1 points.
    Loop over all cells:
        average = sum(phi * pi * (y(i+1)^2 - y(i)^2)) / sum(pi * (y(i+1)^2 - y(i)^2))
    y(i) is the slice intersection location for the cell with the lowest y value
    """
    # TODO include flag to do average average or mass flux average
    # mass flux average = sum(phi * rho * vel_x * pi * (y2^2 - y1^2)) / sum(rho * vel_x * pi*(y2^2 - y1^2))
    
    pass