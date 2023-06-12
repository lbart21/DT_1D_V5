"""
Function:
Author: Luke Bartholomew
Edits:
"""
class EmptyMeshObject():
    """
    Initialises mesh object attributes for manipulation later.
    """
    __slots__ = ["cell_array", "interface_array", "map_cell_id_to_west_interface_idx", \
                    "map_cell_id_to_east_interface_idx", "map_interface_id_to_west_cell_idx", \
                    "map_interface_id_to_east_cell_idx", "component_labels", "boundary_conditions", \
                    "boundary_interface_ids", "cell_idx_to_track", "interface_idx_to_track"]
    def __init__(self):
        self.cell_array = None
        self.interface_array = None
        self.map_cell_id_to_west_interface_idx = None
        self.map_cell_id_to_east_interface_idx = None
        self.map_interface_id_to_west_cell_idx = None
        self.map_interface_id_to_east_cell_idx = None
        self.boundary_conditions = None
        self.boundary_interface_ids = None
        self.component_labels = None
        self.cell_idx_to_track = None
        self.interface_idx_to_track = None
