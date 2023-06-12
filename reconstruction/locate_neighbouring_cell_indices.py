"""
Function:
Author: Luke Bartholomew
Edits:
"""
def find_idx_of_cell_recursively(interface_id, recursion_depth, direction, \
                                map_interface_id_to_east_cell_idx, \
                                map_cell_id_to_east_interface_idx, \
                                map_interface_id_to_west_cell_idx, \
                                map_cell_id_to_west_interface_idx):
        """ Given interface ID, direction you want to look in and how far you want to look, return the index of the cell that you end at
        Eg. recursion_depth = 1, direction = "East" will give you the index of the East cell of the interface [ | | | | | !0| | | | | | ]
            recursion_depth = 3, direction = "West" will give you the index of the West cell of the West cell of the West cell of the interface [ | | | |0| | ! | | | | | ]
        """
        current_interface_id = interface_id
        current_recursion_depth = 0
        hit_bad_interface = False
        
        while current_recursion_depth < recursion_depth:
            if direction == "East":
                new_cell_idx = map_interface_id_to_east_cell_idx[current_interface_id]
                if new_cell_idx is None:
                    print("No cell in the East direction")
                    hit_bad_interface = True
                    return new_cell_idx, hit_bad_interface
                else:
                    current_recursion_depth += 1
                    if current_recursion_depth != recursion_depth:
                        if len(map_cell_id_to_east_interface_idx[new_cell_idx]) > 1:
                            print("Reconstuction entering multi-interface cell boundary is not implemented yet.")
                            hit_bad_interface = True
                            return new_cell_idx, hit_bad_interface
                        else:
                            current_interface_id = map_cell_id_to_east_interface_idx[new_cell_idx][0]
            
            if direction == "West":
                new_cell_idx = map_interface_id_to_west_cell_idx[current_interface_id]
                if new_cell_idx is None:
                    print("No cell in the West direction")
                    hit_bad_interface = True
                    return new_cell_idx, hit_bad_interface
                else:
                    current_recursion_depth += 1
                    if current_recursion_depth != recursion_depth:
                        if len(map_cell_id_to_west_interface_idx[new_cell_idx]) > 1:
                            print("Reconstuction entering multi-interface cell boundary is not implemented yet.")
                            hit_bad_interface = True
                            return new_cell_idx, hit_bad_interface
                        else:
                            current_interface_id = map_cell_id_to_west_interface_idx[new_cell_idx][0]

        return new_cell_idx, hit_bad_interface
