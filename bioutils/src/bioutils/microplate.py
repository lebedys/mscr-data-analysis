"""
Utilities for handling data from protocols with microplates.
"""


# 96-well plate: 8 rows (A-H) x 12 columns (1-12)
ROWS_96 = 8
COLS_96 = 12

# 384-well plate: 16 rows (A-P) x 24 columns (1-24)
ROWS_384 = 16
COLS_384 = 24


def well_coords_to_index_96(row: int, col: int) -> int:
    """
    (0,0)  -> 0
    (0,1)  -> 1
    ...
    (1,0)  -> 12
    ...
    (7,11) -> 95
    """
    return row * COLS_96 + col


def well_index_to_coords_96(index: int) -> tuple[int, int]:
    """
    0  -> (0,0)
    1  -> (0,1)
    ...
    12 -> (1,0)
    ...
    95 -> (7,11)
    """
    return (index // COLS_96, index % COLS_96)


def well_coords_to_label_96(row: int, col: int) -> str:
    """
    (0,0)  -> "A1"
    (0,1)  -> "A2"
    ...
    (1,0)  -> "B1"
    ...
    (7,11) -> "H12"
    """
    return f"{chr(ord('A') + row)}{col + 1}"


def well_label_to_coords_96(label: str) -> tuple[int, int]:
    """
    "A1"  -> (0,0)
    "A2"  -> (0,1)
    ...
    "B1"  -> (1,0)
    ...
    "H12" -> (7,11)
    """
    row = ord(label[0].upper()) - ord('A')
    col = int(label[1:]) - 1
    return row, col


def well_index_to_label_96(index: int) -> str:
    """
    0  -> "A1"
    1  -> "A2"
    ...
    12 -> "B1"
    ...
    95 -> "H12"
    """
    row, col = well_index_to_coords_96(index)
    return well_coords_to_label_96(row, col)


def well_label_to_index_96(label: str) -> int:
    """
    "A1"  -> 0
    "A2"  -> 1
    ...
    "B1"  -> 12
    ...
    "H12" -> 95
    """
    row, col = well_label_to_coords_96(label)
    return well_coords_to_index_96(row, col)


def well_coords_to_index_384(row: int, col: int) -> int:
    """
    (0,0)   -> 0
    (0,1)   -> 1
    ...
    (1,0)   -> 24
    ...
    (15,23) -> 383
    """
    return row * COLS_384 + col


def well_index_to_coords_384(index: int) -> tuple[int, int]:
    """
    0   -> (0,0)
    1   -> (0,1)
    ...
    24  -> (1,0)
    ...
    383 -> (15,23)
    """
    return index // COLS_384, index % COLS_384


def well_coords_to_label_384(row: int, col: int) -> str:
    """
    (0,0)   -> "A1"
    (0,1)   -> "A2"
    ...
    (1,0)   -> "B1"
    ...
    (12,4)  -> "M5"
    ...
    (15,23) -> "P24"
    """
    return f"{chr(ord('A') + row)}{col + 1}"


def well_label_to_coords_384(label: str) -> tuple[int, int]:
    """
    "A1"  -> (0,0)
    "A2"  -> (0,1)
    ...
    "B1"  -> (1,0)
    ...
    "M5"  -> (12,4)
    ...
    "P24" -> (15,23)
    """
    row = ord(label[0].upper()) - ord('A')
    col = int(label[1:]) - 1
    return row, col


def well_index_to_label_384(index: int) -> str:
    """
    0   -> "A1"
    1   -> "A2"
    ...
    24  -> "B1"
    ...
    383 -> "P24"
    """
    row, col = well_index_to_coords_384(index)
    return well_coords_to_label_384(row, col)


def well_label_to_index_384(label: str) -> int:
    """
    "A1"  -> 0
    "A2"  -> 1
    ...
    "B1"  -> 24
    ...
    "P24" -> 383
    """
    row, col = well_label_to_coords_384(label)
    return well_coords_to_index_384(row, col)


def get_adjacent_wells(row: int, col: int, width: int) -> list[tuple[int, int]]:
    """
    (0, 0, 3) -> [(0,0), (0,1), (0,2)]
    (7, 16, 4) -> [(7,16), (7,17), (7,18), (7,19)]
    """
    return [(row, col + i) for i in range(width)]


def get_circuit_wells(circuit: int) -> list[tuple[int, int]]:
    """
    Circuit 1-92 -> its four wells, one per inducer condition (None, aTc, IPTG, Both).

    1  -> [(0,0),   (0,1),  (0,2),  (0,3)]    A1-A4
    2  -> [(1,0),   (1,1),  (1,2),  (1,3)]    B1-B4
    46 -> [(7,16),  (7,17), (7,18), (7,19)]   H17-H20
    92 -> [(15,12), (15,13), (15,14), (15,15)] P13-P16
    """
    index = circuit - 1
    row = 2 * (index // 12) + (index % 2)
    col = 4 * ((index % 12) // 2)
    return get_adjacent_wells(row, col, 4)
