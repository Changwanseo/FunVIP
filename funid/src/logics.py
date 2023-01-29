import logging
import re
from matplotlib import colors as mcolors  # for color debugging

# Move all simple logics function deciding True all False here


def isnan(value) -> bool:
    if type(value) is float or type(value) is np.float64:
        if np.isnan(value):
            return True
    return False


def isvalidcolor(color: str) -> bool:

    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    # Check if color is named colors
    if color in colors.keys():
        return True
    # Check if color is available hex color
    elif re.fullmatch(r"^#(?:[0-9a-fA-F]{3}){1,2}$", color):
        return True
    else:
        return False


def isnewicklegal(string: str) -> bool:
    """
    if string is newick legal -> return True
    if string is newick illegal -> return False
    """
    NEWICK_ILLEGAL = (
        "(",
        '"',
        "[",
        ":",
        ";",
        "/",
        "[",
        "]",
        "{",
        "}",
        "(",
        ")",
        ",",
        "]",
        "+",
        '"',
        ")",
        " ",
    )

    if any(x in string for x in NEWICK_ILLEGAL):
        return True
    else:
        return False


def isuniquecolumn(
    list_column: list, column: tuple, table_name: str, check_none: bool = True
) -> bool:

    """
    if given column is unique in list_column -> return True
    elif give column is not in list_column -> return False (Error when check_none is True)
    if more than one column found -> raise Error
    """

    # Lower and make safe column names
    columns = [x.lower().strip() for x in list_column]

    # Count given columns
    # check ambiguities
    cnt_column = 0
    selected_column = None

    for col in column:
        cnt_column += columns.count(col)
        if columns.count(col) != 0:
            selected_column = col

    # If only one of the given column exists
    if cnt_column == 1:
        return selected_column
    # If column does not exists and column is not mandatory
    elif check_none is False and cnt_column == 0:
        return False
    # If more than one column exists
    elif cnt_column > 1:
        logging.error(
            f'More than 1 column of "{" or ".join(column)}" found in {table_name}'
        )
        raise Exception
    # If column is mandatory and does not exists
    elif check_none is True and cnt_column == 0:
        logging.error(f'Column "{column}" is mandatory, but not found in {table_name}')
        raise Exception
    else:
        logging.error(f"Unknown error occured while checking table column names")
        logging.debug(columns)
        logging.debug(column)
        logging.debug(cnt_column)
        raise Exception
