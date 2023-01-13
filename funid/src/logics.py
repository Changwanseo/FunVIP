import logging


def isnewicklegal(string):
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
    list_column: list, column: str, table_name: str, check_none: bool = True
) -> bool:
    """
    if given column is unique in list_column -> return True
    elif give column is not in list_column -> return False (Error when check_none is True)
    if more than one column found -> raise Error
    """

    columns = [x.lower().strip() for x in list_column]
    cnt_column = columns.count(column)

    if cnt_column == 1:
        return True
    elif check_none is False and cnt_column == 0:
        return False
    elif cnt_column > 1:
        logging.error(f'More than 1 column of "{column}" found in {table_name}')
        raise Exception
    elif check_none is True and cnt_column == 0:
        logging.error(f'Column "{column}" is mandatory, but not found in {table_name}')
        raise Exception
    else:
        print(columns)
        print(column)
        print(cnt_column)
        raise Exception
