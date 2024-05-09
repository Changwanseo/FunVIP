from datetime import datetime
import shelve

# I think this is deprecated
"""
def initialize(log):
    global log_file
    log_file = log


# Logging Functions
def Time_now():
    now = datetime.now()
    Date_time = "%s-%s-%s %s:%s:%s " % (
        now.year,
        now.month,
        now.day,
        now.hour,
        now.minute,
        now.second,
    )
    return Date_time


def Mes(text):
    text = str(text)
    text = Time_now() + text

    try:
        log_file.write(text + "\n")
    except:
        pass

    print(text)
"""


# should be deleted
"""
def Save(path):
    save = shelve.open(path.save, "n")
    for key in dir():
        try:
            save[key] = globals()[key]
        except TypeError:
            print(f"ERROR shelving: {key}")
    save.close()
"""


import logging


class CustomFormatter(logging.Formatter):
    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    white = "\x1b[38;97m"
    reset = "\x1b[0m"
    fmt = "%(asctime)s [%(levelname)s] %(message)s"
    fmt_err = "%(asctime)s [%(levelname)s] %(message)s (%(filename)s:%(lineno)d)"

    FORMATS = {
        logging.DEBUG: grey + fmt + reset,
        logging.INFO: white + fmt + reset,
        logging.WARNING: yellow + fmt_err + reset,
        logging.ERROR: bold_red + fmt_err + reset,
        logging.CRITICAL: bold_red + fmt_err + reset,
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, datefmt="%m/%d/%Y %I:%M:%S %p")
        return formatter.format(record)


def setup_logging(list_info, list_warning, list_error, path, opt, tool):
    # setup logging
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        encoding="utf-8",
        level=tool.get_level(opt.verbose),
        handlers=[logging.FileHandler(path.log), logging.StreamHandler()],
    )

    # Set the custom formatter
    formatter = CustomFormatter()
    for handler in logging.root.handlers:
        handler.setFormatter(formatter)

    # Delayed logging for option parsing
    for info in list_info:
        logging.info(info)

    for warning in list_warning:
        logging.warning(warning)

    for error in list_error:
        logging.error(error)
