from datetime import datetime
import os
import shelve
import subprocess
import logging


class CustomFormatter(logging.Formatter):
    # Initialize ansicode on windows
    if os.name == "nt":  # Only if we are running on Windows
        from ctypes import windll

        k = windll.kernel32
        k.SetConsoleMode(k.GetStdHandle(-11), 7)

    # ANSI \x1b
    # UNICODE \u001b

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

    # I don't know why, but in some environment, ANSI color works only after first subprocess.call was done
    _ = subprocess.call("", shell=True)

    for info in list_info:
        logging.info(info)

    for warning in list_warning:
        logging.warning(warning)

    for error in list_error:
        logging.error(error)
