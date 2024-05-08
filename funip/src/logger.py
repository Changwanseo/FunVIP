from datetime import datetime
import shelve


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
