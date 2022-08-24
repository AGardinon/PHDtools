# -------------------------------------------------- #
# Computes - miscellaneous stuffs
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

# --- custom exception
class WrongConfigFormat(Exception):
    pass

# --- get todays date string
from datetime import date

def todayDate(date_format="%d%b%Y"):
    """
    Gets you todays date in according to the format choosen
    Default: DayMonthYear (ex. 13Jan1992)
    """
    today = date.today()
    return today.strftime(date_format)

# --- parse dict file .json .toml
import json
import toml
from types import SimpleNamespace

def parse_config(filename, toNameSpace=False):
    with open(filename, 'r') as f:
        if filename.endswith('.json'):
            _config = json.load(f)
        elif filename.endswith('.toml'):
            _config = toml.load(f)
        else:
            raise WrongConfigFormat(
                f"Format of the '{filename}' is not supported. "
                "Available formats: .json, .toml"
            )
        if toNameSpace:
            return SimpleNamespace(**_config)
        else:
            return _config