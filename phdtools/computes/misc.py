# -------------------------------------------------- #
# Computes - miscellaneous stuffs
#
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

from datetime import date

# --- get todays date string
def todayDate(date_format="%d%b%Y"):
    """
    Gets you todays date in according to the format choosen
    Default: DayMonthYear (ex. 13Jan1992)
    """
    today = date.today()
    return today.strftime(date_format)