# -------------------------------------------------- #
# QUIP tools - Handler for quipd descriptors and 
#              related objects
#
# AUTHOR: Andrea Gardin
# -------------------------------------------------- #

import copy

# -------------------------------------------------- #
# --- QUIP descriptor tools

QUIPTYPES = ['soap', 'soap_turbo']

def quipSOAP_string(params_dict, type):
    if type == 'soap':
        return quipSOAP_to_str(params_dict)
    elif type == 'soap_turbo':
        return quipSOAPTURBO_to_str(params_dict)
    else:
        raise NameError(f"Type not implemented, types available: {QUIPTYPES}.")

def quipSOAP_to_str(params_dict):
    dscr_str = 'soap '
    for param,value in params_dict.items():
        if isinstance(value, list):
            entry = param+'={'+' '.join(map(str,value))+'} '
        else:
            entry = param+'='+str(value)+' '
        dscr_str += entry
    return dscr_str

def quipSOAPTURBO_to_str(params_dict):
    dscr_str = ''
    for Z in params_dict["Zs"]:
        params_dict["central_index"] = params_dict["species_Z"].index(Z)+1
        qs = _params2qs(params_dict, params_dict["multi"])
        dscr_str += qs+' '
    return dscr_str

def _params2qs(orig_params, multi):
    """convert param dicts to quippy string"""
    params = copy.deepcopy(orig_params)
    params["n_species"] = len(params["species_Z"])
    Zstr = "{"
    for Z in params["species_Z"]:
        Zstr += str(Z) + " "
    Zstr = Zstr[:-1]+"}"
    params["species_Z"] = Zstr
    for key, value in multi.items():
        s = "{"
        for i in range(0, params["n_species"]):
            s += str(value) + " "
        params[key] = s[:-1] + "}"
    s = 'soap_turbo '
    for key, value in params.items():
        if key != "Zs":
            s += key + "=" + str(value) + " "
    return s


# -------------------------------------------------- #
# --- DScribe descriptor tools

# empty