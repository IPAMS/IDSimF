import json

def sim_parameter_list(in_filename):
    with open(in_filename) as fh:
        
        data = json.load(fh)        
        result = ''
        
        for key in data:
            result += '  * - ``'+str(key)+"``\n    - \n    - \n"
        
    return result

def sim_parameter_def_list(in_filename):
    with open(in_filename) as fh:
        
        data = json.load(fh)        
        result = ''
        
        for key in data:
            result += '``'+str(key)+"`` : \n\n"
        
    return result