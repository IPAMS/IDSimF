import os
import re

app_utils_methods = (
    'runTrajectoryIntegration',
    'isIonCloudDefinitionPresent',
    'readIonDefinitionFromIonCloudFile',
    'getStartZoneFromIonDefinition',
    'setIonsKineticEnergy',
    'readIonDefinition'
)

def parse_cpp_file(fname): 

    pattern_param = re.compile(r""".*simConf->(?P<method>.*)\(\"(?P<key>.*)\"\).*""")
    pattern_app_util_method = re.compile(r""".*AppUtils::(?P<method>.*)\(.*""")


    result = []
    with open(fname) as fh:
        for line in fh:
            match = pattern_param.match(line)
            if match:
                result.append(
                    (match.group('key'), match.group('method'))
                )
            match_app_util = pattern_app_util_method.match(line)
            if match_app_util:
                method = match_app_util.group('method')
                if method in app_utils_methods:
                    result.append(
                        (method,)
                    )

    return(result)


def parse_rst_file(fname): 

    pattern_param = re.compile(r"""\s*``(?P<key>.*)``\s?:\s?(?P<type>.*)""")


    result = []
    with open(fname) as fh:
        for line in fh:
            match = pattern_param.match(line)
            if match:
                result.append(
                    (match.group('key'), match.group('type'))
                )
    return(result)


replacements = {
    'isParameter'     : '(is parameter?)',
    'boolParameter'   : 'boolean',
    'intParameter'    : 'integer',
    'intVectorParameter'    : 'vector of integers',
    'unsignedIntParameter' : 'unsigned int',
    'unsignedIntVectorParameter' : 'vector of unsigned int',
    'double3dBox'          : 'box vectors (double)',
    'doubleParameter' : 'double',
    'doubleVectorParameter' : 'vector of doubles' ,
    'stringParameter' : 'string',
    'stringVectorParameter' : 'vector of strings' 
}

def print_param_rst(param, print_mode):
    if len(param) == 2:
        if param[1] in replacements:
            type_str = replacements[param[1]]
        else:
            type_str = f'undef -> {param[1]}'

        print(f' ``{param[0]}`` : {type_str}')
        if print_mode == 'rst_with_notes':
            print('    (Fixme: Explanation)\n')
    else:
        print(f' **** -> AppUtil method present: {param[0]}\n')    

def extract_parameters_from_source(sourcepath, print_mode='list'):

    result = []
    for root, dirs, files in os.walk(sourcepath):
        for name in files:
            if name.endswith((".cpp", ".hpp")):
                result.append( (name, parse_cpp_file(os.path.join(root,name))))

    for res in result:
        if len(res[1]) > 0:
            if print_mode == 'list':
                print(res[0])
                for param in res[1]:
                    if len(res[1] == 2):
                        print(f'     ->  {param[0]:<30}        {param[1]}')

            elif print_mode == 'rst' or  print_mode == 'rst_with_notes':
                print(res[0])
                for param in res[1]:
                    print_param_rst(param, print_mode)

                    
            print('\n')


def compare_params(rst_file, source_file):
    rst_params_list = parse_rst_file(rst_file)
    rst_params = {rp[0]: rp[1] for rp in rst_params_list}

    cpp_params = parse_cpp_file(source_file)

    for cpp_p in cpp_params:
        if cpp_p[0] not in rst_params:
            #print(f'  -> not found in doc {cpp_p[0]}')

            print_param_rst(cpp_p, 'rst_with_notes')



def main():
    #extract_parameters_from_source('../../applications', print_mode='list')
    extract_parameters_from_source('../../applications', print_mode='rst_with_notes')
    
    #compare_params(
    #    '../../documentation/usersguide/applications/QITSim.rst',
    #    '../../applications/ionTraps/QITSIm/QITSim.cpp')

if __name__ == "__main__":
    main()

