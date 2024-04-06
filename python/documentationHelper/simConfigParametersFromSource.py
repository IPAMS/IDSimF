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

def parse_file(fname): 

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
    

def extract_parameters_from_source(sourcepath, print_mode='list'):

    result = []
    for root, dirs, files in os.walk(sourcepath):
        for name in files:
            if name.endswith((".cpp", ".hpp")):
                result.append( (name, parse_file(os.path.join(root,name))))

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
            print('\n')


def main():
    #extract_parameters_from_source('../../applications', print_mode='list')
    extract_parameters_from_source('../../applications', print_mode='rst_with_notes')

if __name__ == "__main__":
    main()

