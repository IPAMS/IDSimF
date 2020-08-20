//
// Created by dystopic on 15.01.18.
//

#ifndef IONSIMULATION_CPP_PARAMETERPARSING_HPP
#define IONSIMULATION_CPP_PARAMETERPARSING_HPP

#define doubleConfParameter(jsonName,varName)    \
if (confRoot.isMember(jsonName)==true){        \
    (varName) = confRoot.get(jsonName,0).asDouble(); \
    std::cout << (jsonName)<<":"<<(varName)<<std::endl;    \
}else{                                           \
    std::cout << "missing configuration value: "<< (jsonName)<<std::endl; \
    return(0); \
}

#define intConfParameter(jsonName,varName)    \
if (confRoot.isMember(jsonName)==true){        \
    (varName) = confRoot.get(jsonName,0).asInt(); \
    std::cout << (jsonName)<<":"<<(varName)<<std::endl;    \
}else{                                           \
    std::cout << "missing configuration value: "<< (jsonName)<<std::endl; \
    return(0); \
}

#define strConfParameter(jsonName,varName)    \
if (confRoot.isMember(jsonName)==true){        \
    varName = confRoot.get(jsonName,0); \
    std::cout << jsonName<<":"<<varName<<std::endl;    \
}else{                                           \
    std::cout << "missing configuration value: "<< jsonName<<std::endl; \
    return(0); \
}

#endif //IONSIMULATION_CPP_PARAMETERPARSING_HPP
