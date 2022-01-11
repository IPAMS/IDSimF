/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2020 - Physical and Theoretical Chemistry /
 Institute of Pure and Applied Mass Spectrometry
 of the University of Wuppertal, Germany

 IDSimF is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 IDSimF is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with IDSimF.  If not, see <https://www.gnu.org/licenses/>.
 ****************************/

#include "RS_ConfigFileParser.hpp"

using sMap = std::map<RS::Substance*,int>;
using sPair= sMap::value_type;



/**
 *  Parse Substances: parses the substance part of the configuration file
 * @param simConf a simulation configuration object to add the substances to
 * @param input the input string to parse (substances section of a config file in a string)
 * @return if the parsing was completed sucessfully
 */
bool RS::ConfigFileParser::parseSubstances(SimulationConfiguration *simConf, const std::string& input) const{

    std::istringstream iss(input);
    std::regex pattern(R"(([\w-]+)[ \t]*(\w+)[ \t]*([\deE._\+\-]*)[ \t]*([\deE._\-]*)[ \t]*([\deE._\-]*)[ \t]*([\deE._\-]*))");

    std::smatch matches;

    for (std::string line; std::getline(iss, line); )
    {
        if(std::regex_search(line, matches, pattern)) {
            std::string name = matches[1];
            std::string sType = matches[2];
            std::string optNumber1 = matches[3];
            std::string optNumber2 = matches[4];
            std::string optNumber3 = matches[5];
            std::string optNumber4 = matches[6];

            std::unique_ptr<RS::Substance> subst = std::unique_ptr<RS::Substance>(
                    new RS::Substance(name,sType));
            if (subst->type() == RS::Substance::substanceType::isotropic && optNumber1 != ""){
                subst->staticConcentration(std::strtof(optNumber1.c_str(), nullptr));
            }
            else if (subst->type() == RS::Substance::substanceType::isotropic){
                throw(RS::ConfigurationFileException("Isotropic substance without concentration value found"));
            }

            if (subst->type() == RS::Substance::substanceType::discrete &&
                optNumber1 != "" &&
                optNumber2 != "" &&
                optNumber3 != "" &&
                optNumber4 != "")
            {
                subst->mass(std::strtof(optNumber1.c_str(), nullptr));
                subst->charge(std::strtof(optNumber2.c_str(), nullptr));
                subst->mobility(std::strtof(optNumber3.c_str(), nullptr));
                subst->collisionDiameter(std::strtof(optNumber4.c_str(), nullptr));

            } else if (subst->type() == RS::Substance::substanceType::discrete) {
                throw (
                        RS::ConfigurationFileException(
                                "Discrete substance: Mass, charge, mobility or collision diameter value missing"));
            }

            simConf->addSubstance(subst);
        }
    }
    return true;
}

bool RS::ConfigFileParser::parseReactions(SimulationConfiguration *simConf, const std::string& input) const {
    std::istringstream iss(input);
    std::regex patternLine(R"(([\+\-_ \w]+)\s*=>\s*([\+\-_ \w]+)\s*\|\s*(\w+)\s*((?:\s*;\s*(?:[\deE._\+\-]+))+)\s*#?\s*(\w*))");
    //std::regex patternLine(R"(([\+\-_ \w]+)\s*=>\s*([\+\-_ \w]+)\s*;\s*(\w+)\s*(.*)\s*#?\s*(\w*))");
    std::regex patternParameters(R"(\s*;\s*([\deE._\+\-]+))");
    std::smatch matchesLine;
    std::smatch matchesParameters;

    for (std::string line; std::getline(iss, line); )
    {
        if(std::regex_search(line, matchesLine, patternLine)) {
            std::string eductFStr = matchesLine[1]; //the educt string of the formula part (left of the ";")
            std::string productFStr = matchesLine[2]; //the product string of the formula part (left of the ";")
            std::string typeStr = matchesLine[3]; //the reaction type keyword (after the first ";")
            std::string parametersStr = matchesLine[4]; //the parameters of the reaction (";" delimited list of numbers)
            std::string labelStr = matchesLine[5]; //the reaction label

            //strip whitespace from formula strings:
            eductFStr.erase(remove(eductFStr.begin(), eductFStr.end(), ' '), eductFStr.end());
            productFStr.erase(remove(productFStr.begin(), productFStr.end(), ' '), productFStr.end());

            //parse the educt side:
            std::map<Substance*,int> educts;
            std::vector<std::string> eductsStr = splitString(eductFStr,R"(\+)").first;
            for (auto edStr: eductsStr){
                std::pair<int,std::string> ed =  parseReactionPartnerString(edStr);
                RS::Substance* subst = simConf->substanceByName(ed.second);

                if (educts.count(subst) > 0){ //the substance is already present on the educt side
                    educts.at(subst) += ed.first; //the present substance multiplier is increased
                }else{
                    educts.insert(std::pair<Substance*,int>(subst,ed.first)); //new substance entry on the educt side
                }
            }

            //parse the product side:
            std::map<Substance*,int> products;
            std::vector<std::string> producsStr = splitString(productFStr,R"(\+)").first;
            for (auto proStr: producsStr){
                std::pair<int,std::string> pro =  parseReactionPartnerString(proStr);
                RS::Substance* subst = simConf->substanceByName(pro.second);

                if (products.count(subst) > 0){ //the substance is already present on the product side
                    products.at(subst) += pro.first; //the present substance multiplier is increased
                }else{
                    products.insert(std::pair<Substance*,int>(subst,pro.first)); //new substance entry on the product side
                }
            }

            // now parse reaction parameters (rate constant etc.) string
            // (details of the given parameters are specific to the individual reaction types)
            std::sregex_iterator iter(parametersStr.begin(), parametersStr.end(), patternParameters);
            std::sregex_iterator end; //defaults to an empty iterator, which equal to an iterator out of words

            std::vector<double> parsedParams;
            while(iter != end)
            {
                std::string matchResult = (*iter)[1];
                parsedParams.push_back(std::strtof( matchResult.c_str(), nullptr ));
                ++iter;
            }

            /*if(std::regex_search(parametersStr, matchesParameters, patternParameters)) {
                for (auto match: matchesParameters){
                    std::cout << match <<std::endl;
                }
            }*/
            /*double rateConstant = std::strtof(constantStr.c_str(), nullptr);
            std::unique_ptr<RS::AbstractReaction> reaction(
                    new RS::StaticReaction(educts,products,rateConstant,labelStr));
            simConf->addReaction(reaction);*/


            // now parse the reaction type string and prepare actual reaction
            std::unique_ptr<RS::AbstractReaction> reaction = nullptr;

            if (typeStr == "static"){
                reaction = std::make_unique<RS::StaticReaction>(
                        educts,
                        products,
                        parsedParams.at(0),
                        labelStr);
            }
            if (typeStr == "static_thermalizing"){
                reaction = std::make_unique<RS::StaticThermalizingReaction>(
                        educts,
                        products,
                        parsedParams.at(0),
                        labelStr);
            }
            else if(typeStr == "vanthoff"){
                reaction = std::make_unique<RS::VantHoffReaction>(
                        educts,
                        products,
                        parsedParams.at(0),
                        parsedParams.at(1),
                        parsedParams.at(2),
                        labelStr);
            }
            else if(typeStr == "simple_step"){
                reaction = std::make_unique<RS::SimpleCollisionStepReaction>(
                        educts,
                        products,
                        parsedParams.at(0),
                        labelStr);
                if ( !(reaction->discreteEducts()->size() == 1 && reaction->educts()->size() == 2) ){
                    throw (RS::ConfigurationFileException(
                            "Simple step reaction requires exactly one discrete and one isotropic educt"));
                }
            }
            else if(typeStr == "vanthoff_field"){
                //set mobility of the reaction as mobility of first discrete educt species:
                double mobility = 0.0;
                for (const auto& ed:educts){
                    if (ed.first->type() == RS::Substance::substanceType::discrete){
                        mobility = ed.first->mobility();
                    }
                }

                reaction = std::make_unique<RS::FieldDependentVantHoffReaction>(
                        educts,
                        products,
                        parsedParams.at(0),
                        parsedParams.at(1),
                        parsedParams.at(2),
                        mobility,
                        parsedParams.at(4),
                        parsedParams.at(3),
                        labelStr);
            }



            if (reaction == nullptr){
                throw (RS::ConfigurationFileException("Illegal reaction in configuration file"));
            }

            simConf->addReaction(reaction);
        }
    }
    return true;
}

std::pair<int,std::string> RS::ConfigFileParser::parseReactionPartnerString(const std::string& input) const{
    std::regex patternFull(R"((\d+)([\w_\-]+))");
    std::regex patternSymbol(R"([\w_\-]+)");
    std::smatch matches;


    int mult = 0;
    std::string speciesSymbol = "";

    if (std::isdigit(input[0])){
        if( std::regex_search(input, matches, patternFull) ) {
            mult = std::atoi(matches[1].str().c_str());
            speciesSymbol = matches[2];
        }
    }
    else{
        if( std::regex_search(input, matches, patternSymbol) ) {
            mult = 1;
            speciesSymbol = matches[0];
        }
    }
    return std::pair<int,std::string>(mult,speciesSymbol);
}


std::pair<std::vector<std::string>,std::vector<std::string>> RS::ConfigFileParser::splitString(
        std::string str, const std::string& patterntxt) const{

    std::vector<std::string> outMatches = std::vector<std::string>();
    std::vector<std::string> outResults = std::vector<std::string>();

    std::regex pattern(patterntxt);
    std::smatch m;
    //std::regex_search (str, m, pattern);

    while (std::regex_search (str, m, pattern)) {
        if (m.position(0)>0) {
            outResults.push_back(str.substr(0, static_cast<std::size_t>(m.position(0))));
        }
        outMatches.push_back(m[0].str());
        str = m.suffix().str();
    }
    if (str.length() > 0){
        outResults.push_back(str);
    }
    return std::make_pair(outResults,outMatches);
}

std::unique_ptr<RS::SimulationConfiguration> RS::ConfigFileParser::parseFile(const std::string& filename) const{
    //open stream:
    std::ifstream in;
    in.open(filename);

    if (in.good()) {
        //parse: read line by line, tokenize the lines
        std::stringstream buffer;
        buffer << in.rdbuf();
        std::string confStr = buffer.str();

        return (parseText(confStr));
    } else {
        throw (RS::ConfigurationFileException("file not found"));
    }
}

std::unique_ptr<RS::SimulationConfiguration> RS::ConfigFileParser::parseText(const std::string& confStr) const {

    std::unique_ptr<RS::SimulationConfiguration> simConf =
            std::unique_ptr<RS::SimulationConfiguration>(new RS::SimulationConfiguration());

    //split into the substances part and the reactions part
    std::string pattern = R"(\[\w+\])";
    std::vector<std::string> parts (splitString(confStr,pattern).first);

    if (parts.size() < 2){
        throw(RS::ConfigurationFileException("Configuration file parts missing"));
    }

    //parse substances
    parseSubstances(simConf.get(),parts[0]);

    //parse reactions
    parseReactions(simConf.get(),parts[1]);

    return(simConf);
}

/**
 *
 * @return
 */
std::unique_ptr<RS::SimulationConfiguration> RS::ConfigFileParser::getTestConfigWaterClusters() const {
    std::unique_ptr<RS::SimulationConfiguration> simConf =
            std::unique_ptr<RS::SimulationConfiguration>(new RS::SimulationConfiguration());

    //std::vector<RS::Substance> substances;

    std::unique_ptr<RS::Substance> HydrogenIons = std::unique_ptr<RS::Substance>(
            new RS::Substance("[H3]+",RS::Substance::substanceType::discrete));

    std::unique_ptr<RS::Substance> HydroxyIons = std::unique_ptr<RS::Substance>(
            new RS::Substance("[OH]-",RS::Substance::substanceType::discrete));

    std::unique_ptr<RS::Substance> Cluster_1 = std::unique_ptr<RS::Substance>(
            new RS::Substance("[H3O]+",RS::Substance::substanceType::discrete));

    std::unique_ptr<RS::Substance> Cluster_2 = std::unique_ptr<RS::Substance>(
            new RS::Substance("[H2O+H3O]+",RS::Substance::substanceType::discrete));
    ;
    std::unique_ptr<RS::Substance> Nitrogen = std::unique_ptr<RS::Substance>(
            new RS::Substance("N2",RS::Substance::substanceType::isotropic));

    std::unique_ptr<RS::Substance> Water = std::unique_ptr<RS::Substance>(
            new RS::Substance("H2O",RS::Substance::substanceType::isotropic));


    sMap e_Cl1_Cl2;
    sMap p_Cl1_Cl2;
    e_Cl1_Cl2.insert(sPair(Cluster_1.get(),1));
    e_Cl1_Cl2.insert(sPair(Water.get(),1));
    e_Cl1_Cl2.insert(sPair(Nitrogen.get(),1));

    p_Cl1_Cl2.insert(sPair(Cluster_2.get(),1));
    p_Cl1_Cl2.insert(sPair(Nitrogen.get(),1));

    std::unique_ptr<RS::AbstractReaction> r_Cl1_Cl2( new RS::StaticReaction(e_Cl1_Cl2,p_Cl1_Cl2,6.98e-29,"Cl1->Cl2"));
    simConf->addReaction(r_Cl1_Cl2);

    sMap e_recomb;
    sMap p_recomb;
    e_recomb.insert(sPair(HydroxyIons.get(),1));
    e_recomb.insert(sPair(Cluster_1.get(),1));
    e_recomb.insert(sPair(Nitrogen.get(),1));

    p_recomb.insert(sPair(Water.get(),2));
    p_recomb.insert(sPair(Nitrogen.get(),1));
    std::unique_ptr<RS::AbstractReaction> r_recomb(new RS::StaticReaction(e_recomb,p_recomb,6.98e-29,"Recombination"));
    simConf->addReaction(r_recomb);

    RS::SimulationConfiguration config;

    //substances are added at last because the add / move operation invalidates the pointer to the
    //substances in this scope:
    simConf->addSubstance(HydrogenIons);
    simConf->addSubstance(HydroxyIons);
    simConf->addSubstance(Cluster_1);
    simConf->addSubstance(Cluster_2);
    simConf->addSubstance(Nitrogen);
    simConf->addSubstance(Water);

    return (simConf);
}

/**
 *
 * @return
 */
std::unique_ptr<RS::SimulationConfiguration> RS::ConfigFileParser::getTestConfigSimple() const{
    std::unique_ptr<RS::SimulationConfiguration> simConf =
            std::unique_ptr<RS::SimulationConfiguration>(new RS::SimulationConfiguration());

    //std::vector<RS::Substance> substances;

    std::unique_ptr<RS::Substance> s_A= std::unique_ptr<RS::Substance>(
            new RS::Substance("A",RS::Substance::substanceType::discrete));

    std::unique_ptr<RS::Substance> s_B= std::unique_ptr<RS::Substance>(
            new RS::Substance("B",RS::Substance::substanceType::isotropic));

    std::unique_ptr<RS::Substance> s_C= std::unique_ptr<RS::Substance>(
            new RS::Substance("C",RS::Substance::substanceType::discrete));

    s_A->mobility(0.1);
    s_B->mobility(0.2);
    s_C->mobility(0.3);

    s_A->collisionDiameter(1e-5);
    s_B->collisionDiameter(2e-5);
    s_C->collisionDiameter(3e-5);

    s_A->charge(1);
    s_B->charge(2);
    s_C->charge(3);

    s_A->mass(100);
    s_A->mass(200);
    s_A->mass(300);


    sMap e_R1;
    sMap p_R1;
    e_R1.insert(sPair(s_A.get(),1));
    e_R1.insert(sPair(s_B.get(),1));

    p_R1.insert(sPair(s_C.get(),1));
    p_R1.insert(sPair(s_B.get(),1));

    std::unique_ptr<RS::AbstractReaction> r_R1( new RS::StaticReaction(e_R1,p_R1,1.0,"Reaction 1"));
    simConf->addReaction(r_R1);

    simConf->addSubstance(s_A);
    simConf->addSubstance(s_B);
    simConf->addSubstance(s_C);

    return (simConf);
}