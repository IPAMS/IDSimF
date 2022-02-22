/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2021 - Physical and Theoretical Chemistry /
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

 ------------
 appUtils_signalHandler.hpp

 Basic handling of OS signals, primarily SIGTERM

 ****************************/
#ifndef IDSIMF_APPUTILS_SIGNALHANDLER_HPP
#define IDSIMF_APPUTILS_SIGNALHANDLER_HPP

#include "Integration_abstractTimeIntegrator.hpp"
#include <stdexcept>

namespace AppUtils{
    class SignalException : public std::runtime_error {
    public:
        SignalException(const std::string& _message)
                :std::runtime_error(_message) { }
    };

    /**
     * Simple management class which provides basic handling of POSIX signals
     */
    class SignalHandler {

    public:
        static void setReceiver(Integration::AbstractTimeIntegrator& receiver);
        static void registerSignalHandler();
        static void sendTerminateToReceiver();
        static bool isTerminationSignaled();

    private:
        inline static Integration::AbstractTimeIntegrator* receiver_;
        inline static bool terminationSignaled_ = false;

    };
}

#endif //IDSIMF_APPUTILS_SIGNALHANDLER_HPP