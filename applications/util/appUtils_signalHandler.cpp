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
 appUtils_signalHandler.cpp

 Description

 ****************************/
#include "appUtils_signalHandler.hpp"
#include <csignal>

/**
 * Sigterm handler: Calls the static terminate receiver method
 * @param _ignored
 */
void exitSignalHandler_(int _ignored) {
    AppUtils::SignalHandler::sendTerminateToReceiver();
}

/**
 * Sets a new receiver time integrator which should be terminated by SIGTERM
 */
void AppUtils::SignalHandler::setReceiver(ParticleSimulation::AbstractTimeIntegrator& receiver) {
    receiver_ = &receiver;

    if (signal((int) SIGTERM, exitSignalHandler_) == SIG_ERR)
    {
        throw SignalException("Error setting up signal handlers");
    }
}

/**
 * Send a termination signal to the receiving time integrator
 */
void AppUtils::SignalHandler::sendTerminateToReceiver() {
    receiver_->setTerminationState();
}

