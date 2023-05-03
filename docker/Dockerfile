FROM gcc:latest

RUN apt-get update

RUN apt-get -y install git && apt-get -y install cmake
RUN apt-get -y install hdf5-tools && apt-get -y install h5utils && apt-get -y install libhdf5-dev 
RUN git clone https://github.com/IPAMS/IDSimF.git
WORKDIR /IDSimF
RUN mkdir build 
WORKDIR /IDSimF/build 
RUN cmake .. && make -j 4 


