#!/bin/sh

Help()
{
   # Display Help
   echo "Executes IDSimF applications inside docker container. Output files are saved to current working directory."
   echo "Docker container is called "idsimf-cli"."
   echo
   echo "Syntax: ./IDSimF.sh <directory of the executable inside build directory> 
          <directory where results should be saved to inside container>
          <executable> <config file for application> <output file name> <number of threads>"
   echo "options:"
   echo "-h     Print syntax information."
  
   echo
}

while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done

exec sudo docker run \
  -t \
  -v "$PWD":/IDSimF/build/"$1"/"$2" \
  -w /IDSimF/build/"$1" \
  "idsimf-cli" ./"$3" "$2"/"$4" "$2"/"$5" --n_threads "$6"


