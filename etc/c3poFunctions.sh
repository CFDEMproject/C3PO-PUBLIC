#!/bin/sh
###############################################
# Set of shell functions for CPPPO
###############################################

printLogo()
{
echo '                  ___   _____   _____   _____     ___                            '
echo '                 / ___\/\  __`\/\  __`\/\  __`\  / __`\                          '
echo '                /\ \__/\ \ \_\ \ \ \_\ \ \ \_\ \/\ \_\ \                         '
echo '                \ \____\\ \  __/\ \  __/\ \  __/\ \____/                         '
echo '                 \/____/ \ \ \/  \ \ \/  \ \ \/  \/___/                          '
echo '                          \ \_\   \ \_\   \ \_\                                  '
echo '                          \/_/    \/_/    \/_/                                   '
echo '                                                                                 '
echo '         A Compilation for Fluid-Particle Data Post PrOcessing                   '
echo '                                                                                 '
echo ' Copyright (C): 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria   '
echo "     2014 Graz University of Technology (ippt.tugraz.at), Graz, Austria          "
}

#========================================#
#- function to check if a directory exits
checkDir()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    filePath="$1"
    #--------------------------------------------------------------------------------#
    if [ -d "$filePath" ]; then
         echo "true" 
    else
        echo "false" 
    fi
}

#- function to check if a directory exits
checkDirComment()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    filePath="$1"
    varName="$2"
    critical="$3"
    #--------------------------------------------------------------------------------#
    if [ $(checkDir $filePath) == "true" ]; then
         echo "valid:yes critical:$critical - $varName = $filePath" 
    else
        echo "valid:NO  critical:$critical - $varName = $filePath does not exist" 
    fi
}

#========================================#
#- function to check if a variable exits
checkEnv()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    var="$1"
    #--------------------------------------------------------------------------------#
    if [[ $var == "" ]]; then
        echo "false"
    else
        echo "true"
    fi
}

#========================================#
#- function to check if a variable exits
checkEnvComment()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    var="$1"
    varName="$2"
    critical="$3"
    #--------------------------------------------------------------------------------#
    if [ $(checkEnv $var) == "true" ]; then
         echo "valid:yes critical:$critical - $varName = $var" 
    else
        echo "valid:NO  critical:$critical - $varName = $var variable not set!" 
    fi
}


#- function to check if a file exits
checkFile()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    filePath="$1"
    #--------------------------------------------------------------------------------#
    if [ -f "$filePath" ]; then
         echo "true" 
    else
        echo "false" 
    fi
}


