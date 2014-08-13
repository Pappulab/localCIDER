#!/bin/bash





update_version() {
    
    if [ ! "$#" == "1" ]
    then
	echo "FAIL - PASS A VERSION"
	return
    fi
    
    # first print everything which matches
    echo " ALL MATCHING ---> "
    grep "!    Version 0.1.0" *.py
    echo "This will update the version of all python files in this directory to" $2

    echo "Enter Y to continue"
    read ok

    if [ "$ok" == "Y" ]
    then

	for i in *.py
	do
	    echo $i
	    echo $1
	    sed -i.bak s/"!    Version 0.1.0"/"!    Version $1"/g $i
	done
    else
	echo "ABORT"
    fi

}
