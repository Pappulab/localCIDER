#!/bin/bash





update_version() {
    
    if [ ! "$#" == "2" ]
    then
	echo "FAIL - PASS PREVIOUS AND NEW VERSION"
	return
    fi
    
    PREV=$1
    NEW=$2
    # first print everything which matches
    echo " ALL MATCHING ---> "
    grep "!    Version ${PREV}" *.py
    echo "This will update the version of all python files in this directory to" $NEW

    echo "Enter Y to continue"
    read ok

    if [ "$ok" == "Y" ]
    then

	for i in *.py
	do
	    echo $i
	    echo $1
	    sed -i.bak s/"!    Version ${PREV}"/"!    Version ${NEW}"/g $i
	done
    else
	echo "ABORT"
    fi

}
