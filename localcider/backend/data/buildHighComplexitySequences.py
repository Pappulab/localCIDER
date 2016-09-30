""" 
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Version 0.1.9                                                         !
   !                                                                          !
   !    Copyright (C) 2014 - 2016                                             !
   !    The localCIDER development team (current and former contributors)     !
   !    Alex Holehouse, James Ahad, Rahul K. Das.                             !
   !                                                                          !
   !    localCIDER was developed in the lab of Rohit Pappu at Washington      !
   !    University in St. Louis. Please see the website for citation          !
   !    information:                                                          !
   !                                                                          !
   !    http://pappulab.github.io/localCIDER/                                 !
   !                                                                          !
   !    For more information please see the Pappu lab website:                !
   !                                                                          !
   !    http://pappulab.wustl.edu/                                            !
   !                                                                          !
   !    localCIDER is free software: you can redistribute it and/or modify    !
   !    it under the terms of the GNU General Public License as published by  !
   !    the Free Software Foundation, either version 3 of the License, or     !
   !    (at your option) any later version.                                   !
   !                                                                          !
   !    localCIDER is distributed in the hope that it will be useful,         !
   !    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
   !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
   !    GNU General Public License for more details.                          !
   !                                                                          !
   !    You should have received a copy of the GNU General Public License     !
   !    along with localCIDER.  If not, see <http://www.gnu.org/licenses/>.   !
   !--------------------------------------------------------------------------!
   ! AUTHORSHIP INFO:                                                         !
   !--------------------------------------------------------------------------!
   !                                                                          !
   ! MAIN AUTHOR:    Alex Holehouse                                           !
   !                                                                          !
   !--------------------------------------------------------------------------!

   Stand alone program to compute the maximally complex sequences associated with
   a string of length *n* made of a 20 letter alphabet. Note that while proteins
   have some bias towards particular amino acids (i.e. few Trp, lots of Leu) the
   maximally complexity sequence is a theoretical max which selects uniformly from
   the 20 amino acids. This is a distinction related to how we define complexity -
   i.e. here complexity is defined in an absolute entropic sense and not in a
   way biased towards the maximally complex sequence seen in nature. In practice
   this probably means that you will never see a natural sequence with a complexity  
   of 1, so one pratical post-processing approach might be to normalize sequences
   by the most complex natural sequence you observe to get a sense of relative 
   complexity.

   We felt that by providing a raw, entropic/information theoretic measure of 
   complexity and then suggesting that users use that in whatever means they 
   want is a much safer way, rather than assuming some biasing will apply to
   all users' needs.

   This isn't really part of localCIDER, but may be useful for theoretical calculatoins
   associated with complexity. Be wary of depleting the entropy pool, though!
   
"""


import random
import zlib

NUM_ITERATIONS = 200000
MAX_LENGTH     = 1000+1
LETTERS        = "qwertyipasdfghklcvnm"


# build a random flat sequence
maxComp=-1
maxCompDict={}
CONVERGENCE_TARGET=600
# convergence test
for i in xrange(0,0):
    maxComp=-1
    maxCompDict={}

    # seed the random number generator!
    random.seed()

    for convergenceIteration in xrange(1,200000):
        length=MAX_LENGTH

        # build random sequence
        seq=""
        for pos in xrange(0,length):
            seq=seq+random.choice(LETTERS)
        
        compressed = len(zlib.compress(seq))
        
        # update
        if maxComp < compressed:
            #print "Updatding at iter %i with %i" % (convergenceIteration,compressed)
            maxCompDict[compressed] = convergenceIteration
            maxComp = compressed

    #print maxCompDict.keys()

    if CONVERGENCE_TARGET in maxCompDict.keys():
        print "Found max after %i" % maxCompDict[CONVERGENCE_TARGET]

        
        with open("convergenceDist.dat","a") as fh:
            fh.write("%i \n" % maxCompDict[CONVERGENCE_TARGET])
    else:
        print "Did not find max"
        with open("convergenceDist.dat","a") as fh:
            fh.write("-1\n")


    
maxComplexity = {}

# for each sequence of some length between
# 5 and MAX_LENGTH
for length in xrange(5,MAX_LENGTH):

    print "On sequences of length %i " % length

    # intialize 
    maxComplexity[length] = -1

    for iteration in xrange(0, NUM_ITERATIONS):
        
        # bulild flat, random string
        seq=""
        for pos in xrange(0,length):
            seq=seq+random.choice(LETTERS)
        
        # determine sequence complexity as the length of the
        # compressed string
        compressed = len(zlib.compress(seq))
        
        # update
        if maxComplexity[length] < compressed:
            maxComplexity[length] = compressed


# once complete write to a csv

with open("sequence_complexity.dat", "w") as fh:
    # write each max complexity to file!
    for length in xrange(5, MAX_LENGTH):
        fh.write("%i \t %i\n" % (length, maxComplexity[length]))
        
print "Complexity done for sequences between 5 and %i (using %i iterations per sequence)" % (MAX_LENGTH, NUM_ITERATIONS)


    
    

