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
   ! MAIN AUTHOR:   James Ahad and Alex Holehouse                             !
   !                                                                          !
   !--------------------------------------------------------------------------!



   File Description:
   ================

   Class which carrys out all permutation and density of states calculations.

   This class also implements all the Wang-Landau algorithms, both normal
   WL, and HistogramZoom WL.

   This functionality is NOT YET OPERATIONAL and will be introduced in
   version 0.2.0

"""


import sys
import numpy as np
import random as rng
import time as t
import time
import os

from backendtools import running_dotdotdot
from sequence import Sequence

from localciderExceptions import WLException


class WangLandauMachine:
    """
    The WangLandauMachine is the main object which carries out Wang-Landau Monte Carlo to generate random permutants and estimates the density of states.

    It should be interfaced through the sequencePermutants API in the main localCIDER directory, not here

    """
    #...................................................................................#

    def __init__(
            self,
            seq,
            writedir,
            frozenResidues=set(
                []),
            nbins=10,
            binmin=0,
            binmax=1,
            flatchk=10000,
            flatcrit=.7,
            convergence=np.exp(.000001),
            WL_type="NORMAL"):

        # the self.seq must be a sequence object
        if seq.__class__.__name__ == "Sequence":
            self.seq = seq                     # Input was already a sequence object
        else:
            # Convert a string into a sequence object
            self.seq = Sequence(seq)

        # set some general variables
        self.writeDir = writedir               # Writable output dir
        # Frequency at which a flatcheck is performed
        self.nflatchk = int(flatchk)
        self.flatcrit = float(flatcrit)        # Criterion for flatness
        self.convergence = float(convergence)  # Convergence threshold
        self.WL_type = WL_type

        # Set the bin resolution for the flatness criterion region required
        # In both HistogramZoom and NormalWL
        # --> This will remain constant through the search for both normal
        #     and HistogramZoom modes
        self.nbins_target = int(nbins)

        # binmin and binmax:
        # Set the bin min and bin max values (kappa space) in the full kappa histogram
        # In HistogramZoom:
        # --> Start at 0 and 1, and we slowly converge on the region of kappa
        #     density
        # In NormalWL:
        # --> Start at whatever values defined in the input
        #
        # nbins_actual: actual number of bins the 0-1 histogram is divided into
        # In HistogramZoom
        # --> This will change through the search to accomodate a constant
        #     $target_bins in a relevant region of kappa space
        # In NormalWL
        # --> Like the actual number of bins this remains constant
        #
        # relevantmin and relevantmax: Define the index values into the G/H
        # histograms which correspond to the binmin and binmax values
        # In HistogramZoom:
        # --> These always start as the first and last index in the histogram
        # In NormalWL:
        # --> Find the bin centers which contain the binmin and binmax centers
        #

        if WL_type == "ZOOM":
            # Set the initial bin mim
            self.binmin = 0.0  # Max bin value
            self.binmax = 1.0  # Min bin value

            self.nbins_actual = self.nbins_target

            self.relevant_min = 0
            self.relevant_max = int(nbins) - 1

        # do the same for the non-adaptive WL histogram
        else:
            self.binmin = float(binmin)
            self.binmax = float(binmax)

            diff = self.binmax - self.binmin
            binWidth = diff / float(self.nbins_target)

            self.nbins_actual = int(round(1 / binWidth))

            bincts = self.getBinCenters()

            self.relevant_min = np.argmin(
                abs(bincts - (self.binmin + (binWidth / 2))))
            self.relevant_max = (self.relevant_min + self.nbins_target) - 1

        # set a freezefile if we have one
        self.frozen = set(frozenResidues)
        """
        Freezefile will be introduced in version 0.2.0 to allow for executable reading
        of a freeze file for the localCIDER commandline interface
        if(not frzfile == ''):
            self.frozen = self.parseFreezeFile(frzfile)
        else:
            self.frozen = set([])
        """

        # for ... during output
        self.setDotFreq()

        # print the sequence in a sensible way...
        if len(self.seq) > 30:
            print 'Sequence                  : ' + str(self.seq[0:10]) + "..." + str(self.seq[-10:len(self.seq)])
        else:
            print 'Sequence                  : ' + str(self.seq)

        if len(self.frozen) == 0:
            print 'Frozen Residues           : [NONE]'
        else:
            print 'Frozen Residues           : ' + str(self.frozen)
        print 'initial nbins             : ' + str(self.nbins_actual)
        print 'bins range in kappa space : ' + str(self.binmin) + " to " + str(self.binmax)
        print 'Writing output to         : ' + str(os.path.abspath(self.writeDir))
        print "+-------------------------------------------------+"
        print '|   WANG-LANDAU MACHINE INITIALIZED SUCCESFULLY   |'
        print "+-------------------------------------------------+"

    #...................................................................................#
    def parseFreezeFile(self, frzfile):
        """
        Function which parses freezefiles.

        Freeze files define a region of sequence which is *not* mutated
        during the WL process. The freezefile defines a region to ignore
        using two sequence positions, which are inclusive of the first
        e.g.
        10 15 would ignore positions 10,11,12,13,14

        =========================================== START OF FILE
        # Comments can be added following a hash symbol
        1 1     # comments can also be inline

        # spaces are allowed

        10 15


        """
        frozen = set()
        with open(frzfile) as f:

            # for each line in the fole
            for line in f:

                line = line.strip()
                if line[0] == "#":
                    # if we're on a comment line
                    continue

                # lets you have inline comments too..
                line = line.split("#")[0]

                # clean up a bit more and split by space
                line = line.strip()
                columns = line.split(' ')

                # check we can make a sensible range
                if columns[0] > columns[1]:
                    raise SequenceError(
                        "FREEZEFILE ERROR: First index position in freeze file is greater than second (" +
                        str(line) +
                        ")")

                # sequence indices start at 0 in python, but 1 in real world
                # so we offset the first position by -1
                frozen |= set(np.arange(int(columns[0]) - 1, int(columns[1])))

        return frozen

    #...................................................................................#
    def run(self):
        """ Main function which launches the WL job as initialized """
        if self.WL_type == "NORMAL":
            return self.run_normal_WL()
        else:
            return self.run_histogramZoom_WL()

    #...................................................................................#
    def sanity_check(self):
        """
        These algorithms are really complicated so we have a set of sanity checks here...
        """
        pass  # have to write specific types of checks for each type of WL
        # if not self.nbins_target =

    #...................................................................................#
    def setDotFreq(self):
        """
        Function which determines the dotdotdot frequency (makes ... as the program progresses!)
        """

        self.dotdotfreq = int(self.nflatchk / 20)
        if self.dotdotfreq < 1:
            self.dotdotfreq = 1

    #...................................................................................#
    def getBinSize(self):
        """
        Function which returns the size of each bin in the system as defined by the
        actual number of bins.
        """

        # this works because kappa space is ALWAYS divided into some number of bins between
        # 0 and 1, but we often only care about a subset of those bins. However, this does
        # mean that the bin width is always 1/the actual number of bins
        return (1.0 / float(self.nbins_actual))

    #...................................................................................#
    def indexInsideRelevantRegion(self, idx):
        """
        Function which evaluates if some integer index (idx) is inside the relevant
        region on the full kappa histogram (i.e. inside Hlocal), where Hlocal is
        defined by self.relevant_min and self.relevant_max

        """
        if (idx <= self.relevant_max and idx >= self.relevant_min):
            return True
        else:
            return False

    #...................................................................................#
    def getBinCenters(self):
        """
        Function which returns the center of the actual bin vector
        """
        binsz = self.getBinSize()

        return binsz / 2 + binsz * np.arange(0, self.nbins_actual)

    #...................................................................................#
    def run_histogramZoomWL(self, zoomcount=1):
        """
        This is the main WL script which is used to estimate the density of states of kappa using
        the histogramZoom algorithm

        It is *heavily* annotated, to try and be as transparent and easy to understand as possible.

        """

        restart = True
        globalStartTime = t.time()

        # This main while loop runs until restart is set to false
        # which only happens when we reach convergence
        while(restart):

            # Initialize WL histogram information
            bincts = self.getBinCenters()    # bincts = histogram bin centers
            # g = density of state (DOS) histogram
            g = [0] * self.nbins_actual
            H = [0] * self.nbins_actual        # H = kappa histogram
            f = np.exp(1)                    # f = convergence test
            seqcount = 0                     # freq at which sequences are written out
            nstep = 0                        # initialize number of steps
            # initialize number of iterations with succesfull flatchecks
            niter = 0
            flatcount = 0                    # initialize number of cycles past the  value which have
            # failed to reach flatness

            # (re)seed the RNG
            rand = rng.Random()
            rand.seed(t.time())

            # Create output files for recording progess (overwrites if
            # the file already exists - this happens on a restart)
            hlog = self.mklog(os.path.join(self.writeDir, "hlog.txt"))
            glog = self.mklog(os.path.join(self.writeDir, "glog.txt"))
            seqlog = self.mklog(os.path.join(self.writeDir, "seqlog.txt"))
            hblog = self.mklog(
                os.path.join(
                    self.writeDir,
                    "histogram_bins.txt"))
            hlocalblog = self.mklog(
                os.path.join(
                    self.writeDir,
                    "local_histogram_bins.txt"))

            # Write initial header information to those files
            self.writeLog(hlog, "flchk #\tHistograms\n")
            self.writeLog(hlog, "iter 1:\n")
            self.writeLog(glog, "iter\tdensity of states\n")
            self.writeLog(seqlog, "Kappa\tSequence\n")

            # Write histogram centers
            self.writeLog(hblog, self.fprintVertVector(bincts))
            self.writeLog(
                hblog,
                self.fprintVertVector(
                    bincts[
                        self.relevant_min:self.relevant_max +
                        1]))

            # initial starting conditions
            oseq = self.seq
            kold = oseq.kappa()

            # get histogram bin index of this sequence's kappa value
            idx_old = np.argmin(abs(bincts - kold))

            startTime = t.time()

            while(f > self.convergence):

                if nstep % self.dotdotfreq == 0:
                    running_dotdotdot()

                # There are two possible Monte Carlo moves, which we alternate
                # between at a 80:20 split
                # - Random swapping of two residues which changes the kappa (80%)
                # - Complete unbiased shuffle of sequence (20%)

                if rand.random() > 0.8:
                    nseq = oseq.swapRandChargeRes(self.frozen)
                else:
                    nseq = oseq.full_shuffle(self.frozen)

                # calculate the kappa of that new sequence
                knew = nseq.kappa()

                # get the histogram index of that new kappa
                idx_new = np.argmin(abs(bincts - knew))

                ###############################################################
                #                WL Acceptance Criterion                      #
                ###############################################################

                # From http://www.pages.drexel.edu/~cfa22/msim/node53.html

                # we accept IF the
                # A) there are fewer sequences in the new bin than the old bin (DOS bins)
                # B) with an exponential decay factor relative the difference if there are more in the
                #    the old bind than the new bin
                #
                # Note if we're outside the region of interest we evaluate based on the Histogram densities
                # instead of the densities of states

                # print self.relevant_min, self.relevant_max, idx_new
                if self.indexInsideRelevantRegion(idx_new):
                    acceptProb = min([1, np.exp(g[idx_old] - g[idx_new])])
                else:
                    acceptProb = min([1, np.exp(H[idx_old] - H[idx_new])])

                # print(acceptProb)
                # if new sequence kappa is less visited than old sequence
                # kappa, visit it

                # =============================================================================
                # ACCEPTANCE REGION
                # if we accept the move
                if(rand.random() < acceptProb):

                    # We write every seqlen^2  sequence
                    # to the sequence log
                    if(seqcount == 0):
                        self.writeLog(
                            seqlog, '%0.3f\t%s\n' %
                            (oseq.kappa(), oseq))
                        seqcount = int(oseq.len**2)
                    else:
                        seqcount = seqcount - 1

                    # set the oldsequence to the sequence we just accepted and
                    # simmilarly update the kappa old and bin index old
                    oseq = Sequence(nseq.seq, nseq.dmax, nseq.chargePattern)

                    # TO DO test if oseq.kappa() is ever not equal to knew??
                    kold = oseq.kappa()
                    idx_old = np.argmin(abs(bincts - kold))

                    # reset the new sequence and new sequence histogram index
                    nseq = None
                    idx_new = 0
                else:
                    nseq = None
                    idx_new = 0

                # END OF ACCEPTANCE REGION
                # =============================================================================

                g[idx_old] = g[idx_old] + np.log(f)
                H[idx_old] = H[idx_old] + 1

                # increment the number of steps taken
                nstep = nstep + 1

                # if we're at a flatcheck
                if(nstep == self.nflatchk):

                    print ""
                    print "Interest indicies = " + str(self.relevant_min) + " and " + str(self.relevant_max)

                    Hlocal = H[self.relevant_min:self.relevant_max + 1]

                    print "On flatcheck number:  " + str(flatcount)
                    flatcount += 1
                    self.writeLog(hlog, str(flatcount) + "\t" +
                                  self.fprintHVector(Hlocal) + "\n")

                    (H, f, niter, nstep) = self.__run_flatcheck(
                        H, Hlocal, niter, f, hlog, glog, g)

                    # to timing stats for this set of cycles
                    endTime = t.time()
                    print("Time for this flat-check cycle: " +
                          str(endTime - startTime))
                    startTime = t.time()

                # if this is the 20th flatcheck and we've not had a single succesfull iteration
                # then time to ZOOM
                if (flatcount >= 20 and niter == 0):
                    break

            # This is OUTSIDE the convergence while loop
            if (flatcount >= 20 and niter == 0):

                # create the local kappa histogram
                Hlocal = H[self.relevant_min:self.relevant_max + 1]

                # run the histogram zoom algorithm and restart
                zoomcount = self.histogram_zoom(Hlocal, H, zoomcount)
                flatcount = 0

            else:
                # if we're here we're done, so turn off the
                # restart
                restart = False

        # finalize and write up
        dos = open(os.path.join(self.writeDir, "DOS.txt"), 'w')
        dos.write('bincts\tlog(omega)\n')
        for i in range(len(bincts)):
            dos.write("%0.3f\t%5.6f\n" % (bincts[i], g[i]))
        dos.close()
        dos = open(os.path.join(self.writeDir, "DOS_local.txt"), 'w')
        dos.write('bincts\tlog(omega)\n')
        for i in range(self.relevant_min, self.relevant_max + 1):
            dos.write("%0.3f\t%5.6f\n" % (bincts[i], g[i]))
        dos.close()

        globalEndTime = t.time()
        print("Total Run Time in Seconds")
        print(globalEndTime - globalStartTime)
        return np.vstack((bincts, g))

    #...................................................................................#
    def run_normal_WL(self):
        """
        This is the main WL script which is used to estimate the density of states of kappa using the normal
        Wang Landau algorithm

        It is *heavily* annotated, to try and be as transparent and easy to understand as possible.
        """

        globalStartTime = t.time()

        # WL histogram information
        bincts = self.getBinCenters()    # bincts = histogram bin centers
        # g = density of state (DOS) histogram
        g = [0] * self.nbins_actual
        H = [0] * self.nbins_actual        # H = kappa histogram
        f = np.exp(1)                    # f = convergence test
        # initialize seqcount (determines sequence writing)
        seqcount = 0
        nstep = 0                        # initialize number of steps
        # initialize number of iterations (with succesfull flatchecks)
        niter = 0
        flatcount = 0                    # initialize number of cycles past the
        # flatchk value which have failed to reach flatness
        self.sanity_check()

        # RNG
        rand = rng.Random()
        rand.seed(t.time())

        # Create output files for recording progess (this overwrites if
        # the file already exists, which happens on a restart)
        hlog = self.mklog(os.path.join(self.writeDir, "hlog.txt"))
        glog = self.mklog(os.path.join(self.writeDir, "glog.txt"))
        seqlog = self.mklog(os.path.join(self.writeDir, "seqlog.txt"))
        hblog = self.mklog(os.path.join(self.writeDir, "histogram_bins.txt"))

        # Write initial header information to those files
        self.writeLog(hlog, "flchk #\tHistograms\n")
        self.writeLog(hlog, "iter 1:\n")
        self.writeLog(glog, "iter\tdensity of states\n")
        self.writeLog(seqlog, "Kappa\tSequence\n")

        # Write histogram centers
        self.writeLog(hblog, self.fprintVertVector(bincts))
        self.writeLog(
            hblog,
            self.fprintVertVector(
                bincts[
                    self.relevant_min:self.relevant_max +
                    1]))

        # initial starting conditions
        oseq = self.seq
        kold = oseq.kappa()

        # get histogram bin index of original (old) sequence kappa
        idx_old = np.argmin(abs(bincts - kold))

        startTime = t.time()
        reject = 0

        # This main while loop runs until we reach convergence. Note that depending on how long your sequence is
        # this might run for a while...
        while(f > self.convergence):

            if nstep % self.dotdotfreq == 0:
                running_dotdotdot()

            # There are two possible Monte Carlo moves, which we alternate
            # between at a 70:30 split
            # - Random swapping of two residues which changes the kappa (80%)
            # - Complete unbiased shuffle of sequence (20%)

            if rand.random() > 0.8:
                nseq = oseq.swapRandChargeRes(self.frozen)
            else:
                nseq = oseq.full_shuffle(self.frozen)

            # calculate the kappa of that new sequence
            knew = nseq.kappa()

            # get the histogram index of that new kapp
            idx_new = np.argmin(abs(bincts - knew))

            ###############################################################
            #                WL Acceptance Criterion                      #
            ###############################################################

            # From http://www.pages.drexel.edu/~cfa22/msim/node53.html

            # we accept IF the
            # A) there are fewer sequences in the new bin than the old bin (DOS bins)
            # B) with an exponential decay factor relative the difference if there are more in the
            #    the old bind than the new bin
            #
            # Note if we're outside the region of interest we automatically reject, so we're confining
            # states of interest to be inside the kappa region of interest and will NEVER accept a move
            # outside that region. This has a relfective boundary property.

            skip = False

            # print self.relevant_min, self.relevant_max, idx_new
            if self.indexInsideRelevantRegion(idx_new):
                acceptProb = min([1, np.exp(g[idx_old] - g[idx_new])])
            else:
                reject = reject + 1
                acceptProb = 0
                skip = True

            # print(acceptProb)
            # if new sequence kappa is less visited than old sequence kappa,
            # visit it

            # =============================================================================
            # ACCEPTANCE REGION
            # if we accept the move
            if(rand.random() < acceptProb):

                #
                # I think the following if/else means we write every seqlen^2
                # sequence to the sequence log, though its honestly not clear to me
                # *why* we should do this...
                #

                if(seqcount == 0):
                    self.writeLog(seqlog, '%0.3f\t%s\n' % (oseq.kappa(), oseq))
                    seqcount = int(oseq.len**2)
                else:
                    seqcount = seqcount - 1

                    # set the oldsequence to the sequence we just accepted and
                    # simmilarly update the kappa old and bin index old

                # update the old sequence with the info from the new sequence
                oseq = Sequence(nseq.seq, nseq.dmax, nseq.chargePattern)
                kold = oseq.kappa()
                idx_old = np.argmin(abs(bincts - kold))

                # reset the new sequence and new sequence histogram index
                nseq = None
                idx_new = 0

            # if we do not accept the move
            else:
                nseq = None
                idx_new = 0
            # END OF ACCEPTANCE REGION
            # =============================================================================

            # increment the DOS histogram and the energy histogram (note if we rejected the move
            # we're staying where we are in kappa space)
            if not skip:
                g[idx_old] = g[idx_old] + np.log(f)
                H[idx_old] = H[idx_old] + 1

            # increment the number of steps taken
            nstep = nstep + 1

            # if we're at a flatcheck
            # TODO add some kind of 'this is taking WAY TOO LONG' check
            if(nstep % self.nflatchk == 0):

                print ""
                print "Interest indicies = " + str(self.relevant_min) + " and " + str(self.relevant_max)

                print "Rejected [outside histogram] = " + str(reject / float(self.nflatchk) * 100) + "%"
                reject = 0

                Hlocal = H[self.relevant_min:self.relevant_max + 1]

                print "On flatcheck number:  " + str(flatcount)
                flatcount += 1
                self.writeLog(
                    hlog,
                    str(flatcount) +
                    "\t" +
                    self.fprintHVector(Hlocal) +
                    "\n")

                (H, f, niter, nstep) = self.__run_flatcheck(
                    H, Hlocal, niter, f, hlog, glog, g)
                # to timing stats for this set of cycles
                endTime = t.time()
                print("Time for this flat-check cycle: " +
                      str(endTime - startTime))
                startTime = t.time()

        # Finalize and clean up
        dos = open(os.path.join(self.writeDir, "DOS.txt"), 'w')
        dos.write('bincts\tlog(omega)\n')
        for i in range(len(bincts)):
            dos.write("%0.3f\t%5.6f\n" % (bincts[i], g[i]))
        dos.close()
        dos = open(os.path.join(self.writeDir, "DOS_local.txt"), 'w')
        dos.write('bincts\tlog(omega)\n')
        for i in range(self.relevant_min, self.relevant_max + 1):
            dos.write("%0.3f\t%5.6f\n" % (bincts[i], g[i]))
        dos.close()

        globalEndTime = t.time()
        print("Total Run Time in Seconds")
        print(globalEndTime - globalStartTime)
        return np.vstack((bincts, g))
        # close while loop

    #...................................................................................#
    def run_permutant_generator(self, numberOfSequences):
        """
        Modified normal WL which instead of writing density of state information and logfiles will simply generate a number of kappa permutants for some bin values

        """

        # This main while loop runs until restart is set to false
        # which only happens when we reach convergence
        globalStartTime = t.time()

        # WL histogram information
        bincts = self.getBinCenters()    # bincts = histogram bin centers
        # g = density of state (DOS) histogram
        g = [0] * self.nbins_actual
        H = [0] * self.nbins_actual        # H = kappa histogram
        f = np.exp(1)                    # f = convergence test
        nstep = 0                        # initialize number of steps

        # construct an empty list of lists, which we'll populate with $numberOfSequences sequences at
        # each bin position
        sequence_vector = []
        for i in xrange(0, self.nbins_actual):
            sequence_vector.append([])

        # sequence vector is now a list of empty lists equal to the full kappa
        # histogram

        # RNG
        rand = rng.Random()
        rand.seed(t.time())

        # initial starting conditions
        oseq = self.seq
        kold = oseq.kappa()

        # get histogram bin index of original (old) sequence kappa
        idx_old = np.argmin(abs(bincts - kold))

        startTime = t.time()
        reject = 0
        binsfull = False

        while(binsfull == False):

            if nstep % self.dotdotfreq == 0:
                running_dotdotdot()

            # There are two possible Monte Carlo moves, which we alternate
            # between at a 70:30 split
            # - Random swapping of two residues which changes the kappa (80%)
            # - Complete unbiased shuffle of sequence (20%)

            if rand.random() > 0.8:
                nseq = oseq.swapRandChargeRes(self.frozen)
            else:
                nseq = oseq.full_shuffle(self.frozen)

            # calculate the kappa of that new sequence
            knew = nseq.kappa()

            # get the histogram index of that new kapp
            idx_new = np.argmin(abs(bincts - knew))

            ###############################################################
            #                WL Acceptance Criterion                      #
            ###############################################################

            # From http://www.pages.drexel.edu/~cfa22/msim/node53.html

            # we accept IF the
            # A) there are fewer sequences in the new bin than the old bin (DOS bins)
            # B) with an exponential decay factor relative the difference if there are more in the
            #    the old bind than the new bin
            #
            # Note if we're outside the region of interest we automatically reject, so we're confinig
            # states of interest to be inside the kappa region of interest and will NEVER accept a move
            # outside that region

            skip = False

            # print self.relevant_min, self.relevant_max, idx_new
            if self.indexInsideRelevantRegion(idx_new):
                acceptProb = min([1, np.exp(g[idx_old] - g[idx_new])])

                if len(sequence_vector[idx_new]) < numberOfSequences:
                    sequence_vector[idx_new].append(nseq.seq)
                else:
                    # randomly replace one of the sequences in the bins
                    sequence_vector[idx_new][rand.randint(
                        0, numberOfSequences)] = nseq.seq
            else:
                reject = reject + 1
                acceptProb = 0
                skip = True

            # print(acceptProb)
            # if new sequence kappa is less visited than old sequence kappa,
            # visit it

            # =============================================================================
            # ACCEPTANCE REGION
            # if we accept the move
            if(rand.random() < acceptProb):

                # update the old sequence with the info from the new sequence
                oseq = Sequence(nseq.seq, nseq.dmax, nseq.chargePattern)
                kold = oseq.kappa()
                idx_old = np.argmin(abs(bincts - kold))

                # reset the new sequence and new sequence histogram index
                nseq = None
                idx_new = 0

            # if we do not accept the move
            else:
                nseq = None
                idx_new = 0
            # END OF ACCEPTANCE REGION
            # =============================================================================

            # increment the DOS histogram and the energy histogram (note if we rejected the move
            # we're staying where we are in kappa space)
            if not skip:
                g[idx_old] = g[idx_old] + np.log(f)
                H[idx_old] = H[idx_old] + 1

            # increment the number of steps taken
            nstep = nstep + 1

            # check bin full-ness every 1000 steps
            if nstep % 1000 == 0:

                # print number rejected outside region of interest
                print "Rejected [outside histogram] = " + str(float(reject) / float(1000)) + "%"
                reject = 0

                # calculate the total number of sequences we have
                totalNumSeq = 0
                for i in sequence_vector:
                    totalNumSeq = totalNumSeq + len(i)

                # if the total number of sequences is equal to
                # the target number then we're done
                if ans == numberOfSequences * self.nbins_target:
                    binsfull = True

    #...................................................................................#
    def __run_flatcheck(self, H, Hlocal, niter, f, hlog, glog, g):
        """
        Function which analyzes the local histogram region of interest to evaluate flatness.
        """

        print "H local     = " + str(Hlocal)
        print "bin regions = [" + str(self.getBinCenters()[self.relevant_min]) + "] to [" + str(self.getBinCenters()[self.relevant_max]) + "]"

        print "Full Kappa histo: " + str(H)
        print "Full DOS        : " + str(g)

        # If everything is flat (i.e. all bins are equally full then
        # each bin in H (the vector of the number of sequences in each bin)
        # divided by the mean will be 1
        #
        # As you have increasing lack of flatness you have more and more which
        # are further away from the mean in a symetric fashion
        #

        flatness_number = len(
            np.where(
                Hlocal /
                np.mean(Hlocal) >= self.flatcrit)[0])

        print(
            "Number of histograms which meet flat check criterion: " +
            str(flatness_number))

        if flatness_number == self.nbins_target:
            # If bins are all flat!
            print ""
            print "+========================================================+"
            print "| Flatcheck criterion was met - resetting the histograms |"
            print "+========================================================+"
            print ""

            # update the convergence criterion
            f = f**0.5

            # reset ALL the bins
            H = [0] * self.nbins_actual

            # incremenet the iteration counter and print
            niter = niter + 1
            print("Now on iteration " + str(niter))
            print("f = " +
                  str(f) +
                  " convergence set at " +
                  str(self.convergence))
            self.writeLog(
                glog,
                str(niter) +
                "\t" +
                self.fprintGVector(g) +
                "\n")
            if(f > self.convergence):
                self.writeLog(hlog, "\niter %d:\n" % (niter + 1))

        return(H, f, niter, 0)

    #...................................................................................#
    def histogram_zoom(self, Hlocal, H, zoomcount):
        """

           Histogram Zoom: Introduction
           ===========================

           Histogram zoom is an algorithm I developed to try and ensure arbitrarily good resolution
           for the kappa density of states, while avoiding artefacts introduced by poor binsize
           selection.

           It centers around two key ideas

           1) The distribution of kappa values accross all sequence permutants of a given sequence
              must be continous. That is to say, you cannot generate a historgram which looks like this

               |
               |             #
          DOS  | #           #
               | #           #      #
               | # # # # # # #      #
               |_#_#_#_#_#_#_#_#____#_#___
                     KAPPA

              The important point to take away is that we have a region of 0 density surrounded by
              two regions of high density. This cannot happen - i.e. you cannot have a total absence
              of DOS but instead there must be a smoothe transition.

              You can have 0 (or near 0) density on the edges, but not in the middle. Because kappa
              is a self-consistent parameter there is by definition never any density histogram (assuming
              a wide enough bin) which has a zero-density at a given value assuming full permutant space
              sampling. If you DO see zero space this is a consequence of bin widths being too narrow
              meaning numerically you can't achieve certain values of kappa. This non-continous
              behaviour is a numerical artefact and not representative of the actuall density of states.

          2) When assessing histogram flatness, the histogram bins of interest MUST be the same size


          Algorithm design
          ================
          * Note: variables discussed in prose are prefixed by a '$' symbol - so, "increment the value" would
            mean increment some sentence-context value, whereas "increment the $value" would mean increment
            the variable value
          *

          Considering this, the algorithm has the following approach

          1) After 20 iterations of an appropriate number of steps if flatness was not reached we call this
             function, including how many times without flatness success we've had ($zoomcount)

          2) Scan through each of the histogram bins and identify those which are 1% * $zoomcount of the
             most full bind. These are our designated "bad_indexes"

          3) We then shrink the local region of interest to halfway between the bad index bins and the
             previous bin max and bin min

          4) There are also a bunch of heuristics like extending the binds if there's a lot of density
             just outside the max or min bin of interest


        Essentially, this lets us slowly zoom in on the region of the full histogram where the majority
        of the density is, importantly, where the degree of zoom is inversly proportional to the number
        of steps before a flatcheck. This means you can essentially estimate the region of density to
        arbitrary precision, giving a tradeoff of computation time vs. true density.

        This is very useful to rapidly assess if your sequence's kappa and the [effective] maximum likelihood
        kappa of that sequence composition are close, without needing to estimate the FULL density of states.

        It also lets you guage where you might want to set density region boundaries for carrying out true
        Wang Landau

        """

        # IF WE'RE DOING NORMAL WL JUST RETURN. NOTHING TO SEE HERE...
        if self.WL_type == "NORMAL":
            return

        # first thing - is the local region failing to capture all the density? We use a nearest neighbour
        # check...
        newzoom = self.__check_local_skew(zoomcount, H)

        if newzoom == zoomcount:
            print "Shifting region of interest in histogram"
            return zoomcount

        high = max(Hlocal)

        # so initialize the local histogram indices to the extremes of
        # the local histogram vector. new_end and new_start refer to position in
        # the Hlocal, NOT the full H
        new_start = 0
        new_end = self.nbins_target - 1

        # search for bins that had less than 1%*zoomcount of the highest bin
        # evertime we meet the flatcheck criterion we reset the zoomcount to 1
        bad_indices = []
        index = 0
        for i in Hlocal:
            if i < (zoomcount * 0.01) * high:
                bad_indices.append(index)
            index = index + 1

        if len(bad_indices) == 0:
            "Local vector appropriatly full - trying again..."
            return zoomcount + 1

        # determine the first and last index position in Hlocal where we have density
        # if if Hlocal = [0,0,3,4,5,3,0,0]
        # then
        # ---> new_start = 1
        # ---> new_end   = 6
        # pushing_min/max are set to true if there's density up to the edge of the HLocal
        # histogram ("pushing against the maximimum or pushing against the
        # minimum)
        new_start, pushing_min = self.__get_new_Hlocal_start(bad_indices)
        new_end, pushing_max = self.__get_new_Hlocal_end(bad_indices)

        print ">> bad indicies = " + str(bad_indices)
        print ">> New start    = " + str(new_start)
        print ">> New end      = " + str(new_end)

        # new_start and new_end are relevant for the local histogram. Next
        # lets move out of the local histogram into the full histogram
        # to determine where in the full histogram those indicies lie

        full_histogram_new_start = self.relevant_min + new_start
        full_histogram_new_end = self.relevant_max - \
            (self.nbins_target - new_end)

        # we have to do some sanity checking here - specifically if we're pushing
        # up againts the edges of the histogram...
        if full_histogram_new_start < 0:
            full_histogram_new_start = 0

        if full_histogram_new_end > (self.nbins_actual - 1):
            full_histogram_new_end = self.nbins_actual - 1

        # so now we have the index value into the FULL histogram, lets translate
        # that index value into kappa-space
        new_min_kappa = full_histogram_new_start * self.getBinSize()
        new_max_kappa = full_histogram_new_end * self.getBinSize()

        old_start_kappa = self.relevant_min * self.getBinSize()
        old_end_kappa = self.relevant_max * self.getBinSize()

        # now figure out 1/2 way beween the new kappa and the old kappa min
        # and max values
        new_start_kappa = (
            (new_min_kappa - old_start_kappa) / 2) + old_start_kappa
        new_end_kappa = old_end_kappa - ((old_end_kappa - new_max_kappa) / 2)

        print ">> old start kappa  = " + str(old_start_kappa)
        print ">> old end kappa    = " + str(old_end_kappa)
        print "---------------------------------------"
        print ">> new start kappa  = " + str(new_start_kappa)
        print ">> new end kappa    = " + str(new_end_kappa)

        # how big is this new region of interest in kappa space?
        kappa_distance = new_end_kappa - new_start_kappa

        # how wide must bins be to get self.nbins_target in this region?
        new_binWidth = kappa_distance / float(self.nbins_target)

        # so now many bins do we have to reset the FULL histogram to for
        # this resolution?
        self.nbins_actual = int(1.0 / new_binWidth)
        bin_centers = self.getBinCenters()

        print "new bin centers" + str(bin_centers)

        # finally determine the index values for the relevant regions based on the
        # new_start_kappa2 and new_end_kappa2 values
        self.relevant_min = np.argmin(abs(bin_centers - new_start_kappa))
        self.relevant_max = np.argmin(abs(bin_centers - new_end_kappa))

        print "Index for min = " + str(self.relevant_min)
        print "Index for max = " + str(self.relevant_max)

        while (
                not len(
                    bin_centers[
                self.relevant_min:self.relevant_max +
                1]) == self.nbins_target):

            # this function corrects the bin size in a sensnible way
            self.__correct_for_numerical_errors(
                new_start_kappa,
                new_end_kappa,
                pushing_min,
                pushing_max,
                bin_centers)

        print "After bin index update bin indices are: " + str(self.relevant_min) + " and " + str(self.relevant_max)
        return zoomcount

    #...................................................................................#
    def __correct_for_numerical_errors(
            self,
            new_start_kappa,
            new_end_kappa,
            pushing_min,
            pushing_max,
            bin_centers):
        # numerica errors can be introduced when determining the number of bins and whatnot
        # So we use this loop to ensure the HLocal is always the correct number
        # of bins
        print "Before bin index update bin indices are: " + str(self.relevant_min) + " and " + str(self.relevant_max)

        lowerDeviation = bin_centers[self.relevant_min] - new_start_kappa
        upperDeviation = bin_centers[self.relevant_max] - new_end_kappa

        # if we have more bins than the target number of bins
        if len(
            bin_centers[
                self.relevant_min:self.relevant_max +
                1]) > self.nbins_target:

            # Whatever happens we want to either incrase the relevant_min or decrease
            # the relevant max. To choose which we find out which is bigger;
            #
            # min_kappa - bin_centers[self.relevant_min]
            # bin_centers[self.relevant_max] - max_kappa
            #
            # if the former then we increase self.relevant_min
            # if the later we decrease self.relevant_max

            if pushing_max:
                # want to skew the shift to make the max (relativly) bigger, so
                # making the min larger is the lesser of two evils
                self.relevant_min = self.relevant_min + 1
            elif pushing_min:
                # want to make the min smaller, so making the max smaller is
                # better than making the min bigger
                self.relevant_max = self.relevant_max - 1
            else:
                # if the lower kappa bin deviates from the ideal kappa min the
                # most
                if -lowerDeviation > upperDeviation:
                    self.relevant_min = self.relevant_min + 1
                else:
                    self.relevant_max = self.relevant_max - 1
        # we have too few bins
        else:

            # Whatever happens we want to either decrease the relevant_min or increase
            # the relevant max. To choose which we find out which is bigger;
            #
            # min_kappa - bin_centers[self.relevant_min]
            # bin_centers[self.relevant_max] - max_kappa
            #
            # if the former then we increase self.relevant_min
            # if the later we decrease self.relevant_max

            if pushing_min:
                self.relevant_min = self.relevant_min - 1
            elif pushing_max:
                self.relevant_max = self.relevant_max + 1
            else:
                if lowerDeviation > -upperDeviation:
                    self.relevant_min = self.relevant_min - 1
                else:
                    self.relevant_max = self.relevant_max + 1

    #...................................................................................#
    def __check_local_skew(self, zoomcount, H):
        """
            This function takes the full kappa histogram (H). It examines the status of the histogram bins either
            side of the region of interest. If there is density in those bins then rather than carrying out a zoom
            operation, it might be better to just move the region of interest, which is what the function then
            does (by changing the position of the self.relevant_min and self.relevant_max

            The function returns the zoomcount if this re-sizing operation occurs, or -1 if not

        """

        if self.relevant_min - 1 > -1:
            # ensure we don't step beyond the 0th element
            if not H[self.relevant_min - 1] == 0:
                # if the element preceding the 0th element in the HLocal is not empty then skewshift
                # Hlocal down a set of bins and retry
                self.relevant_min = self.relevant_min - 1
                self.relevant_max = self.relevant_max - 1
                return zoomcount

        if self.relevant_max + 1 < len(H):
            # ensure we don't step beyond the end of the histogram vector
            if not H[self.relevant_max + 1] == 0:
                # if the element after the final element in the HLocal is not empty then skewshift
                # Hlocal down a set of bins and retry
                self.relevant_min = self.relevant_min + 1
                self.relevant_max = self.relevant_max + 1
                return zoomcount

        return -1

    #...................................................................................#
    def __get_new_Hlocal_end(self, bad_indices):

        # if the last bad index is the end of the histogram, work out
        # how many contigous bins there are (going backwards) which are
        # also bad
        #
        # The new_end index is the INDEX into the first empty position
        # in Hlocal from the RHS (or is one outside to the right if that
        # final position isn't empty

        if bad_indices[len(bad_indices) - 1] == (self.nbins_target - 1):
            index = self.nbins_target - 1
            for i in reversed(bad_indices):
                if i != index:
                    break
                else:
                    index = index - 1
            new_end = index + 1
            pushing_max = False
        else:
            new_end = self.nbins_target
            pushing_max = True

        return (new_end, pushing_max)

    #...................................................................................#
    def __get_new_Hlocal_start(self, bad_indices):

        # if the first local bin was basically empty, we determine
        # how many contigous bins from the first locabin were also empty
        #
        # the new_start is the INDEX into the first empty position in Hlocal
        # from the LHS (or -1 if the first position isn't empty)

        if bad_indices[0] == 0:
            index = 0
            for i in bad_indices:
                if i != index:
                    break
                else:
                    index = index + 1
            new_start = (index - 1)
            pushing_min = False
        else:
            new_start = -1
            pushing_min = True
        return (new_start, pushing_min)

    #...................................................................................#
    def mklog(self, logfile, initial=''):
        """
        Create a logfile
        """
        log = open(logfile, 'w')
        log.write(initial)
        log.close()
        return logfile

    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    #                        LOGGING FUNCTIONALITY
    #
    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #...................................................................................#
    def writeLog(self, logfile, output):
        """
            Write to a logfile
        """
        log = open(logfile, 'a')
        log.write(output)
        log.close()

    #...................................................................................#
    def fprintHVector(self, vector):
        """
           Convert a vector of integers into a string
           of integers seperated by a tab
        """

        out = ""
        for num in vector:
            out += "%d\t" % num
        return out

    #...................................................................................#
    def fprintGVector(self, vector):
        """
           Convert a vector of floats into a string
           of floats seperated by a tab
        """

        out = ""
        for num in vector:
            out += "%5.4f\t" % num
        return out

    #...................................................................................#
    def fprintVertVector(self, vector):
        out = ""
        for num in vector:
            out += "%5.4f\n" % num
        return out


# The code below is going into a stand alone WL executable. Ignore for now...


"""
## ===================================================================================================
##                              Main Script - hold onto your hat!
## ===================================================================================================

if __name__=="__main__":
    import argparse
    from keyfile import KeyFile as KF

    parser = argparse.ArgumentParser()

    parser.add_argument("--sequence-file","-s", help="File containing sequence protein sequence")
    parser.add_argument("--keyfile", "-k", help="Wang-Landau permutant generator keyfile")
    parser.add_argument("--test", help="Run Testfunction")

    #parser.add_argument("--freeze-file","-f", help="File containing regions of sequence which should not be moved during permutation")
    #parser.add_argument("--output-directory", "-o", help="Directory for logfiles to be written to")
    #parser.add_argument("--bins", "-", help="Which column is are numbers being extracted from [default=2]")
    #parser.add_argument("--target", "-t", help="Threshold value - will search for values, equal to, greater than or less than this value depending on the target-type flag")
    #parser.add_argument("--target-type","-tt", help="What kind of target are you searching for (gt, e, lt, ge, le, ne)")

    args = parser.parse_args()

    print "Done"

    if args.keyfile:
        keyfile_object=KF(args.keyfile)
        print ""
        print "Keyfile was parsed and validated succesfully"
        print ""

        bin_min_k             = keyfile_object.get_bin_min()
        bin_max_k             = keyfile_object.get_bin_max()
        nbins_k               = keyfile_object.get_num_bins()
        flatchk_k             = keyfile_object.get_flatcheck_freq()
        convergence_k         = keyfile_object.get_convergence()
        flatness_criterion_k  = keyfile_object.get_flatness_criterion()
        seq_k                 = keyfile_object.get_sequence()
        outdir_k              = keyfile_object.get_outdir()
        freezefile_k          = keyfile_object.get_freezefile()
        wl_type_k             = keyfile_object.get_WL_type()


        wl = WangLandauMachine(seq=seq_k,
                               writedir=outdir_k,
                               frzfile=freezefile_k,
                               nbins=nbins_k,
                               binmin=bin_min_k,
                               binmax = bin_max_k,
                               flatchk=flatchk_k,
                               flatcrit=flatness_criterion_k,
                               convergence = convergence_k,
                               WL_type = wl_type_k)


        wl.run()

    if args.test:
        testingMethod()

"""
