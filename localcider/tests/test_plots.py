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
   ! MAIN AUTHOR:   Alex Holehouse                                            !
   !                                                                          !
   !--------------------------------------------------------------------------!


   Unit tests for functions in plots


"""

import os.path
import time
import unittest
import random
from localcider import plots, sequenceParameters
import testTools


class TestPlotsFunctions(unittest.TestCase):

    def setUp(self):
        self.rseq = sequenceParameters.SequenceParameters(
            testTools.generate_random_sequence(minLen=20, maxLen=500))

    def test_save_single_phasePlot(self):
        fp = self.rseq.get_fraction_positive()
        fn = self.rseq.get_fraction_negative()

        plots.save_single_phasePlot(fp, fn, 'tmpfiles/single_PP.png', 'TEST TITLE')

        plots.save_single_phasePlot(
            fp, fn, 'tmpfiles/single_PP_NO_LEGEND.png', 'TEST TITLE', False)

        plots.save_single_phasePlot(
            fp, fn, 'tmpfiles/single_PP_NO_LEGEND.pdf', 'TEST TITLE', False, saveFormat='pdf')

    def test_save_multiple_phasePlot(self):

        rseq2 = sequenceParameters.SequenceParameters(
            testTools.generate_random_sequence(minLen=20, maxLen=500))

        fp1 = self.rseq.get_fraction_positive()
        fn1 = self.rseq.get_fraction_negative()
        fp2 = rseq2.get_fraction_positive()
        fn2 = rseq2.get_fraction_negative()

        plots.save_multiple_phasePlot(
            [fp1, fp2], [fn1, fn2], 'tmpfiles/mult_PP', ['a', 'b'], 'TEST TITLE')
        plots.save_multiple_phasePlot([fp1,
                                       fp2],
                                      [fn1,
                                       fn2],
                                      'tmpfiles/mult_PP_NO_LEGEND.png',
                                      ['a',
                                       'b'],
                                      title='TEST TITLE',
                                      legendOn=False)

        plots.save_multiple_phasePlot([fp1,
                                       fp2],
                                      [fn1,
                                       fn2],
                                      'tmpfiles/mult_PP_NO_LEGEND_xLim1.png',
                                      ['a',
                                       'b'],
                                      title='TEST TITLE',
                                      legendOn=False,
                                      xLim=1)

        plots.save_multiple_phasePlot([fp1,
                                       fp2],
                                      [fn1,
                                       fn2],
                                      'tmpfiles/mult_PP_NO_LEGEND_xlim2.png',
                                      ['a',
                                       'b'],
                                      title='TEST TITLE',
                                      legendOn=False,
                                      xLim=2)

        plots.save_multiple_phasePlot([fp1,
                                       fp2],
                                      [fn1,
                                       fn2],
                                      'tmpfiles/mult_PP_NO_LEGEND_xlim2_ylim2.png',
                                      ['a',
                                       'b'],
                                      title='TEST TITLE',
                                      legendOn=False,
                                      xLim=2,
                                      yLim=2)

        plots.save_multiple_phasePlot([fp1,
                                       fp2],
                                      [fn1,
                                       fn2],
                                      'tmpfiles/mult_PP_NO_LEGEND_xlim2_ylim2_fs_20.png',
                                      ['a',
                                       'b'],
                                      title='TEST TITLE',
                                      legendOn=False,
                                      xLim=2,
                                      yLim=2,
                                      fontSize=20)

        plots.save_multiple_phasePlot([fp1,
                                       fp2],
                                      [fn1,
                                       fn2],
                                      'tmpfiles/mult_PP_NO_LEGEND_xlim2_ylim2_fs_20.pdf',
                                      ['a',
                                       'b'],
                                      title='TEST TITLE',
                                      legendOn=False,
                                      xLim=2,
                                      yLim=2,
                                      fontSize=20,
                                      saveFormat='pdf')

    def test_save_multiple_phasePlot2(self):

        rseq2 = sequenceParameters.SequenceParameters(
            testTools.generate_random_sequence(minLen=20, maxLen=500))

        plots.save_multiple_phasePlot2(
            [self.rseq, rseq2], 'tmpfiles/mult_PP', ['a', 'b'], 'TEST TITLE')
        plots.save_multiple_phasePlot2(
            [self.rseq, rseq2], 'tmpfiles/mult_PP_NO_LEGEND', ['a', 'b'], 'TEST TITLE', False)
        plots.save_multiple_phasePlot2([self.rseq,
                                        rseq2],
                                       'tmpfiles/mult_PP_NO_LEGEND_xLim1.png',
                                       ['a',
                                        'b'],
                                       title='TEST TITLE',
                                       legendOn=False,
                                       xLim=1)

        plots.save_multiple_phasePlot2([self.rseq,
                                        rseq2],
                                       'tmpfiles/mult_PP_NO_LEGEND_xlim2.png',
                                       ['a',
                                        'b'],
                                       title='TEST TITLE',
                                       legendOn=False,
                                       xLim=2)

        plots.save_multiple_phasePlot2([self.rseq,
                                        rseq2],
                                       'tmpfiles/mult_PP_NO_LEGEND_xlim2_ylim2.png',
                                       ['a',
                                        'b'],
                                       title='TEST TITLE',
                                       legendOn=False,
                                       xLim=2,
                                       yLim=2)

        plots.save_multiple_phasePlot2([self.rseq,
                                        rseq2],
                                       'tmpfiles/mult_PP_NO_LEGEND_xlim2_ylim2_fs_20.png',
                                       ['a',
                                        'b'],
                                       title='TEST TITLE',
                                       legendOn=False,
                                       xLim=2,
                                       yLim=2,
                                       fontSize=20)

        plots.save_multiple_phasePlot2([self.rseq,
                                        rseq2],
                                       'tmpfiles/mult_PP_NO_LEGEND_xlim2_ylim2_fs_20.pdf',
                                       ['a',
                                        'b'],
                                       title='TEST TITLE',
                                       legendOn=False,
                                       xLim=2,
                                       yLim=2,
                                       fontSize=20,
                                       saveFormat='pdf')

    def test_save_single_uverskyPlot(self):
        hydro = self.rseq.get_uversky_hydropathy()
        mnc = self.rseq.get_mean_net_charge()

        plots.save_single_uverskyPlot(
            hydro, mnc, 'tmpfiles/single_UV', 'TEST TITLE')

        plots.save_single_uverskyPlot(
            hydro,
            mnc,
            'tmpfiles/single_UV_NO_LEGEND.png',
            'TEST TITLE',
            False)

        plots.save_single_uverskyPlot(
            hydro,
            mnc,
            'tmpfiles/single_UV_NO_LEGEND.pdf',
            'TEST TITLE',
            False,
            saveFormat='pdf')

    def test_save_multiple_uverskyPlot(self):

        rseq2 = sequenceParameters.SequenceParameters(
            testTools.generate_random_sequence(minLen=20, maxLen=500))

        hy1 = self.rseq.get_uversky_hydropathy()
        mnc1 = self.rseq.get_mean_net_charge()
        hy2 = rseq2.get_uversky_hydropathy()
        mnc2 = rseq2.get_mean_net_charge()

        plots.save_multiple_uverskyPlot(
            [hy1, hy2], [mnc1, mnc2], 'tmpfiles/mult_UV', ['a', 'b'], 'TEST TITLE')
        plots.save_multiple_uverskyPlot([hy1,
                                         hy2],
                                        [mnc1,
                                         mnc2],
                                        'tmpfiles/mult_UV_NO_LEGEND.png',
                                        ['a',
                                         'b'],
                                        title='TEST TITLE',
                                        legendOn=False)

        plots.save_multiple_uverskyPlot([hy1,
                                         hy2],
                                        [mnc1,
                                         mnc2],
                                        'tmpfiles/mult_UV_NO_LEGEND_xLim1.png',
                                        ['a',
                                         'b'],
                                        title='TEST TITLE',
                                        legendOn=False,
                                        xLim=1)

        plots.save_multiple_uverskyPlot([hy1,
                                         hy2],
                                        [mnc1,
                                         mnc2],
                                        'tmpfiles/mult_UV_NO_LEGEND_xlim2.png',
                                        ['a',
                                         'b'],
                                        title='TEST TITLE',
                                        legendOn=False,
                                        xLim=2)

        plots.save_multiple_uverskyPlot([hy1,
                                         hy2],
                                        [mnc1,
                                         mnc2],
                                        'tmpfiles/mult_UV_NO_LEGEND_xlim2_ylim2',
                                        ['a',
                                         'b'],
                                        title='TEST TITLE',
                                        legendOn=False,
                                        xLim=2,
                                        yLim=2)

        plots.save_multiple_uverskyPlot([hy1,
                                         hy2],
                                        [mnc1,
                                         mnc2],
                                        'tmpfiles/mult_UV_NO_LEGEND_xlim2_ylim2_fs_20.png',
                                        ['a',
                                         'b'],
                                        title='TEST TITLE',
                                        legendOn=False,
                                        xLim=2,
                                        yLim=2,
                                        fontSize=20)

        plots.save_multiple_uverskyPlot([hy1,
                                         hy2],
                                        [mnc1,
                                         mnc2],
                                        'tmpfiles/mult_UV_NO_LEGEND_xlim2_ylim2_fs_20.pdf',
                                        ['a',
                                         'b'],
                                        title='TEST TITLE',
                                        legendOn=False,
                                        xLim=2,
                                        yLim=2,
                                        fontSize=20,
                                        saveFormat='pdf')

    def test_save_multiple_uverskyPlot2(self):

        rseq2 = sequenceParameters.SequenceParameters(
            testTools.generate_random_sequence(minLen=20, maxLen=500))

        hy1 = self.rseq.get_uversky_hydropathy()
        mnc1 = self.rseq.get_mean_net_charge()
        hy2 = rseq2.get_uversky_hydropathy()
        mnc2 = rseq2.get_mean_net_charge()

        plots.save_multiple_uverskyPlot2(
            [self.rseq, rseq2], 'tmpfiles/mult_UV.png', ['a', 'b'], 'TEST TITLE')
        
        plots.save_multiple_uverskyPlot2([self.rseq,
                                          rseq2],
                                         'tmpfiles/mult_UV2_NO_LEGEND.png',
                                         ['a',
                                          'b'],
                                         title='TEST TITLE',
                                         legendOn=False)
        
        plots.save_multiple_uverskyPlot2([self.rseq,
                                          rseq2],
                                         'tmpfiles/mult_UV2_NO_LEGEND_xLim1.png',
                                         ['a',
                                          'b'],
                                         title='TEST TITLE',
                                         legendOn=False,
                                         xLim=1)

        plots.save_multiple_uverskyPlot2([self.rseq,
                                          rseq2],
                                         'tmpfiles/mult_UV2_NO_LEGEND_xlim2.png',
                                         ['a',
                                          'b'],
                                         title='TEST TITLE',
                                         legendOn=False,
                                         xLim=2)

        plots.save_multiple_uverskyPlot2([self.rseq,
                                          rseq2],
                                         'tmpfiles/mult_UV2_NO_LEGEND_xlim2_ylim2.png',
                                         ['a',
                                          'b'],
                                         title='TEST TITLE',
                                         legendOn=False,
                                         xLim=2,
                                         yLim=2)

        plots.save_multiple_uverskyPlot2([self.rseq,
                                          rseq2],
                                         'tmpfiles/mult_UV2_NO_LEGEND_xlim2_ylim2_fs_20.png',
                                         ['a',
                                          'b'],
                                         title='TEST TITLE',
                                         legendOn=False,
                                         xLim=2,
                                         yLim=2,
                                         fontSize=20)

        plots.save_multiple_uverskyPlot2([self.rseq,
                                          rseq2],
                                         'tmpfiles/mult_UV2_NO_LEGEND_xlim2_ylim2_fs_20.pdf',
                                         ['a',
                                          'b'],
                                         title='TEST TITLE',
                                         legendOn=False,
                                         xLim=2,
                                         yLim=2,
                                         fontSize=20,
                                         saveFormat='pdf')
