""" 
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of LocalCider.                                      !
   !                                                                          !
   !    Version 0.1.0                                                           !
   !                                                                          !
   !    Copyright (C) 2014, The LocalCIDER development team (current and      !
   !                        former contributors): Alex Holehouse, James       !
   !                        Ahad, Rahul K. Das.                               !
   !                                                                          !
   !    LocalCIDER is a PappuLab tool. Please see the website for citation    !
   !    information.                                                          !
   !                                                                          !
   !    Website: http://pappulab.github.io/LocalCIDER/                        !
   !                                                                          !
   !    LocalCIDER is free software: you can redistribute it and/or modify    !
   !    it under the terms of the GNU General Public License as published by  !
   !    the Free Software Foundation, either version 3 of the License, or     !
   !    (at your option) any later version.                                   !
   !                                                                          !
   !    LocalCIDER is distributed in the hope that it will be useful,         !
   !    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
   !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
   !    GNU General Public License for more details.                          !
   !                                                                          !
   !    You should have received a copy of the GNU General Public License     !
   !    along with LocalCIDER.  If not, see <http://www.gnu.org/licenses/>.   !
   !--------------------------------------------------------------------------!
   ! AUTHORSHIP INFO:                                                         !
   !--------------------------------------------------------------------------!
   !                                                                          !
   ! MAIN AUTHOR:   Alex Holehouse                                            !
   !                                                                          !
   !--------------------------------------------------------------------------!

   
   File Description:
   ================
   
   The WL files contain code for running Wang-Landau Monte Carlo. Specifically,
   this codes allows you to;

   - Calculate the density of states using standard Wang-Landau ("Normal" Wang-Landau)

   - Employ the HistogramZoom algorithm to rapidly estimate where, for a sequence, the
     primary regions of high state density lie (for a complete description of 
     Histogram please see the localCIDER doumentation 

   - Generate kappa permutants of a range of kappa values for a specific sequence


"""

import wl
