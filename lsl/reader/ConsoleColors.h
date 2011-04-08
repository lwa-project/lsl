/*
  MCS-DROS: Monitor and Control System - Data Recorder Operating Software
  Copyright (C) 2009-2010  Virginia Tech, Christopher Wolfe <chwolfe2@vt.edu>

  This file is part of MCS-DROS.
 
  MCS-DROS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  MCS-DROS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with MCS-DROS.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
 * ConsoleColors.h
 *
 *  Created on: Dec 18, 2009
 *      Author: chwolfe2
 */

#ifndef CONSOLECOLORS_H_
#define CONSOLECOLORS_H_
#define USE_ANSI_COLORS
#ifdef USE_ANSI_COLORS
	#define COLOR_NORMAL "\033[1;37m"
	#define COLOR_RED "\033[0;31m"
	#define COLOR_GREEN "\033[0;32m"
	#define COLOR_BLUE "\033[0;34m"
	#define COLOR_GREY "\033[1;30m"
	#define HLG(a) COLOR_GREEN a COLOR_NORMAL
	#define HLR(a) COLOR_RED a COLOR_NORMAL
	#define HLB(a) COLOR_BLUE a COLOR_NORMAL
#else
	#define COLOR_NORMAL ""
	#define COLOR_RED ""
	#define COLOR_GREEN ""
	#define COLOR_BLUE ""
	#define COLOR_GREY ""
	#define HLG(a) a
	#define HLR(a) a
	#define HLB(a) a
#endif


#endif /* CONSOLECOLORS_H_ */
