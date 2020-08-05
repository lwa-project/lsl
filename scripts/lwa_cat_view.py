#! /usr/bin/env python

"""
Simple LWA1 astronomical source catalogue display application.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import math
import argparse
import Tkinter

from lsl import astro
from lsl import transform
from lsl import catalog
from lsl.common import stations

from lsl.misc import telemetry
telemetry.track_script()


__version__   = "0.1"
__author__	= "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class CatalogViewer(object):
    
    def __init__(self, root, args):
        
        self.root = root
        self.period = args.period * 1000
        self.catalog = None
        self.cursor = self.root.cget("cursor")
        site = args.site.lower().replace('-', '')
        if site == 'lwa1':
            station = stations.lwa1
        elif site == 'lwasv':
            station = stations.lwasv
        elif site == 'ovrolwa':
            station = stations.lwa1
            station.lat, station.lon, station.elev = ('37.2397808', '-118.2816819', 1183.4839)
        else:
            raise RuntimeError("Unknown site name: %s" % site)
        self.site = transform.GeographicalPosition((station.long*180.0/math.pi, station.lat*180.0/math.pi), name=station.name)
        
        self.root.title('LWA Catalog Viewer')
        
        frame = Tkinter.Frame(self.root)
        frame.pack(side = Tkinter.TOP, fill = Tkinter.BOTH, expand = True)
        
        self.catalog_name = Tkinter.StringVar()
        frame1 = Tkinter.Frame(frame)
        frame1.pack(side = Tkinter.LEFT, fill = Tkinter.Y)
        catNames = catalog.CatalogFactory.get_names()
        catNames.sort()
        for name in catNames:
            button = Tkinter.Radiobutton(frame1, text = name, value = name, anchor = Tkinter.W,
                variable = self.catalog_name, command = self.select_catalog, indicatoron = False)
            button.pack(fill = Tkinter.X)
        button = Tkinter.Radiobutton(frame1, text = 'PLANETS', value = 'PLANETS', anchor = Tkinter.W,
            variable = self.catalog_name, command = self.select_catalog, indicatoron = False)
        button.pack(fill = Tkinter.X)
        button = Tkinter.Button(frame1, text = 'Select', command = self.select_source)
        button.pack(fill = Tkinter.X, side = Tkinter.BOTTOM)
        
        
        frame2 = Tkinter.Frame(frame)
        frame2.pack(side = Tkinter.RIGHT, fill = Tkinter.BOTH, expand = True)
        yscroll = Tkinter.Scrollbar(frame2, orient = Tkinter.VERTICAL)
        xscroll = Tkinter.Scrollbar(frame2, orient = Tkinter.HORIZONTAL)
        self.source_list = Tkinter.Listbox(frame2, height = 40, width = 60,
            yscrollcommand = yscroll.set, xscrollcommand = xscroll.set, 
            font = "fixed 12 bold", selectmode = Tkinter.SINGLE, exportselection = False,
            selectbackground = 'yellow')
        yscroll.config(command = self.source_list.yview)
        xscroll.config(command = self.source_list.xview)
        self.source_list.bind("<Double-Button-1>", self.select_source)
        yscroll.pack(side = Tkinter.RIGHT, fill = Tkinter.Y)
        xscroll.pack(side = Tkinter.BOTTOM, fill = Tkinter.X)
        self.source_list.pack(side = Tkinter.LEFT, fill = Tkinter.BOTH, expand = True)
        
        frame = Tkinter.Frame(self.root)
        frame.pack(side = Tkinter.BOTTOM, fill = Tkinter.BOTH)
        
        self.utc_time = Tkinter.Label(frame, relief = Tkinter.SUNKEN,
            anchor = Tkinter.W)
        self.last_time = Tkinter.Label(frame, relief = Tkinter.SUNKEN,
            anchor = Tkinter.W)
        self.utc_time.pack(side = Tkinter.LEFT, fill = Tkinter.BOTH, expand = True)
        self.last_time.pack(side = Tkinter.RIGHT, fill = Tkinter.BOTH, expand = True)
        
        self.root.after(0, self.timer)
        
    def timer(self):
    
        currentTime = transform.Time.from_system()
        self.utc_time.config(text = currentTime.utc_str)
        
        last = astro.deg_to_dms(self.site.sidereal(currentTime))
        self.last_time.config(text = "%s LAST: %s" % (self.site.name, last)) 
        
        self.root.after(self.period, self.timer)
        
    def open_catalog(self):
        
        self.source_list.delete(0, Tkinter.END)
        
        if self.catalog_name.get() == 'PLANETS':
        
            self.select_planet()
            self.root.config(cursor = self.cursor)
            return
            
        self.catalog = catalog.CatalogFactory.get_catalog(self.catalog_name.get())
        
        for name in sorted(self.catalog.keys()):
        
            source = self.catalog[name]
            name = source.name.strip().ljust(18)
            s = "%s | " % name
            equ = source.position.j2000_equ
            sra = "%.03f" % equ.ra
            sdec = "%.03f" % equ.dec
            if equ.dec >= 0:
                sdec = '+' + sdec
            equ = equ.format()
            sra = str(equ[0]) + ' (' + sra.zfill(7) + ')'
            sdec = str(equ[1]) + ' (' + sdec.zfill(7) + ') '
            s += "%s | %s " % (sra.rjust(20), sdec.rjust(20))
            for alias in sorted(source.alias_list):
                s += "%s " % alias.strip() 
            self.source_list.insert(Tkinter.END, s)
            
        self.root.config(cursor = self.cursor)
        
    def select_catalog(self):
    
        self.root.config(cursor = "watch")
        self.root.after(0, self.open_catalog)
            
    def select_planet(self):
    
        self.catalog = None
    
        planets = sorted(transform.PlanetaryPosition.known_names)
        for name in planets:
            self.source_list.insert(Tkinter.END, name)
            
    def select_source(self, *args):
        
        select = self.source_list.curselection()
        if not len(select):
            return
        select = self.source_list.get(select[0])
        name = select.split()[0]
        
        SourceWindow(self.catalog, name, self.site, self.period)


class SourceWindow(Tkinter.Toplevel):
    
    def __init__(self, cat, name, site, period):
        
        Tkinter.Toplevel.__init__(self)
        
        self.lift()
        self.catalog = cat
        self.site = site
        self.period = period
        self.name = name
        
        if self.catalog is None:
            catName = 'PLANETS'
        else:
            catName = self.catalog.name
        self.title("%s/%s" % (catName, name))
        
        sourceFrame = Tkinter.Frame(self)
        sourceFrame.pack(expand = True, fill = Tkinter.X)
        
        labelFrame = Tkinter.Frame(sourceFrame)
        labelFrame.pack(side = Tkinter.LEFT, fill = Tkinter.X)
        
        SourceLabel(labelFrame, 'Visible')
        SourceLabel(labelFrame, 'J2000 RA')
        SourceLabel(labelFrame, 'J2000 DEC')
        SourceLabel(labelFrame, 'Current RA')
        SourceLabel(labelFrame, 'Current DEC')
        SourceLabel(labelFrame, 'Galactic Long.')
        SourceLabel(labelFrame, 'Galactic Lat.')
        SourceLabel(labelFrame, 'Azimuth')
        SourceLabel(labelFrame, 'Altitude')
        SourceLabel(labelFrame, 'Current Time')
        SourceLabel(labelFrame, 'Rise Time')
        SourceLabel(labelFrame, 'Transit Time')
        SourceLabel(labelFrame, 'Set Time')
        
        self.values = {}
        
        valueFrame = Tkinter.Frame(sourceFrame)
        valueFrame.pack(side = Tkinter.RIGHT, fill = Tkinter.X, expand = True)
        
        self.values['visible'] = SourceValue(valueFrame)
        self.values['visible'].config(anchor = Tkinter.CENTER, justify = Tkinter.CENTER)
        self.values['mean_ra'] = SourceValue(valueFrame)
        self.values['mean_dec'] = SourceValue(valueFrame)
        self.values['cur_ra'] = SourceValue(valueFrame)
        self.values['cur_dec'] = SourceValue(valueFrame)
        self.values['gal_lng'] = SourceValue(valueFrame)
        self.values['gal_lat'] = SourceValue(valueFrame)
        self.values['az'] = SourceValue(valueFrame)
        self.values['alt'] = SourceValue(valueFrame)
        self.values['date'] = SourceValue(valueFrame)
        self.values['rise'] = SourceValue(valueFrame)
        self.values['transit'] = SourceValue(valueFrame)
        self.values['set'] = SourceValue(valueFrame)		
        
        if self.catalog is None:
            self.position = transform.PlanetaryPosition(self.name)
        else:
            source = self.catalog[self.name]
            self.position = source.position
        self.direction = transform.PointingDirection(self.position, self.site)
        
        self.after(0, self.timer)
        
    def timer(self):
        
        currentTime = transform.Time.from_system()
        
        cur_equ = self.position.apparent_equ(currentTime)
        if self.catalog is None:
            mean_equ = astro.get_equ_prec2(cur_equ, currentTime.utc_jd, astro.J2000_UTC_JD)
        else:
            mean_equ = self.position.j2000_equ
            
        (mean_ra, mean_dec) = mean_equ.format()
        (cur_ra, cur_dec) = cur_equ.format()
        
        if self.catalog is None:
            gal = astro.get_gal_from_equ2000(mean_equ)
        else:
            gal = self.position.j2000_gal
        (gal_lng, gal_lat) = gal.format()
        
        hrz = self.direction.hrz(currentTime)
        rst = self.direction.rst(currentTime)
        if rst is None:
            rise = None
            transit = None
            set = None
        else:
            (rise, set, transit) = rst.format()
            
        if hrz.alt > 10.0:
            self.values['visible'].config(text = "YES", bg = "green")
        else:
            self.values['visible'].config(text = "NO", bg = "red")
            
        self.values['mean_ra'].config(text = mean_ra)
        self.values['mean_dec'].config(text = mean_dec)
        self.values['cur_ra'].config(text = cur_ra)
        self.values['cur_dec'].config(text = cur_dec)
        self.values['gal_lng'].config(text = gal_lng)
        self.values['gal_lat'].config(text = gal_lat)
        self.values['az'].config(text = "%0.3f %s" % (hrz.az, astro.hrz_to_nswe(hrz)))
        self.values['alt'].config(text = "%0.3f" % hrz.alt)
        self.values['date'].config(text = currentTime.utc_str)
        self.values['rise'].config(text = rise)
        self.values['transit'].config(text = transit)
        self.values['set'].config(text = set)
        
        self.after(self.period, self.timer)


class SourceLabel(Tkinter.Label):
    
    def __init__(self, frame, text):
        
        Tkinter.Label.__init__(self, frame, text = text, anchor = Tkinter.W,
            justify = Tkinter.LEFT, relief = Tkinter.GROOVE, font = "fixed 12 bold")
        self.pack(fill = Tkinter.X)


class SourceValue(Tkinter.Label):
    
    def __init__(self, frame):
        
        Tkinter.Label.__init__(self, frame, anchor = Tkinter.E, font = "fixed 12 bold",
            justify = Tkinter.RIGHT, relief = Tkinter.SUNKEN, foreground = 'blue')
        self.pack(fill = Tkinter.X)


if __name__ == '__main__':
    # parse command line
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('-s', '--site', type=str, default='LWA-1',
                        help='site name')
    parser.add_argument('-p', '--period', type=int, default=5, 
                        help='update period in seconds')
    args = parser.parse_args()
    
    # create the GUI
    tk = Tkinter.Tk()
    CatalogViewer(tk, args)
    
    # execute the GUI
    try:
        tk.mainloop()
    except KeyboardInterrupt:
        pass
        
