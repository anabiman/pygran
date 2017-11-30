#!/usr/bin/python
# -*- coding: utf8 -*- 

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# -------------------------------------------------------------------------


'''
Created on Sep 8, 2016
@author: Andrew Abi-Mansour
'''
'''
Created on Sep 11, 2013

@author: Andrew Abi-Mansour
'''

__SofName__ = 'PyGran GUI'

class French:
    def __init__(self):
        __SofName__ = __SofName__
        self.File = '&Dossier'
        self.Open = '&Ouvrir'
        self.Close = '&Terminer'
        self.Configure = '&Configurer'
        self.Analysis = '&Analyse'
        self.Help = '&Assistance'
        self.ChangeDir = '&Changer Dir'
        
        self.description = """Analyse macromoléculaire et la modélisation est une boîte à 
        outils la construction et la manipulation des capsides virales. Les caractéristiques comprennent puissant haut-
générateur de capside d'une unité monomère, des capacités avancées de recherche sur le Web, la production oligomères
et la manipulation, l'intégration avec un puissant logiciel de MD de visualisation comme VMD, Chimera,
etc, etc."""


class Arabic:
    def __init__(self):
        __SofName__ = '<<كرامر>>'
        self.File = '&ملف'
        self.Open = '&افتح'
        self.Close = '&أغلق'
        self.Configure = '&تكوين'
        self.Analysis = '&تحليل'
        self.Help = '&مساعدة'
        self.ChangeDir = 'غير دليل'
        
        self.description = """ تحليل الجزيئات والنمذجة هو مجموعة أدوات لبناء والتلاعب الفيروسية وتشمل الميزات 
        القوية المدمج في مولد قفيصة من وحدة مونومر، قدرات متطورة لتصفح الانترنت جيل قليل 
        وحدات والتلاعب، والتكامل مع قوة التصور البرمجيات مثل ,  وغيرها، وأكثر من 
        .ذلك
        """
        
class English:
    def __init__(self):
        self.File = '&File'
        self.Open = '&Open'
        self.Close = '&Close'
        self.Configure = '&Configure' 
        self.Analysis = '&Analysis'
        self.Help = '&Help'
        self.ChaneDir ='Change &Dir'
        
        self.description = """
        DEM simulation engine, released by 
        DCS Computing Gmbh, Linz, Austria
        www.dcs-computing.com, office@dcs-computing.com

        LIGGGHTS® is open-source, distributed under the 
        terms of the GNU Public License, version 2 or later.

        LIGGGHTS® is part of CFDEM®project: 
        www.liggghts.com | www.cfdem.com

        Core developer and main author:
        Andrew Abi-Mansour
        """
