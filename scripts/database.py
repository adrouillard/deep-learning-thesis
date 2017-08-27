#! /usr/bin/env python
"""
Author: Antoine Drouillard
student nr : 930421200020
Date : 2017-03-13
This program will parse data from mass bank to only store useful data into a database.
"""
from __future__ import division
import mysql.connector
from sys import  argv
import re


def parse(data, origin):
    """
        OUTPUT : Return a matrix with peak values and translate a spectrum
        INPUT  : the txt file containing reads of molecules.
    """

    # Initialization
    formula = "NA" # store the formula of the molecule
    InChi = "NA"  # store the InChi code of the molecule
    mass = 0.0
    val = []
    annotation = []
    sciname = ''
    comname = ''
    idname = 0
    peaks = []
    spectra = False
    # Create and write label
    #print data, maxi
    l = 0
    for line in data:

        #print 'inside'
        if line.startswith('//'):
            l += 1
            #create the final string
            if len(annotation) != 0 and len(annotation) == len(val):
                for i in range(0, len(val)) :
                    peaks.append(val[i] + annotation[i])
            else :
                for i in range(0, len(val)) :
                    peaks.append(val[i] + [None,None,None,None])
            mol = (access,InChi, smiles, sciname, comname, nbpeak, origin, formula)
            imp_mol(mol)
            imp_peak(peaks, access)
            #print ids
            #print(ids.encode('string-escape'))
            #print val
            print l, "Molecule Saved", access
            #print round(l / 3723137, 3) * 100

            formula = "NA"  # store the formula of the molecule
            InChi = "NA"  # store the InChi code of the molecule
            mass = 0.0
            val = []
            annotation = []
            peaks = []
            spectra = False
            idname = 0

        elif line.startswith('CH$IUPAC'):
            if "=" in line:
                InChi = line.split("=")[1]
                InChi = InChi.strip("\n")
		if len(InChi)>500:
			InChi = ""
            #print InChi
        elif line.startswith('CH$FORMULA'):
            formula = line.split()[1]
        elif line.startswith('ACCESSION:'):
            access = line.split()[1]
        elif line.startswith('CH$EXACT_MASS'):
	    if ',' in line :
		line = line.replace(",",".")
            mass = line.split()[1]
        elif line.startswith('PK$PEAK'):
            spectra = True
        elif line.startswith('CH$SMILES'):
            smiles = line.split()[1]
        elif line.startswith(' ') and spectra == False:
	    if ',' in line :
		line = line.replace(",",".")
	    #print line.split()
            try :
		annotation.append([line.split()[1], float(line.split()[2]), \
                        float(line.split()[3]), float(line.split()[4])])
	    except :
		pass	
        elif line.startswith(' ') and spectra == True:
	    if ',' in line :
		line = line.replace(",",".")
            val.append(map(float,line.split()))
        elif line.startswith('CH$NAME'):
            if idname == 0 :
                comname = line.split()[1]
                idname += 1
            else:
                sciname = line.split()[1]
        elif line.startswith('PK$NUM_PEAK'):
            nbpeak = int(line.split()[1])
        #elif line.startswith('AUTHORS'):
        #    author = line.split(" ", 1)[1]
        #    author = author.strip("\n")
    print "done"


def imp_mol(alist):
    """
        OUTPUT : implement values in database
        INPUT  : A list of values to decribe the molecule to put inside the database
    """
    database = mysql.connector.connect(user='root', password='Macrobio21!', database='mass_spectra')

    cursor = database.cursor()

    #cursor.execute('SHOW DATABASES;')
    add_mol = """INSERT INTO molecule
              VALUES (%s,%s,%s,%s,%s,"%s",%s,%s)"""

    data_mol = (alist[0], "test", "test2", "test3", "test4", 0,
               "database", "C4")
    #print alist
    cursor.execute(add_mol, alist)
    #cursor.execute("""INSERT INTO molecule
    #          VALUES ("%s", "test", "test2", "test3", "test4", 0,"database", "C4")""" % alist[0])

    #emp_no = cursor.lastrowid

    database.commit()
    #print 'ok'
    cursor.close()
    database.close()

def imp_peak(alist, accession):
    """
        OUTPUT : implement values in database
        INPUT  : A list of values of peaks to put inside the database and the corresponding access
        for the molecule
    """
    database = mysql.connector.connect(user='root', password='Macrobio21!', database='mass_spectra')
    cursor = database.cursor()


    for i in range(0, len(alist)):
        # id = cursor.lastrowid
        # if id == None :
        #    id = 1
        # else :
        #    id = cursor.lastrowid
        alist[i] = alist[i] + [accession]
        #print tuple(alist[i])
        if None in alist[i] :
            add_peak = """INSERT INTO peaks
                           VALUES (NULL,"%s","%s",%s,%s,%s,%s,%s,%s)"""
        else :
            add_peak = """INSERT INTO peaks
                           VALUES (NULL,"%s","%s",%s,%s,"%s","%s","%s",%s)"""
        cursor.execute(add_peak, tuple(alist[i]))
        database.commit()


    cursor.close()
    database.close()

if __name__ == "__main__":
    data = open(argv[1], 'r')
    parse(data, "MassBank")
