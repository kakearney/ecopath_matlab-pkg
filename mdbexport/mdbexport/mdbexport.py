"""
.mdb export module

Kelly Kearney
"""

import pyodbc
import csv
import sys
import os

def mdbexport(outdir, mdbfile):
	"""
	Export tables from a MS Access Database (.mdb) to csv files

	Input arguments:
	outdir:		name of folder where .csv files will be saved.  
				If it doesn't exisit, it will be created
	mdbfile:	name of database file

	The resulting files will be named according to the table
	names in the database.

	"""

	# Create output folder if necessary

	if not os.path.isdir(outdir):
		os.mkdir(outdir)

    # Look for mdb driver
    
    DRV = ['{{{}}}'.format(s) for s in pyodbc.drivers() if 'Microsoft Access Driver' in s]
    if not DRV:
        raise ValueError('mdbexport did not find a valid Microsoft Access Driver on this computer')

	# Connect to database

	PWD = 'pw'

	con = pyodbc.connect('DRIVER={};DBQ={};PWD={};CHARSET=UTF8'.format(drv,mdbfile,PWD))            
	cur = con.cursor()

	# List tables in file, filtering out system tables

	tables = list(cur.tables())

	tablenames = [x[2] for x in tables if x[3]=='TABLE']

	# For now, filter out tables with spaces in the name
	# because I have no idea how to query those

	problemnames = [x for x in tablenames if ' ' in x]
	tablenames = list(set(tablenames) - set(problemnames))

	# Dump tables to csv

	for tbl in tablenames:

		# Column names
		
		cols = [colrow.column_name for colrow in cur.columns(table=tbl)]
		
		# Main data
		
		SQL = 'SELECT * FROM {}'.format(tbl)
		rows = cur.execute(SQL).fetchall()
		
		# Write to file
		
		csvname = os.path.join(outdir, tbl + '.csv')
		with open(csvname, 'w') as fou:
			c = csv.writer(fou, lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC)
			c.writerow(cols)
			c.writerows(rows)

	cur.close()
	con.close()
