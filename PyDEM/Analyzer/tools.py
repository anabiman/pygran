try:
	from xlrd import open_workbook
except:
	pass

from numbers import Number

def readExcel(fname):
	"""
	reads an excel sheet(s) and appends each data column to 
	a dictionary
	"""
	book = open_workbook(fname, on_demand=True)
	data = {}

	for name in book.sheet_names():
		sheet = book.sheet_by_name(name)

		for coln in range(sheet.ncols):
			for celln, cell in enumerate(sheet.col(coln)):
				if celln == 0:
					dname = cell.value
					data[dname] = []
				else:
					if isinstance(cell.value, Number):
						data[dname].append(cell.value)
					else:
						print cell.value, ' ignored'

	for key in data.keys():
		data[key] = array(data[key])

	return data