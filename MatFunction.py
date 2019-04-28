import numpy as np
def finddt(dt,dIndex):
	# dt order was ascending
	# dtcmp = dt
	length = len(dt)
	midlength = int(length/2)
	# print(dt)
	# print(dIndex)
	# print('dt[length] = '+str(dt[length-1]))
	# print('dt[0] = '+str(dt[0]))
	if dIndex>dt[length-1] or dIndex<dt[0]:
		print('error dIndex>max or dIndex<min')
		pass
	if length<3:
		# print(dt[0])
		# print(dt[1])
		return dt
	else:
		if dIndex>dt[midlength] or dIndex==dt[midlength]:
			# print(dt[midlength:length-1])
				# return dt[midlength]
				# pass
				return finddt(dt[midlength:length],dIndex)
		elif dIndex<dt[midlength]:
			# print(dt[0:midlength])

				return finddt(dt[0:midlength+1],dIndex)
		else:
			pass

def findIndex(lt,dIndex):
	i = finddt(lt,dIndex)
	index = np.argwhere(lt == i[1])[0][0]
	return index


def test():
	testList = np.arange(1,12,0.1)
	testIndex = 3.51
	i = findIndex(testList,testIndex)
	print(i)

