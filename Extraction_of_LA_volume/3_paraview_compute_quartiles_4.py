import os
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

def compute_write_Quartiles(vtu_path, file, avg_list, tstep_starting_avg, counter):
	# create a new 'XML Unstructured Grid Reader'
	c_laa_ = XMLUnstructuredGridReader(registrationName='C_laa_', FileName=[vtu_path])
	# Properties on c_laa_
	array_name = "nu_relative"
	c_laa_.PointArrayStatus = ['{}'.format(str(array_name))]

	# get active view
	renderView1 = GetActiveViewOrCreate('RenderView')
	# show data in view
	c_laa_Display = Show(c_laa_, renderView1, 'UnstructuredGridRepresentation')
	# update the view to ensure updated data information
	renderView1.Update()
	# set active view
	SetActiveView(renderView1)
	# set active source
	SetActiveSource(c_laa_)

	# get animation scene
	animationScene1 = GetAnimationScene()
	# update animation scene based on data timesteps
	animationScene1.UpdateAnimationUsingDataTimeSteps()

	#tsteps= vtu_path.split("/")[-1].split("_")[-1]
	tsteps = c_laa_.TimestepValues # List with the different time steps
	print("tsteps= ", tsteps)

	c_laa_.UpdatePipeline()

	#----------------------------------------------------------------------------------------
	# Compute quartiles
	computeQuartiles1 = ComputeQuartiles(registrationName='ComputeQuartiles1', Input=c_laa_)
	computeQuartiles1.UpdatePipeline()
	#----------------------------------------------------------------------------------------
	# Compute mean value, creating a new 'Python Annotation'
	pythonAnnotation1 = PythonAnnotation(registrationName='PythonAnnotation1', Input=c_laa_)
	# Properties modified on pythonAnnotation1
	pythonAnnotation1.ArrayAssociation = 'Point Data'
	pythonAnnotation1.Expression = 'mean({})'.format(str(array_name))
	pythonAnnotation1.UpdatePipeline()
	#----------------------------------------------------------------------------------------

	# Fetch the data from the the vtkColumn -> compute mean value of the array_name
	column=paraview.servermanager.Fetch(pythonAnnotation1)
	c=column.GetColumn(0)
	mean_value = c.GetValue(0) # mean_value is the type -> "string"

	#----------------------------------------------------------------------------------------
	# Fetch the data from the the vtkTable -> Get the 5 statistical data
	p=paraview.servermanager.Fetch(computeQuartiles1)
	# It is saved as Row Data
	row=p.GetRowData()
	# It is the only one array that we are looking for
	arr = row.GetArray(0)

	for i in range(arr.GetNumberOfTuples()):
		print(arr.GetTuple(i)[0])
		if tsteps[0] >= tstep_starting_avg:
			avg_list[i] = avg_list[i] + arr.GetTuple(i)[0]
	print("mean_value = {}".format(str(mean_value)), "\n")

	if tsteps[0] >= tstep_starting_avg:
		avg_list[5] = avg_list[5] + float(mean_value)
		counter = counter + 1
	#----------------------------------------------------------------------------------------
	file.write( str(tsteps[0])+"\t"+str(arr.GetTuple(0)[0])+"\t"+str(arr.GetTuple(1)[0])+"\t"+str(arr.GetTuple(2)[0])+"\t"+str(arr.GetTuple(3)[0])+"\t"+str(arr.GetTuple(4)[0])+"\t"+str(mean_value)+"\n" )

	return avg_list, counter


if __name__ == "__main__":
	N_tsteps = 200
	#N_tsteps = 5
	t_init=0

	file_LAA = open('/home/sergionio/Downloads/test_nu/200_tsteps/folder_quartiles_new/quartiles_LAA.txt', "a")
	file_LAA.write("#t (ms)"+"\t\t"+"min"+"\t\t"+"25%"+"\t\t"+"50%"+"\t\t"+"75%"+"\t\t"+"max"+"\t\t"+"mean_value"+"\n")
	# Save the average: min, 1st_quartile, 2nd_quartile, 3rd_quartile, max, mean_value
	avg_LAA = [0, 0, 0, 0, 0, 0]
	tstep_starting_avg_LAA = 990
	counter_LAA = 0

	file_LA_body = open('/home/sergionio/Downloads/test_nu/200_tsteps/folder_quartiles_new/quartiles_LA_body.txt', "a")
	file_LA_body.write("#t (ms)"+"\t\t"+"min"+"\t\t"+"25%"+"\t\t"+"50%"+"\t\t"+"75%"+"\t\t"+"max"+"\t\t"+"mean_value"+"\n")
	# Save the average: min, 1st_quartile, 2nd_quartile, 3rd_quartile, max, mean_value
	avg_LA_body = [0, 0, 0, 0, 0, 0]
	tstep_starting_avg_LA_body = 990
	counter_LA_body = 0

	for i in range(N_tsteps):
		tstep = t_init + i
		vtu_path_LAA =  "/home/sergionio/Downloads/test_nu/200_tsteps/C_laa_{}.vtu".format(str(tstep))
		vtu_path_LA_body =  "/home/sergionio/Downloads/test_nu/200_tsteps/C_la_body_{}.vtu".format(str(tstep))

		avg_LAA, counter_LAA = compute_write_Quartiles(vtu_path_LAA, file_LAA, avg_LAA, tstep_starting_avg_LAA, counter_LAA)
		avg_LA_body, counter_LA_body = compute_write_Quartiles(vtu_path_LA_body, file_LA_body, avg_LA_body, tstep_starting_avg_LA_body, counter_LA_body)

	#----------------------------------------------------------------------------------------
	if counter_LAA != 0:
		file_LAA.write( "avg"+"\t"+str(avg_LAA[0]/counter_LAA)+"\t"+str(avg_LAA[1]/counter_LAA)+"\t"+str(avg_LAA[2]/counter_LAA)+"\t"+str(avg_LAA[3]/counter_LAA)+"\t"+str(avg_LAA[4]/counter_LAA)+"\t"+str(avg_LAA[5]/counter_LAA))
		file_LA_body.write( "avg"+"\t"+str(avg_LA_body[0]/counter_LA_body)+"\t"+str(avg_LA_body[1]/counter_LA_body)+"\t"+str(avg_LA_body[2]/counter_LA_body)+"\t"+str(avg_LA_body[3]/counter_LA_body)+"\t"+str(avg_LA_body[4]/counter_LA_body)+"\t"+str(avg_LA_body[5]/counter_LA_body))
	file_LAA.close()
	file_LA_body.close()