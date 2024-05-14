import numpy as np
import h5py

class diagnostic:
    def __init__(self,function, variable):
        self.variable = variable
        self.function = function


    def integrate(self):
        integrate = 0
        for i in range(len(self.variable)-1):
            integrate += (self.variable[i+1] - self.variable[i])*(self.function[i+1] + self.function[i]) / 2.

        return integrate


    def differentiate(self):
        diff = np.zeros(len(self.variable))

        diff[0]  = (self.function[1] - self.function[0]) / (self.variable[1] - self.variable[0])
        diff[-1] = (self.function[-1] - self.function[-2]) / (self.variable[-1] - self.variable[-2])

        for i in range(1, len(diff) - 1):
            diff[i] = (self.function[i+1] - self.function[i-1]) / (self.variable[i+1] - self.variable[i-1])

        return diff

    def depth_average_velocity(self):

        da_field = (self.integrate()) / self.variable[-1]

        return da_field



class file_process:

    def __init__(self, file_name, variable, file_type, variable_name = None):
        self.file_name = file_name
        self.variable  = variable
        self.formate = file_type
        self.variable_name = variable_name
    
    def read(self):
        if self.formate=='hdf5':
            path = self.file_name + '/' + str(self.variable)
            filename  = path.split("/")[-1]
            filename  = filename.strip(".h5")
            file_read = h5py.File(path, 'r')
            data      = file_read[filename]
        return data

    def write(self):
        if self.formate=='hdf5':
            data = self.variable
            path = self.file_name + '/' + self.variable_name
            file_name = path.split("/")[-1]
            file_name = file_name.strip(".h5")
            file_write = h5py.File(path, 'w')
            file_write[file_name] = data

            file_write.close()

        elif self.formate == '.nc':
            
            path = self.file_name
            data_set = Dataset(path + '.nc', mode = 'w')
            shape = len(self.variable.shape)-1

            if shape == 1:
                data_set.createDimension('x', len(coord))
            elif shape == 2:
                data_set.createDimension('x', len(coord[0]))
                data_set.createDimension('y', len(coord[1]))
            elif shape == 3:
                data_set.createDimension('x', len(coord[0]))
                data_set.createDimension('y', len(coord[1]))
                data_set.createDimension('z', len(coord[2]))


            if var_type == 'vec':
                data_set.createDimension('vec', 3)
            elif var_type == 'tensor':
                data_set.createDimension('ten', 9)

            if shape == 1:
                if var_type=='vec':
                    data = data_set.createVariable(self.var_name, np.float64, ('vec', 'x'))
                    data[:] = self.variable
                elif var_type == 'tensor':
                    data = data_set.createVariable(self.var_name, np.float64, ('ten', 'x'))
            elif shape == 2:
                if var_type == 'vec':
                    data = data_set.createVariable(self.var_name, np.float64, ('vec', 'x', 'y'))
                    data[:, :, :] = self.variable
                elif var_type == 'tensor':
                    data = data_set.createVariable(self.var_name, float64, ('ten', 'x', 'y'))
                    data[:, :, :] = self.variable
            elif shape == 3:
                if var_type == 'vec':
                    data = data_set.createVariable(self.var_name, float64, ('vec', 'x', 'y', 'z'))
                    data[:, :, :, :] = self.variable
                elif var_type == 'tensor':
                    data = data_set.createVariable(self.var_name, float64, ('ten', 'x', 'y', 'z'))
                    data[:, :, :, :] = self.variable
            data_set.close()





        
def read_hdf5(filename):
    path      = file_name
    data      = path.split("/")[-1]
    data1     = data.strip(".h5")
    file_read = h5py.File(path,'r')
    Pr        = file_read[data1]

    return Pr

def write_hdf5(filename, variable):
    data = variable
    path = filename
    data_name = path.split("/")[-1]
    data1 = data_name.strip(".h5")
    print (data1)
    file_write = h5py.File(path,'w')
    file_write[data1] = data#[:,:,:]
    file_write.close()

