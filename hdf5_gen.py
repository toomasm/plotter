import numpy as np
from collections import namedtuple
import h5py
from os import path

plot_data = namedtuple('PlottingData', ['nucl_data_16S', 'nucl_data_23S' ,'y_pos_16S', 'y_neg_16S', 'y_pos_23S',
                                        'y_neg_23S', 'colour_16S', 'colour_23S', 'symbol'])

def make_hdf_filename(_input_, scatter_flag, MA_flag, three_prime_flag):
    hdf_filename = '.hdf5'
    for n in (_input_):
        base_name = path.basename(n)
        name = path.splitext(base_name)
        name = name[0]
        hdf_filename = name + hdf_filename
    if scatter_flag:
        hdf_filename = 'Scat_' + hdf_filename
    if MA_flag:
        hdf_filename = 'MA_' + hdf_filename
    if three_prime_flag:
        hdf_filename = '3-prim_' + hdf_filename
    return hdf_filename

def get_hdf_data(hdf_filename):
    h5file = h5py.File(hdf_filename, 'r')
    data_dic = {}
    for key, value in h5file.items():
        for name, data in value.items():
            setattr(plot_data, name, np.array(data).tolist())
        data_dic[key] = plot_data
    h5file.close()
    return data_dic


def add_to_hdf(data_dic, hdf_filename):
    h5_file = h5py.File(hdf_filename, 'w')
    for key, value in data_dic.items():
        for name in value._fields:
            data = getattr(value, name)
            if any([type(item) in (float, np.float64) for item in data]):
                dtype = 'float64'
            elif all([type(item) == int for item in data]):
                dtype = 'float64'
            elif all(type(item) == str for item in data):
                max_length = max([len(item) for item in data])
                dtype = '|S{}'.format(max_length)
            h5_key = '{}/{}'.format(key, name)
            print 'Moving plot data to .hdf file, named {}'.format(hdf_filename)
            h5_file.create_dataset(h5_key, data=np.array(data, dtype=dtype))
    h5_file.close()