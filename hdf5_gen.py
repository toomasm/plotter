import numpy as np
from collections import defaultdict
import h5py


def get_hdf_data(primes, dataset_names,processing_names,_input_):
    compatibility_naming = ['data1', 'data2', 'data3', 'data4']
    with h5py.File('hdf5/plot_data.hdf', 'r') as h5file:
        data_dic = defaultdict(lambda : defaultdict(dict))
        for i in range(len(_input_)):
            req_data = h5file[primes[i]][dataset_names[i]][processing_names[i]]
            for name, data in req_data.items():
                if 'colour' in name:
                    data_temp = np.array(data).tolist()
                    data = []
                    for x in data_temp:
                        if i == 1:
                            if x == 'b':
                                x = 'g'
                        if i == 2:
                            if x == 'b':
                                x = 'y'
                            if x == 'r':
                                x = 'c'
                        if i == 3:
                            if x == 'b':
                                x = 'g'
                            if x == 'r':
                                x = 'c'
                        data.append(x)
                else:
                    data = np.array(data).tolist()
                data_dic[compatibility_naming[i]][name] = data
        return data_dic


def add_to_hdf(data_dic, hdf_filename):
    with h5py.File(hdf_filename, 'a') as h5_file:
        print 'Moving plot data to .hdf file, named {}'.format(hdf_filename)
        for prime, cond in data_dic.items():
            for condition, process in cond.items():
                for key, req_tuple in process.items():
                    for name in req_tuple._fields:
                        data = getattr(req_tuple, name)
                        if any([type(item) in (float, np.float64) for item in data]):
                            dtype = 'float64'
                        elif all([type(item) == int for item in data]):
                            dtype = 'float64'
                        elif all(type(item) == str for item in data):
                            max_length = max([len(item) for item in data])
                            dtype = '|S{}'.format(max_length)
                        try:
                            h5_key = '{}/{}/{}/{}'.format(prime, condition, key, name)
                            h5_file.create_dataset(h5_key, data=np.array(data, dtype=dtype))
                        except (AttributeError, RuntimeError):
                            continue
