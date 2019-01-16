"""
This module contains the BlockFileMap object and a function for initializing blockFile directory structure
"""
import numpy as np
import sys
import os

def init_files(base_directory, nb):
    """ Creates the desired directory structure and files """
    for bi in range(nb):
        # Make folders
        for x in ["_d", "_p", "_n"]:
            folder = base_directory + "/" + str(bi) + x
            if not os.path.exists(folder):
                os.makedirs(folder)
            for bj in range(bi, nb):
                fp = folder + "/" + str(bj) + x + ".block"
                f = open(fp, "w")
                f.close()

class BlockFileMap():
    """ BlockFile object for reading and writing data to block matrices on disk"""
    def __init__(self, path, dtype, block_shape):
        """ Initialized with data, creates a file, writes to the file """
        self.file_path = path
#        self.block_size = block_size
        self.block_shape = block_shape
        self.dtype = dtype
        self.mmap = None

    def __getitem__(self, item):
        """ Getter to directly access mmap """
        return self.mmap[item]

    def __setitem__(self, item, val):
        """ Setter to directly access mmap """
        self.mmap[item] = val

    def _mmap(self):
        """ Create the memory map object """
        self.mmap = np.memmap(self.file_path, dtype=self.dtype, mode='r+', shape=self.block_shape, order='C')
#        self.mmap = np.memmap(self.file_path, dtype=self.dtype, mode='r+', shape=(self.block_size, self.block_size), order='C')

    def _umap(self):
        """ Flush changes to file, remove """
        self.mmap.flush()
        del self.mmap

    def _read_all(self):
        """ Returns whole block data as np.array"""
        return self.mmap

    def _write_all(self, x):
        """ Writes whole block as X """
        self.mmap[:] = x

    def open(self):
        """ Opens a block file to permit read/write """
        self._mmap()

    def close(self):
        """ Closes a block file """
        self._umap()

    def print(self):
        """ Print file for debug """
        print(self.mmap)

    def read_all(self):
        """ Reads a whole block, returns 2d numpy.array """
        return self._read_all()
    
    def write_all(self, arr):
        """ Writes whole 2d numpy array arr into block file """
        self._write_all(arr)

    def init_file(self, arr):
        """ Writes a file normally, since we cannot mmap an empty file """
        fh = open(self.file_path, "wb")
        fh.write(arr.tobytes(order='C'))
        fh.close()

if __name__ == "__main__":
    make_folders("test", 5)
    sys.exit()

    print("/********** Small map test *********/")
    # Test
    print("Testing blocks with small printable array")
    print("Initializing array")
    arr = np.array([[1,0],[0,1]], dtype=np.uint8)
    print(arr)
    print(arr.shape)
    print(type(arr))
    print()
    testblock = BlockFileMap("test_small.block", arr, np.uint8)
    testblock.open()
    print()
    print("Result:")
    testblock.print()
    testblock.close()
    print()

    testblock.open()
    print("Overwriting block with:")
    arr2 = np.array([[9,1,3],[2,14,18],[14,17,88]], dtype=np.int32)
    print(arr2)
    print()
    print("Result:")
    testblock.write_all(arr2)
    testblock.print()
    testblock.close()

    testblock.open()
    testblock.print()
    testblock.close()
    print()

    testblock.open()
    print("Reading block, sorting by column 0, and overwriting")
    arr2 = testblock.read_all()
    ind = np.argsort(arr2[:,0])
    arr2 = arr2[ind]
    print(arr2)
    print()
    print("Result:")
    testblock.write_all(arr2)
    testblock.close()
    testblock.open()
    testblock.print()
    testblock.close()
    print()

    print("Editing elements at:")
    editi = np.array([[0,0],[2,2],[1,1]])
#    print(editi)
    print("[:,0]")
    print("with:")
    editv = np.array([7, 5, 7], dtype=np.int64)
    print(editv)
    print()
    testblock.open()
    testblock.mmap[:,0] = editv
    testblock.close()

    testblock.open()
    print("Result:")
    testblock.read_all()
    testblock.print()
    testblock.close()
    print()

    print("Reading elements at")
    print("[:2,1:]")
    testblock.open()
    a = [i for i in testblock[:2,1:]]
    print(a)
    testblock.close()

    print()
    print("Editing elements with a set of indices:")
    testblock.open()
    testblock[:,np.array([0,2])] = [999,999]
    testblock.close()
    testblock.open()
    testblock.print()
    testblock.close()
    




