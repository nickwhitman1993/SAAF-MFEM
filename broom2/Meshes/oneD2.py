#!/usr/bin/
# -*- coding: ascii -*-
# this file creates a periodic square mesh that mimics
#     a 1-D problem. It is 2 elements tall with variable
#     length. The grid will go from (0,0) to (2,8) because
#     that is the length of the Reed-Hill problem.


#num_cells = raw_input("Enter the number of spacial cells: ")
num_cells = 16

mesh = open('oneD.mesh', 'w')

mesh.write('MFEM mesh v1.0\n')
mesh.write('\n')
mesh.write('#\n')
mesh.write('# MFEM Geometry Types (see mesh/geom.hpp):\n')
mesh.write('#\n')
mesh.write('# POINT       = 0\n')
mesh.write('# SEGMENT     = 1\n')
mesh.write('# TRIANGLE    = 2\n')
mesh.write('# SQUARE      = 3\n')
mesh.write('# TETRAHEDRON = 4\n')
mesh.write('# CUBE        = 5\n')
mesh.write('#\n')
mesh.write('\n')
mesh.write('dimension\n')
mesh.write('2\n')
mesh.write('\n')
mesh.write('# format: <attribute> <geometry type> <vertex 0> <vertex 1> ...\n')
mesh.write('elements\n')
mesh.write(str(num_cells*2) + '\n')

for i in range(0,num_cells):
   mesh.write(str(i+1) + ' 3 ' + str(i) + ' ' + str(1+i) + ' ' + str(num_cells+2+i) + ' ' + str(num_cells+1+i) + '\n')

for i in range(0,num_cells):
   mesh.write(str(num_cells+i+1) + ' 3 ' + str(num_cells+1+i) + ' ' + str(num_cells+2+i) + ' ' + str(1+i) + ' ' + str(i) + '\n')

mesh.write('\n')
mesh.write('boundary\n')
mesh.write('4\n')
mesh.write('1 1 0 ' + str(num_cells+1) + '\n')
mesh.write('1 1 ' + str(num_cells+1) + ' 0\n')
mesh.write('2 1 ' + str(num_cells) + ' ' + str(2*num_cells+1) + '\n')
mesh.write('2 1 ' + str(2*num_cells+1) + ' ' + str(num_cells) + '\n')
mesh.write('\n')
mesh.write('vertices\n')
mesh.write(str(3*(num_cells+1)) + '\n')
mesh.write('\n')
mesh.write('nodes\n')
mesh.write('FiniteElementSpace\n')
mesh.write('FiniteElementCollection: L2_T1_2D_P1\n')
mesh.write('VDim: 2\n')
mesh.write('Ordering: 1\n')
mesh.write('\n')

for i in range(0,num_cells):
   if(i==0):
      mesh.write(str(i) + ' 0\n')
      mesh.write(str(8./num_cells*(i+1)) + ' 0\n')
      mesh.write(str(i) + ' 1\n')
      mesh.write(str(8./num_cells*(i+1)) + ' 1\n')
      mesh.write('\n')
   else:
      mesh.write(str(8./num_cells*i) + ' 0\n')
      mesh.write(str(8./num_cells*(i+1)) + ' 0\n')
      mesh.write(str(8./num_cells*i) + ' 1\n')
      mesh.write(str(8./num_cells*(i+1)) + ' 1\n')
      mesh.write('\n')

for i in range(0,num_cells):
   if(i==0):
      mesh.write(str(i) + ' 1\n')
      mesh.write(str(8./num_cells*(i+1)) + ' 1\n')
      mesh.write(str(i) + ' 2\n')
      mesh.write(str(8./num_cells*(i+1)) + ' 2\n')
      mesh.write('\n')
   else:
      mesh.write(str(8./num_cells*i) + ' 1\n')
      mesh.write(str(8./num_cells*(i+1)) + ' 1\n')
      mesh.write(str(8./num_cells*i) + ' 2\n')
      mesh.write(str(8./num_cells*(i+1)) + ' 2\n')
      mesh.write('\n')


mesh.close()
