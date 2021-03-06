#!/usr/bin/
# -*- coding: ascii -*-
# this file creates a periodic quad mesh that mimics
#     a 1-D problem. It is num_y_cells elements tall and
#     num_x_cells long. The grid will go from x=[0,x_length]
#     and y=[0,y_length]. The Reed-Hill problem is defined on 
#     x=[0,8].
#


num_x_cells = 10.
num_y_cells = 10.
x_length = 1. # problem goes from 0 to x_length
y_length = 1. # problem height from 0 to y_length

mesh = open('WangRagusaUnitPeriodic.mesh', 'w')

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
mesh.write(str(int(num_x_cells*num_y_cells)) + '\n')

for j in range(0,int(num_y_cells)):
   for i in range(0,int(num_x_cells)):
      if(j!=num_y_cells-1 and i!=num_x_cells-1):
         mesh.write(str(int(num_x_cells*j+i+1)) + ' 3 ' + str(int(num_x_cells*j+i)) + ' ' + str(int(num_x_cells*j+i+1)) + ' ' + str(int(num_x_cells*(j+1)+i+1)) + ' ' + str(int(num_x_cells*(j+1)+i)) + '\n')
      elif(j==num_y_cells-1 and i==num_x_cells-1):
         mesh.write(str(int(num_x_cells*j+i+1)) + ' 3 ' + str(int(num_x_cells*j+i)) + ' ' + str(int(num_x_cells*j+i+1-num_x_cells)) + ' ' + str(int(0)) + ' ' + str(int(i)) + '\n')
      elif(i==num_x_cells-1):
         mesh.write(str(int(num_x_cells*j+i+1)) + ' 3 ' + str(int(num_x_cells*j+i)) + ' ' + str(int(num_x_cells*j+i+1-num_x_cells)) + ' ' + str(int(num_x_cells*j+i+1)) + ' ' + str(int(num_x_cells*(j+1)+i)) + '\n')
      elif(j==num_y_cells-1):
         mesh.write(str(int(num_x_cells*j+i+1)) + ' 3 ' + str(int(num_x_cells*j+i)) + ' ' + str(int(num_x_cells*j+i+1)) + ' ' + str(int(i+1)) + ' ' + str(int(i)) + '\n')

mesh.write('\n')
mesh.write('boundary\n')
mesh.write('0\n')
   
mesh.write('\n')
mesh.write('vertices\n')
mesh.write(str(int((num_x_cells)*(num_y_cells))) + '\n')
mesh.write('\n')
mesh.write('nodes\n')
mesh.write('FiniteElementSpace\n')
mesh.write('FiniteElementCollection: L2_T1_2D_P1\n')
mesh.write('VDim: 2\n')
mesh.write('Ordering: 1\n')
mesh.write('\n')

for j in range(0,int(num_y_cells)):
   for i in range(0,int(num_x_cells)):
         mesh.write(str(x_length/num_x_cells*i) + ' ' + str(y_length/num_y_cells*j) + '\n')
         mesh.write(str(x_length/num_x_cells*(i+1)) + ' ' + str(y_length/num_y_cells*j) + '\n')
         mesh.write(str(x_length/num_x_cells*i) + ' ' + str(y_length/num_y_cells*(j+1)) + '\n')
         mesh.write(str(x_length/num_x_cells*(i+1)) + ' ' + str(y_length/num_y_cells*(j+1)) + '\n')
         mesh.write('\n')

mesh.close()
