import numpy as np
import subprocess as sp
import re
import os
import shutil

# Define Regexes
fcRegex = re.compile('\s+\d+\s+\d+\s+(\d+\.\d+)')
f16Regex = re.compile('\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')

# Open some files and read in the pertinent values
with open('fc','r') as file:
    fcdata = file.readlines()

with open('FC_f16_init','r') as file:
    f16 = file.readlines()

f16_header = f16[0]

# generate list of diagonal force constants
fc = np.array([])
for i in fcdata:
    fc = np.append(fc,re.findall(fcRegex,i))

# generate list of lower level diagonal + off-diagonal force constants
f16_arr = np.array([])
for i in f16:
    f16_arr = np.append(f16_arr,re.findall(f16Regex,i))

dim = int(np.sqrt(len(f16_arr)))

f16_mat = np.reshape(f16_arr,(dim,dim))


# This line of code may have to be checked. I'm only using it because it seems that the force constants are closer to the DZ when flipped.
# This can be checked for my own structures.
fc = np.flip(fc,0)

f16_new = f16_mat
for i in range(dim):
    # This is a good debugging check, force constants should be similar
    # print(f16_new[i][i])
    # print(fc[i])
    f16_new[i][i] = fc[i]

f16_flat = f16_mat.flatten()

f16_outputString = f16_header
for i in range(len(f16_flat)):
    if (i+1)%3 == 0:
        f16_outputString += "{:20.10f}".format(float(f16_flat[i])) + '\n'
    else:
        f16_outputString += "{:20.10f}".format(float(f16_flat[i]))

f16_outputString = f16_outputString[:-1]

with open('input_template','r') as file:
    inputString = file.readlines()

# print(inputString[6])
# print(inputString[6][34])
inputString[6] =inputString[6][:34] + str(1) + inputString[6][35:]
# print(inputString[6])

with open('input_template','w') as file:
    file.writelines(inputString)


# This code will be for generating all the necessary intder files for the last frequency computation

if os.path.exists('./IntderFinal'):
    shutil.rmtree('./IntderFinal')

os.mkdir('IntderFinal')
os.chdir('IntderFinal')
with open('file16','w') as file:
    file.write(f16_outputString)
shutil.copy('../input_template','intder.inp')
os.system('/home/vulcan/mel64643/bin/INTDER < intder.inp > intder.out')
# sp.call(['/home/vulcan/mel64643/bin/INTDER','< intder.inp','> intder.out'])

print(os.getcwd())
os.chdir('..')
