import shutil
import os
rootdir = os.getcwd()
# print(rootdir)
# print(os.path.realpath(__file__))

packagepath = os.path.realpath(__file__)
packagepath = packagepath[:-6]

# os.chdir(packagepath + '/subprocess_Scripts')
os.chdir(packagepath)
print('Package path:')
print(os.getcwd())
os.chdir(rootdir)
print('Working path:')
print(os.getcwd())

if os.path.exists(rootdir + '/subprocess_Scripts'):
    shutil.rmtree(rootdir + '/subprocess_Scripts')
if os.path.exists(rootdir + '/templates'):
    shutil.rmtree(rootdir + '/templates')



shutil.copytree(packagepath + '/subprocess_Scripts',rootdir + '/subprocess_Scripts')
shutil.copytree(packagepath + '/templates',rootdir + '/templates')






# shutil.rmtree(rootdir + '/subprocess_Scripts')
# shutil.rmtree(rootdir + '/templates')
