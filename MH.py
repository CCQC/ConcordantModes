import os
rootdir = os.getcwd()
# print(rootdir)
# print(os.path.realpath(__file__))

packagepath = os.path.realpath(__file__)
packagepath = packagepath[:-6]
# print(path)

os.chdir(packagepath + '/subprocess_Scripts')
print(os.getcwd())
os.chdir(rootdir)
print(os.getcwd())

