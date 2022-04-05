import numpy as np
import json
import os
import shutil
import re


class Reap(object):
    def __init__(
        self,
        prog_name,
        zmat,
        disp_cart,
        options,
        n_coord,
        eigs,
        indices,
        energy_regex,
        success_regex,
        deriv_level=0
    ):
        self.prog_name = prog_name
        self.zmat = zmat
        self.disp_cart = disp_cart
        self.options = options
        self.n_coord = n_coord
        self.dispSym = [0]
        # nate
        self.eigs = eigs
        self.energy_regex = energy_regex
        self.success_regex = success_regex
        self.indices = indices
        self.deriv_level = deriv_level

    def run(self):
        # Define energy/gradient search regex
        if not self.deriv_level:
            energy_regex = re.compile(self.energy_regex)
            success_regex = re.compile(self.success_regex)
            self.energies = np.array([])
        else:
            # print(self.energy_regex[0])
            # print(self.energy_regex[1])
            grad_regex1 = re.compile(self.energy_regex[0])
            grad_regex2 = re.compile(self.energy_regex[1])
        eigs = self.eigs
        if type(eigs) == int:
            size = eigs
        else:
            size = len(eigs)
        # n_disp = len(self.disp_cart)

        if not self.deriv_level:
            print(
                "If something looks wrong with the final frequencies, check these energies!"
            )
            print("(Job number 1 == Reference energy) :D")
            if self.options.dir_reap:
                os.chdir("./" + str(1))
                with open("output.dat", "r") as file:
                    data = file.read()
                if not re.search(success_regex, data):
                    print("Energy failed at " + str("ref"))
                    raise RuntimeError
                os.chdir("..")
            else:
                with open("output.1.dat", "r") as file:
                    data = file.read()
                if not re.search(success_regex, data):
                    print("Energy failed at " + str("ref"))
                    raise RuntimeError

            print("this is the energy regex", energy_regex)
            ref_en = float(re.findall(energy_regex, data)[0])
            print("Reference energy: " + str(ref_en))

            indices = self.indices
            # eigs = self.eigs
            # eigs = len(eigs)
            p_en_array = np.zeros((size, size))
            m_en_array = np.zeros((size, size))
            rel_en_p = np.zeros((size, size))
            rel_en_m = np.zeros((size, size))
            relative_energies = []
            absolute_energies = [[("ref", "ref"), "ref", ref_en, 1]]

            direc = 2
            for index in indices:
                i, j = index[0], index[1]
                p_en_array[i, j] = energy = self.reap_energies(
                    direc, success_regex, energy_regex
                )
                print(energy)
                rel = ref_en - energy
                rel_en_p[i, j] = rel
                relative_energies.append([(i, j), "plus", rel, direc])
                absolute_energies.append([(i, j), "plus", energy, direc])
                m_en_array[i, j] = energy = self.reap_energies(
                    direc + 1, success_regex, energy_regex
                )
                print(energy)
                rel = ref_en - energy
                rel_en_m[i, j] = rel
                relative_energies.append([(i, j), "minus", rel, direc + 1])
                absolute_energies.append([(i, j), "minus", energy, direc + 1])
                direc += 2

            self.p_en_array = p_en_array
            self.m_en_array = m_en_array
            self.ref_en = ref_en
            # print_en = np.insert(absolute_energies,0,[("ref", "ref"), "ref", ref_en, 1],axis=0)
            print_en = absolute_energies
            np.set_printoptions(precision=2, linewidth=120)
            print(
                "Relative energies plus-displacements on the diagonal and plus/plus-displacements on the off-diagonal elements"
            )
            print(rel_en_p)
            print(
                "Relative energies minus-displacements on the diagonal and minus/minus-displacements on the off-diagonal elements"
            )
            print(rel_en_m)
            os.chdir("..")
            if self.options.printout_rel_e:
                auxiliary = ""
                header = "Index, relative energy, directory \n"
                print(json.dumps(energy))
                with open("auxiliary", "a") as file:
                    file.seek(0)
                    file.truncate()
                    file.writelines(header)
                for energy in print_en:
                    with open("auxiliary", "a") as file:
                        file.writelines(str(energy) + "\n")
        else:
            indices = self.indices
            p_grad_array = np.array([])
            m_grad_array = np.array([])
            for index in indices:
                grad = self.reap_gradients(2*index+1,grad_regex1,grad_regex2)
                p_grad_array = np.append(p_grad_array,grad,axis=0)
                grad = self.reap_gradients(2*index+2,grad_regex1,grad_regex2)
                m_grad_array = np.append(m_grad_array,grad,axis=0)
            
            self.p_grad_array = p_grad_array.reshape((-1,len(grad)))
            self.m_grad_array = m_grad_array.reshape((-1,len(grad)))
            os.chdir("..")
                

    def reap_energies(self, direc, success_regex, energy_regex):
        if self.options.dir_reap:
            os.chdir("./" + str(direc))
            with open("output.dat", "r") as file:
                data = file.read()
            if not re.search(success_regex, data):
                print("Energy failed at " + os.getcwd())
                raise RuntimeError
            energy = float(re.findall(energy_regex, data)[0])

            os.chdir("..")
        else:
            with open("output." + str(direc) + ".dat", "r") as file:
                data = file.read()
            if not re.search(success_regex, data):
                print("Energy failed at " + os.getcwd())
                raise RuntimeError
            energy = float(re.findall(energy_regex, data)[0])
            
        return energy

    def reap_gradients(self, direc, grad_regex1, grad_regex2):
        os.chdir("./"+str(direc))
        # print(direc)
        grad_array = np.array([])
        with open("output.dat", "r") as file:
            data = file.readlines()
        for i in range(len(data)):
            grad1 = re.search(grad_regex1,data[i])
            if grad1:
                beg_grad = i + 1
                break
        for i in range(len(data) - beg_grad):
            grad2 = re.search(grad_regex2,data[i+beg_grad])
            if grad2:
                end_grad = i + beg_grad
                break
        label_xyz = r"(\s*.*(\s*-?\d+\.\d+){3})+"
        for line in data[beg_grad:end_grad]:
            if re.search(label_xyz,line):
                temp = line.split()[-3:]
                # print(temp)
                grad_array = np.append(grad_array,np.array(temp))
            # regex = grad_regex + label_xyz
            # grad_str = re.search(regex,data).group()
        grad_array = grad_array.astype("float64")
        if not grad1:
            print("Gradient failed at " + os.getcwd())
            raise RuntimeError
        os.chdir("..")
        # raise RuntimeError

        # grad_lines = grad_str.split("\n")
        # twoD_grad_str = [line.split()[-3:] for line in grad_lines[-self.molecule.natom:]]
        # grad_array = np.asarray(twoD_grad_str).astype(float)
        # print(grad_array)
        return grad_array
        # return 0

