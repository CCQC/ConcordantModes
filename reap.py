import numpy as np
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

    def run(self):
        # Define energy search regex
        energy_regex = re.compile(self.energy_regex)
        success_regex = re.compile(self.success_regex)
        eigs = self.eigs
        if type(eigs) == int:
            size = eigs
        else:
            size = len(eigs)
        n_disp = len(self.disp_cart)
        self.energies = np.array([])

        print(
            "If something looks wrong with the final frequencies, check \
                these energies!"
        )
        print("(Job number 1 == Reference energy) :D")

        os.chdir("./" + str(1))
        with open("output.dat", "r") as file:
            data = file.read()
        if not re.search(success_regex, data):
            print("Energy failed at " + str("ref"))
            raise RuntimeError
        os.chdir("..")
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
            m_en_array[i, j] = energy = self.reap_energies(
                direc + 1, success_regex, energy_regex
            )
            print(energy)
            rel = ref_en - energy
            rel_en_m[i, j] = rel
            relative_energies.append([(i, j), "minus", rel, direc + 1])
            direc += 2

        self.p_en_array = p_en_array
        self.m_en_array = m_en_array
        self.ref_en = ref_en
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
            with open("auxiliary", "a") as file:
                file.writelines(header)
            for energy in relative_energies:
                with open("auxiliary", "a") as file:
                    file.writelines(str(energy) + "\n")

    def reap_energies(self, direc, success_regex, energy_regex):
        os.chdir("./" + str(direc))
        with open("output.dat", "r") as file:
            data = file.read()
        if not re.search(success_regex, data):
            print("Energy failed at " + os.getcwd())
            raise RuntimeError
        energy = float(re.findall(energy_regex, data)[0])

        os.chdir("..")
        return energy
        # with open("auxillary.dat", 'w') as file:
        #    for energy in relative_energies:
        #        file.writelines(energy)
        # print(relative_energies)
        # np.set_printoptions(precision=8)


#        self.energiesDict['ref'] = float(self.energies[0])
#        Sum = 0
#        for i in range(n_disp - np.sum(self.dispSym) - 1):
#            j = Sum % 2
#            k = Sum // 2
#            if j == 0:
#                self.energiesDict[str(k+1)+'_plus'] = float(self.energies[i+1])
#            if j == 1:
#                self.energiesDict[str(k+1)+'_minus'] = float(self.energies[i+1])
#            # if self.dispSym[k] == 0:
#                # if j == 0:
#                    # self.energiesDict[str(k+1)+'_plus'] = float(self.energies[i+1])
#                # if j == 1:
#                    # self.energiesDict[str(k+1)+'_minus'] = float(self.energies[i+1])
#            # elif self.dispSym[k] == 1:
#                # self.energiesDict[str(k+1)+'_plus'] = float(self.energies[i+1])
#                # self.energiesDict[str(k+1)+'_minus'] = float(self.energies[i+1])
##            Sum += 1
#            # Sum += self.dispSym[k]
