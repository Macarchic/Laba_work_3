import math
import numpy as np
O = 0.7
C = 0.78
H = 0.35
N = 0.71            #these are the covalent radii of the atoms


def connection_checker(first_array, second_array):     #args are arrays with type of atom and coordinates e.g.: ["C", 23.0, 18.9, 14.3]
    
    if first_array[0] == 'O':
        first_radius = 0.7

    elif first_array[0] == 'C':
        first_radius = 0.78

    elif first_array[0] == 'N':
        first_radius = 0.71
    else:
        first_radius = 0.35

    first_x = float(first_array[1])
    first_y = float(first_array[2])
    first_z = float(first_array[3])


    if second_array[0] == 'O':
        second_radius = 0.7
    elif second_array[0] == 'C':
        second_radius = 0.78
    elif second_array[0] == 'N':
        second_radius = 0.71
    else:
        second_radius = 0.35

    second_x = float(second_array[1])
    second_y = float(second_array[2])
    second_z = float(second_array[3])
    d = math.sqrt((first_x-second_x)**2 + (first_y-second_y)**2 + (first_z-second_z)**2)
    if (first_radius + second_radius >= d) and d != 0:
        return 1
    return 0

def nucleoside_building(data):
    with open (data, "r") as file:
        lines = [line.split() for line in file.readlines()]

    incidence_matrix = [[connection_checker(array1, array2) for array1 in lines] for array2 in lines]
    names = [element[0] for element in lines]
    counter_of_O = names.count('O')
    nucleoside = "DG"
    if len(lines) == 29:
        nucleoside = "DC"
    elif len(lines) == 31:
        if counter_of_O == 5:
            nucleoside = "DT"
        else:
            nucleoside = "DA"
    new_names = ["None"]*len(lines)

    #for C1', C2'
    
    for q in range(0, len(lines)):                                      
        if new_names[q] == "None":
            if names[q] == 'C':
                matrix_for_check = [0] * 4  # 0 - C,  1 - H,  2 - O,  3 - N
                for w in range(0, len(lines)):
                    if incidence_matrix[q][w] == 1:
                        if names[w] == 'C':
                            matrix_for_check[0] += 1
                        elif names[w] == 'H':
                            matrix_for_check[1] += 1
                        elif names[w] == 'O':
                            matrix_for_check[2] += 1
                        else:
                            matrix_for_check[3] += 1
    
                if matrix_for_check[0] == 1 and matrix_for_check[1] == 1 and matrix_for_check[2] == 1 and matrix_for_check[3] == 1:
                    new_names[q] = "C1'"
    
                elif matrix_for_check[0] == 2 and matrix_for_check[1] == 2 and matrix_for_check[2] == 0 and matrix_for_check[3] == 0:
                    new_names[q] = "C2'"
    
    
    #for C3', C4', C5''
    
    for _ in range(0, 3):
        for q in range(0, len(lines)):                                 
            if new_names[q] == "None":
                if names[q] == 'C':
                    for w in range(0, len(lines)):
                        if incidence_matrix[q][w] == 1:
                            if new_names[w] == "C2'":
                                new_names[q] = "C3'"
                            elif new_names[w] == "C3'":
                                new_names[q] = "C4'"
                            elif new_names[w] == "C4'":
                                new_names[q] = "C5'"
                            elif new_names[w] == "C2'":
                                new_names[q] = "C3'"
    
    
    
    #for O4', O5', O3'
    
    for q in range(0, len(lines)):
        if new_names[q] == "None":
            if names[q] == 'O':
                for w in range(0, len(lines)):
                    if incidence_matrix[q][w] == 1:
                        if new_names[w] == "C3'":
                            new_names[q] = "O3'"
                        elif new_names[w] == "C5'":
                            new_names[q] = "O5'"
                        elif new_names[w] == "C4'" or new_names[w] == "C1'":
                            new_names[q] = "O4'"
    if nucleoside == "DC" or nucleoside == "DT":
        for q in range(0, len(lines)):                                       #for C2,C4,C5,C6
            if new_names[q] == "None":
                if names[q] == 'C':
                    matrix_for_check = [0] * 4  # 0 - C,  1 - H,  2 - O,  3 - N
                    for w in range(0, len(lines)):
                        if incidence_matrix[q][w] == 1:
                            if names[w] == 'C':
                                matrix_for_check[0] += 1
                            elif names[w] == 'H':
                                matrix_for_check[1] += 1
                            elif names[w] == 'O':
                                matrix_for_check[2] += 1
                            else:
                                matrix_for_check[3] += 1
                    if nucleoside == "DT":  # for C2, C4, C5, C6
                        if matrix_for_check[0] == 0 and matrix_for_check[1] == 0 and matrix_for_check[2] == 1 and matrix_for_check[3] == 2:
                            new_names[q] = "C2"
                        elif matrix_for_check[0] == 1 and matrix_for_check[1] == 0 and matrix_for_check[3] == 1 and matrix_for_check[2] == 1:
                            new_names[q] = "C4"
                        elif matrix_for_check[0] == 3 and matrix_for_check[1] == 0 and matrix_for_check[2] == 0 and matrix_for_check[3] == 0:
                            new_names[q] = "C5"
                        elif matrix_for_check[0] == 1 and matrix_for_check[1] == 1 and matrix_for_check[2] == 0 and matrix_for_check[3] == 1:
                            new_names[q] = "C6"
    
                    elif nucleoside == "DC":  # for C2, C4, C5, C6
                        if matrix_for_check[0] == 0 and matrix_for_check[1] == 0 and matrix_for_check[2] == 1 and matrix_for_check[3] == 2:
                            new_names[q] = "C2"
                        elif matrix_for_check[0] == 1 and matrix_for_check[1] == 0 and matrix_for_check[2] == 0 and matrix_for_check[3] == 2:
                            new_names[q] = "C4"
                        elif matrix_for_check[0] == 2 and matrix_for_check[1] == 1 and matrix_for_check[2] == 0 and matrix_for_check[3] == 0:
                            new_names[q] = "C5"
                        elif matrix_for_check[0] == 1 and matrix_for_check[1] == 1 and matrix_for_check[2] == 0 and matrix_for_check[3] == 1:
                            new_names[q] = "C6"






        #for N1, O2, O4 and last C5''
        for q in range(0, len(lines)):
            if new_names[q] == "None":
                if names[q] == 'N':
                    for w in range(0, len(lines)):
                        if incidence_matrix[q][w] == 1:
                            if new_names[w] == "C1'":
                                new_names[q] = "N1"
                elif names[q] == 'O':
                    for w in range(0, len(lines)):
                        if incidence_matrix[q][w] == 1:
                            if new_names[w] == "C2":
                                new_names[q] = "O2"
                            if new_names[w] == "C4":
                                new_names[q] = "O4"
                elif names[q] == 'C':
                    for w in range(0, len(lines)):
                        if incidence_matrix[q][w] == 1:
                            if new_names[w] == "C5":
                                new_names[q] = "CC5'"
    
        #for N3
        for q in range(0, len(lines)):
            if new_names[q] == "None":
                if names[q] == 'N':
                    for w in range(0, len(lines)):
                        if incidence_matrix[q][w] == 1:
                            if new_names[w] == "C2":
                                new_names[q] = "N3"
    
        #for N4
        for q in range(0, len(lines)):
            if new_names[q] == "None":
                if names[q] == 'N':
                    for w in range(0, len(lines)):
                        if incidence_matrix[q][w] == 1:
                            if new_names[w] == "C4":
                                new_names[q] = "N4"








    if nucleoside == "DA" or nucleoside == "DG":
        #for DG C2, C4, C5, C6, C8                             DA C5
    
        for q in range(0, len(lines)):
            if new_names[q] == "None":
                if names[q] == 'C':
                    matrix_for_check = [0] * 4  # 0 - C,  1 - H,  2 - O,  3 - N
                    for w in range(0, len(lines)):
                        if incidence_matrix[q][w] == 1:
                            if names[w] == 'C':
                                matrix_for_check[0] += 1
                            elif names[w] == 'H':
                                matrix_for_check[1] += 1
                            elif names[w] == 'O':
                                matrix_for_check[2] += 1
                            else:
                                matrix_for_check[3] += 1
    
                    if nucleoside == "DG":  # for C2, C4, C5, C6, C8
                        if matrix_for_check[3] == 3:
                            new_names[q] = "C2"
                        elif matrix_for_check[0] == 1 and matrix_for_check[1] == 0 and matrix_for_check[2] == 0 and matrix_for_check[3] == 2:
                            new_names[q] = "C4"
                        elif matrix_for_check[0] == 2 and matrix_for_check[1] == 0 and matrix_for_check[2] == 0 and matrix_for_check[3] == 1:
                            new_names[q] = "C5"
                        elif matrix_for_check[0] == 1 and matrix_for_check[1] == 0 and matrix_for_check[2] == 1 and matrix_for_check[3] == 1:
                            new_names[q] = "C6"
                        elif matrix_for_check[0] == 0 and matrix_for_check[1] == 1 and matrix_for_check[2] == 0 and matrix_for_check[3] == 2:
                            new_names[q] = "C8"
    
                    elif nucleoside == "DA":                                       #for C5
                        if matrix_for_check[0] == 2 and matrix_for_check[1] == 0 and matrix_for_check[2] == 0 and matrix_for_check[3] == 1:
                            new_names[q] = "C5"
    
    
            for q in range(0, len(lines)):
                if new_names[q] == "None":
                    if names[q] == 'N':
                        for w in range(0, len(lines)):
                            if incidence_matrix[q][w] == 1:
                                if new_names[w] == "C1'":
                                    new_names[q] = "N9"
        #for Adenine
    
        if nucleoside == "DA":
    
    
    
            checker = 0
            #for C2,C4,C6,C8,N1,N3,N6,N7,N9
            for _ in range(0, 8):
                if checker != 8:
                    for q in range(0, len(lines)):
                        if new_names[q] == "None":
                            if names[q] == 'N' or names[q] == 'C':
                                for w in range(0, len(lines)):
                                    if incidence_matrix[q][w] == 1:
                                        if names[q] == 'N':
                                            if new_names[w] == "C5":
                                                new_names[q] = "N7"
                                                checker += 1
    
                                            elif new_names[w] == "C4":
                                                new_names[q] = "N3"
                                                checker += 1
                                            elif new_names[w] == "C2":
                                                new_names[q] = "N1"
                                                checker += 1
                                            elif new_names[w] == "C6":
                                                new_names[q] = "N6"
                                                checker += 1
    
                                        else:
                                            if new_names[w] == "N7":
                                                new_names[q] = "C8"
                                                checker += 1
                                            elif new_names[w] == "N9":
                                                new_names[q] = "C4"
                                                checker += 1
                                            elif new_names[w] == "N3":
                                                new_names[q] = "C2"
                                                checker += 1
                                            elif new_names[w] == "N1":
                                                new_names[q] = "C6"
                                                checker += 1
    
        else:
            # for N1,N3,N7,N9
            for q in range(0, len(lines)):
                if new_names[q] == "None":
                    if names[q] == 'N':
                        for w in range(0, len(lines)):
                            if incidence_matrix[q][w] == 1:
                                if new_names[w] == "C5":
                                    new_names[q] = "N7"
                                elif new_names[w] == "C4":
                                    new_names[q] = "N3"
                                elif new_names[w] == "C6":
                                    new_names[q] = "N1"
    
    
            # for O6, N2
            for q in range(0, len(lines)):
                if new_names[q] == "None":
                    if names[q] == 'N':
                        for w in range(0, len(lines)):
                            if incidence_matrix[q][w] == 1:
                                if new_names[w] == "C2":
                                    new_names[q] = "N2"
                    elif names[q] == 'O':
                        for w in range(0, len(lines)):
                            if incidence_matrix[q][w] == 1:
                                if new_names[w] == "C6":
                                    new_names[q] = "O6"
    for element_newnames in range(0, len(lines)):
        if new_names[element_newnames] == "None":
            new_names[element_newnames] = "H"

    
    with open('results_for_laba1.pdb', 'w') as file:
        counter = 1
        for index in range(0, len(lines)):
            cord_parts = lines[index]
            file.write("ATOM" + " "*4 + str(counter))
            if counter < 10:
                file.write(" "*4)
            elif counter < 100:
                file.write(" "*3)
            elif counter < 1000:
                file.write(" "*2)
            elif counter < 10000:
                file.write(" ")
    
            file.write(new_names[index])
            if len(new_names[index]) == 1:
                file.write(" " * 4)
            elif len(new_names[index]) == 2:
                file.write(" " * 3)
            elif len(new_names[index]) == 3:
                file.write(" " * 2)
            elif len(new_names[index]) == 4:
                file.write(" ")
            file.write(nucleoside + " " * 3 + "1" + " " * 6)
            for mini_index in range(1, 4):
                if float(cord_parts[mini_index]) == 0:
                    file.write("0.0000   ")
                elif float(cord_parts[mini_index]) < 0:
                    file.write("{:.3f}".format(float(cord_parts[mini_index]), 3) + " " * 3)
                else:
                    file.write("{:.4f}".format(float(cord_parts[mini_index]), 3) + " " * 3)
            file.write("1.00" + " " + "0.00" + " "*11 + names[index] + "\n")
            counter += 1
    return new_names, lines



class Vector:
    def __init__(self, array_of_atoms_for_vectors, needed_atoms, special_atom_names_with_coordinates):
        self.array_atom = array_of_atoms_for_vectors
        self.result = []
        indexes = []
        indexes.append(needed_atoms.index(self.array_atom[0]))
        indexes.append(needed_atoms.index(self.array_atom[1]))
        indexes.append(needed_atoms.index(self.array_atom[2]))
        indexes.append(needed_atoms.index(self.array_atom[3]))

        for i in range(1, 4):
            v = [float(special_atom_names_with_coordinates[indexes[i]][1]) - float(special_atom_names_with_coordinates[indexes[i-1]][1]), float(special_atom_names_with_coordinates[indexes[i]][2]) - float(special_atom_names_with_coordinates[indexes[i-1]][2]),
                 float(special_atom_names_with_coordinates[indexes[i]][3]) - float(special_atom_names_with_coordinates[indexes[i-1]][3])]
            self.result.append(v)







def vector_cross(v1,v2):
    res = []
    res.append(v1[1] * v2[2] - v1[2] * v2[1])
    res.append(v1[2] * v2[0] - v1[0] * v2[2])
    res.append(v1[0] * v2[1] - v1[1] * v2[0])
    return res

def dot(A, B):
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2]
def norm(A):
    return (A[0]**2 + A[1]**2 + A[2]**2)**0.5

def angle(u1, u2):
    return np.arccos(dot(u1, u2)/(norm(u1)*norm(u2)))*180/math.pi
def adjusment(u1, u2, v2):
    if np.linalg.det([u1, u2, v2]) < 0:
        return -1
    return 1
def tan(ns):
    return ((ns[4]+ns[1]) - (ns[3] + ns[0])) / (2*ns[2]*(math.sin(36) + math.sin(72)))

def arctan(tan):
    return np.arctan(tan)

def nmaxes(ns, arctan):
    return ns[2]/math.cos(arctan)

def calculate_angle(list, special_atom_names_with_coordinates, needed_atoms ):
    atom1 = list[0]
    atom2 = list[1]
    atom3 = list[2]
    coordinate1 = special_atom_names_with_coordinates[needed_atoms.index(atom1)][1:4]
    coordinate2 = special_atom_names_with_coordinates[needed_atoms.index(atom2)][1:4]
    coordinate3 = special_atom_names_with_coordinates[needed_atoms.index(atom3)][1:4]
    vector1 = [float(coordinate1[0]) - float(coordinate2[0]), float(coordinate1[1]) - float(coordinate2[1]), float(coordinate1[2]) - float(coordinate2[2])]
    vector2 = [float(coordinate3[0]) - float(coordinate2[0]), float(coordinate3[1]) - float(coordinate2[1]), float(coordinate3[2]) - float(coordinate2[2])]

    dot_product = np.dot(vector1, vector2)
    norm1 = np.linalg.norm(vector1)
    norm2 = np.linalg.norm(vector2)

    cos_angle = dot_product / (norm1 * norm2)
    angle_in_radians = math.acos(cos_angle)*180/math.pi
    return angle_in_radians