#!/usr/bin/python3
"""
2022/09/28
2022/11/08修改：添加了不需要EDIFFG参数的功能
2022/12/01修改：添加了从VASP5格式的POSCAR中读取原子冻结信息，并校正受力的功能
2023/02/19修改：添加了力收敛判断的功能.可从INCAR中读取收敛判据并自动判断
使用力收敛还是能量收敛，无需任何输入参数。
朱洪
用于监控VASP几何优化是否收敛的程序，只适用于ISIF=2的情形。
用法：vaspopt.py
"""
import math

def remove_line_comment(line, comment_chars):
    """
    Get uncommented line from commented or umcommented line
    with ascribed comment chars
    """
    comment_index = len(line)
    for index, char in enumerate(line):
        if char in comment_chars:
            if index < comment_index:
                comment_index = index

    return line[:comment_index]

def read_vasp_keyword(keyword):
    """
    Read input for given VASP keyword from VASP INCAR. If the keyword
    occurs more than once in INCAR, it will return values in the
    first occurence. If the keyword is not specified in INCAR,
    this function will return None.
    This function can only read VASP keyword with only one floating
    number parameter now (such as EDIFF and EDIFFG).
    """
    try:
        incar = open("INCAR", "r")
    except FileNotFoundError:
        print("INCAR not found! Abort.")
        return False, None

    incar_line = incar.readline()
    while incar_line:
        incar_line_uncommented = remove_line_comment(incar_line, ['#', "("])
        if keyword in incar_line_uncommented:
            if incar_line_uncommented.count('=') != 1:
                print("Errorious VASP INCAR for keyword %s!" % keyword)
                incar.close()
                return False, None
            else:
                start_index = incar_line_uncommented.index("=") + 1
                input_args = incar_line_uncommented[start_index:].split()
                for index, item in enumerate(input_args):
                    input_args[index] = float(item)
                incar.close()
                return True, input_args[0]

        incar_line = incar.readline()

    return False, None

def check_energy_convergence(ediffg):
    """
    This function do energy convergence check by reading and printing
    energy and delta energy of every optimization step from OSZICAR.
    if |Delta_E/conv| < 1, the optimization is converged.
    """
    print("Doing energy convergence check.")
    try:
        oszicar = open("OSZICAR", "r")
    except FileNotFoundError:
        print("INCAR not found! Abort.")
        return False, None

    print("%5s%14s%12s%15s" % ("step", "energy", "Delta_E", "|Delta_E/conv|"))
    oszicar_line = oszicar.readline()
    while oszicar_line:
        if "E0" in oszicar_line:
            energy_words = oszicar_line.split()
            step = int(energy_words[0])
            energy = float(energy_words[4])
            delta_energy = float(energy_words[7][1:])
            print("%5d%14.5f%12.5f%15.5f" %
                    (step, energy, delta_energy, math.fabs(delta_energy / ediffg)))

        oszicar_line = oszicar.readline()

    oszicar.close()

def check_force_convergence(ediffg):
    """
    This function do force convergence check by reading and printing
    energy and force of every optimization step from OUTCAR.
    if forces acting on every atom < |EDIFFG|, the optimization is converged.
    For large OUTCAR, this function is slow.
    """
    print("Doing force convergence check.")

    outcar = open("OUTCAR", "r")
    ion_pos_force = []

    line = outcar.readline()
    while line:
        words = line.split()
        if len(words) > 0:
            if words[0] == "POSITION":
                ion_pos_force = []
                line = outcar.readline()
                words2 = line.split()
                while words2[0] != "total":
                    line = outcar.readline()
                    words2 = line.split()
                    ion_pos_force.append(words2)

                if len(ion_pos_force) != 0:
                    ion_pos_force.pop()
                    ion_pos_force.pop()

        line = outcar.readline()

    for i in range(len(ion_pos_force)):
        for j in range(len(ion_pos_force[0])):
            ion_pos_force[i][j] = float(ion_pos_force[i][j])

    ion_pos_totalforce = []
    for i in range(len(ion_pos_force)):
        forcex = ion_pos_force[i][3]
        forcey = ion_pos_force[i][4]
        forcez = ion_pos_force[i][5]
        force = math.sqrt(forcex**2 + forcey**2 + forcez**2)
        pos_totalforce = ion_pos_force[i][:3] + [force]
        ion_pos_totalforce.append(pos_totalforce)

    poscar_content = []
    with open("POSCAR", "r") as fin:
        for lines in fin:
            words = lines.split()
            poscar_content.append(words)

    atom_number = 0
    for item in poscar_content[6]:
        atom_number += int(item)

    if poscar_content[7][0][0] == "S" or poscar_content[7][0][0] == "s":
        coords = poscar_content[9:9+atom_number]

        ion_pos_corrected_force = []
        for i in range(len(ion_pos_force)):
            forcex = 0
            forcey = 0
            forcez = 0
            if coords[i][3] == "T":
                forcex = ion_pos_force[i][3]
            if coords[i][4] == "T":
                forcey = ion_pos_force[i][4]
            if coords[i][5] == "T":
                forcez = ion_pos_force[i][5]
            force = math.sqrt(forcex**2 + forcey**2 + forcez**2)
            ion_pos_corrected_force.append(force)
    else:
        coords = poscar_content[8:8+atom_number]
        ion_pos_corrected_force = []
        for i in range(len(coords)):
            coords[i] += ["T", "T", "T"]
            ion_pos_corrected_force.append(ion_pos_totalforce[i][3])

    print("  #      x        y        z        F      opt?     F_corr   converged?")
    for i in range(len(ion_pos_totalforce)):
        print_line = '%3d ' % (i+1)
        for j in range(4):
            print_line += "%9.5f" % ion_pos_totalforce[i][j]
        print_line += "  %1s  %1s  %1s" % (coords[i][3], coords[i][4], coords[i][5])
        if coords[i][3] == "F" and coords[i][4] == "F" and coords[i][5] == "F":
            print_line += "      0.0"
        else:
            print_line += "%9.5f" % ion_pos_corrected_force[i]
        if ion_pos_corrected_force[i] > math.fabs(ediffg):
            print_line += "   x"
        print(print_line)

ediff_read_flag, ediffg_read_flag = False, False
ediff_default = 1.0e-4
ediffg_default = 1.0e-3
ediff_read_flag, ediff = read_vasp_keyword("EDIFF")
ediffg_read_flag, ediffg = read_vasp_keyword("EDIFFG")

if ediffg_read_flag is False:
    if ediff_read_flag is False:
        ediff = ediff_default
        ediffg = ediffg_default
    else:
        ediffg = ediff * 10
else:
    if ediff_read_flag is False:
        ediff = ediff_default

print("EDIFF=%g, EDIFFG=%g" % (ediff, ediffg))

if ediffg > 0:
    check_energy_convergence(ediffg)
else:
    check_force_convergence(ediffg)
