#!/usr/bin/python3
"""
2022/09/28
2022/11/08修改：添加了不需要EDIFFG参数的功能
2022/12/01修改：添加了从VASP5格式的POSCAR中读取原子冻结信息，并校正受力的功能
2023/02/19修改：添加了力收敛判断的功能.可从INCAR中读取收敛判据并自动判断
使用力收敛还是能量收敛，无需任何输入参数。
2023/04/25修改：为避免使用力收敛时输出内容过多，修改为使用力收敛时仅输出最后10步的收敛情况
2023/04/27修改：添加了-e和-f选项，用于强制使用能量收敛与力收敛。
               添加了--ediff和--ediffg选项，分别用于指定力收敛和能量收敛时的收敛判据。
               使用力收敛时默认忽略固定的原子。
2023/05/04修改：添加了-en选项，用于指定使用力收敛判断时输出最后多少步的能量变化。默认为最后10步。
                添加了-fn选项，用于指定使用能量收敛判断时输出哪一步的受力。默认为最后一步。
                完善了根据命令行参数使用收敛方式及判据的确定方式。
朱洪
用于监控VASP几何优化是否收敛的程序，只适用于ISIF=2的情形。
用法：vaspoptchk.py
"""
import math
import argparse

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

def check_energy_convergence(ediffg, print_steps):
    """
    This function do energy convergence check by reading and printing
    energy and delta energy of every optimization step from OSZICAR.
    if |Delta_E/conv| < 1, the optimization is converged.
    """
    print("Doing energy convergence check with EDIFFG=%g." % ediffg)
    try:
        oszicar = open("OSZICAR", "r")
    except FileNotFoundError:
        print("OSZICAR not found! Abort.")
        return False, None

    print("%5s%14s%12s%15s" % ("step", "energy", "Delta_E", "|Delta_E/conv|"))
    oszicar_line = oszicar.readline()
    infos = []
    while oszicar_line:
        if "E0" in oszicar_line:
            energy_words = oszicar_line.split()
            step = int(energy_words[0])
            energy = float(energy_words[4])
            delta_energy = float(energy_words[7][1:])
            infos.append({'step': step, 'energy': energy, 'delta_energy': delta_energy, 'thresholds': math.fabs(delta_energy / ediffg)})

        oszicar_line = oszicar.readline()
    
    for info in infos[-print_steps:]:
        print("%5d%14.5f%12.5f%15.5f" % (info['step'], info['energy'], info['delta_energy'], info['thresholds']))
    oszicar.close()

def check_force_convergence(ediffg, step_number):
    """
    This function do force convergence check by reading and printing
    energy and force of every optimization step from OUTCAR.
    if forces acting on every atom < |EDIFFG|, the optimization is converged.
    For large OUTCAR, this function is slow.
    """
    print("Doing force convergence check with EDIFFG=%g." % ediffg)

    #Read forces from OUTCAR
    outcar = open("OUTCAR", "r")
    ion_pos_force = []
    ion_pos_forces = []

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
                ion_pos_forces.append(ion_pos_force)

        line = outcar.readline()
    
    #Select and process force for given step
    if step_number > 0:
        ion_pos_force = ion_pos_forces[step_number-1]
        print("Force for step %d" % step_number)
    elif step_number < 0:
        ion_pos_force = ion_pos_forces[step_number]
        print("Force for step %d" % (len(ion_pos_forces) + step_number + 1))
    else:
        print("step number 0 is not allowed.")
        return
    
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
    
    atom_names = []
    for ind, item in enumerate(poscar_content[5]):
        atom_names += [item]*int(poscar_content[6][ind])

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

    print("  # name      x        y        z        F      opt?     F_corr   converged?")
    for i in range(len(ion_pos_totalforce)):
        print_line = '%3d ' % (i+1)
        print_line += '%4s' % atom_names[i]
        for j in range(4):
            print_line += "%9.5f" % ion_pos_totalforce[i][j]
        print_line += "  %1s  %1s  %1s" % (coords[i][3], coords[i][4], coords[i][5])
        if coords[i][3] == "F" and coords[i][4] == "F" and coords[i][5] == "F":
            print_line += "      0.0"
        else:
            print_line += "%9.5f" % ion_pos_corrected_force[i]
        if ion_pos_corrected_force[i] > math.fabs(ediffg):
            print_line += "   x"
        if not (coords[i][3] == "F" and coords[i][4] == "F" and coords[i][5] == "F"):
            print(print_line)

parser = argparse.ArgumentParser()
group1 = parser.add_mutually_exclusive_group()
group1.add_argument("-f", help="use force convergence check", action="store_true")
group1.add_argument("-e", help="use energy convergence check", action="store_true")
group2 = parser.add_mutually_exclusive_group()
group2.add_argument("--energy-ediffg", help="ascribe criteria used by energy convergence check",
                    type=float, default=1e-3)
group2.add_argument("--force-ediffg", help="ascribe criteria used by force convergence check", 
                    type=float, default=-0.02)
parser.add_argument("-en", help="the number of steps to print while using energy convergence check", 
                    type=int, default=10)
parser.add_argument("-fn", help="""the number of step used to print whilt using force convergence check.\n
                                Positive value means usual frame count, and negative value means count from\n
                                last frame to first frame. 0 is not allowed.""", 
                    type=int, default=-1)
args = parser.parse_args()

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

print("INCAR: EDIFF=%g, EDIFFG=%g" % (ediff, ediffg))
print("Script provided EDIFFG args for energy check: %g, for force check: %g" % (args.energy_ediffg, args.force_ediffg))

if args.e:
    if ediff_read_flag is True and ediffg_read_flag is True:
        if ediffg <= 0:
            print("EDIFFG < 0 detected. Using %g from scirpt args instead." % args.energy_ediffg)
            ediffg = args.energy_ediffg
    elif ediff_read_flag is True and ediffg_read_flag is False:
        ediffg = ediff * 10
    elif ediff_read_flag is False and ediffg_read_flag is True:
        if ediffg <= 0:
            print("EDIFFG < 0 detected. Using %g from scirpt args instead." % args.energy_ediffg)
            ediffg = args.energy_ediffg
    else:
        ediffg = ediffg_default
    check_energy_convergence(ediffg, args.en)

elif args.f:
    if ediffg_read_flag is True:
        if ediffg >= 0:
            print("EDIFFG > 0 detected. Using %g from script args instead." % args.force_ediffg)
            ediffg = args.force_ediffg
    else:
        ediffg = args.force_ediffg
    check_force_convergence(ediffg, args.fn)
else:
    if ediffg > 0:
        check_energy_convergence(ediffg, args.en)
    else:
        check_force_convergence(ediffg, args.fn)
